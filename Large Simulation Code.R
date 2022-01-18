################# Simulation Code #########################

#### Written by: Nicholas Weaver 
#### 
#### Date: December 15, 2021
####
#### Purpose: Function to simulate microbiome data given a set of inputs. This code will perturb large 
####            abundance OTUs at timepoint 1 towards a mean relative abundance at timepoint m with 
####            user-supplied steps over time.
####
#### Inputs: Real_Data = Data set of real OTU counts used to simulate initial time point
####         otuLength = Number of OTUs to simulate, if left blank, it will default to length of Real_Data
####         time_length = Number of time points to simulate
####         nsub = Number of subjects to generate data for at each time point. Default is 60
####         lambda = Threshold to separate 'large' OTUs from 'small' OTUs by relative abundance
####         eta = Proportion of large OTUs that will change in this simulation. Defalut is 0.5
####         mean = Default is "". This is the mean number of reads per observation.
####         size = Default is "". This is the dispersion parameter for the number of reads?
####         non_zero_small_prob = Probability that an OTU that was zero at a previous time point will be non_zero at the next time point
####         non_zero_large_prob = Probability that an OTU that was zero at a previous time point will be non_zero at the next time point
####         prop_cover = Vector of length time_length - 1 that says what proportion of the difference (end-start) should occur at that time point
####         bd_metric = Metric/measure to use for beta diversity calculations
####
####
#### Outputs: Data = A matrix of OTU counts at the subject by time point level
####          Perturbed_OTUs = A list of OTUs that were considered to be large abundance and change towards a mean over time
####          Large_OTUs = A list of OTUs that were considered to be large abundance and stay constant over time
####          Small_OTUs = A list of OTUs that were considered to be small abundance and stay constant over time
####          Beta_Div = The average between subject and within subject beta diversity values at each time point
####          Individual_BD = A Matrix of within subject beta diversity values for all consecutive time points.


#############################################################
library(cluster)
library(spam)
library(ggplot2)
library(tidyr)
library(dplyr)
library(MCMCpack)
library(MiSPU)
library(vegan)
library(truncnorm)
library(RColorBrewer)

SimOTUCounts <- function(Real_Data, hyp_vals, otuLength = "", time_length, nsub = 60, 
                         lambda = 0.01, eta = 0.5, mean = "", size = "", non_zero_small_prob = 0.1,
                         non_zero_large_prob = 0.9,
                         prop_cover = "", bd_metric = "bray"){
  
  if(sum(prop_cover > 1 | prop_cover < 0) != 0){
    stop("Values less than 0 or greater than 1 may lead to negative proportions or proportions greater than 1 (which are impossible).")
  }
  if(mean == ""){
    mean <- round(mean(apply(Real_Data,1,sum)))
  }
  
  if(size == ""){
    #size <- (var(apply(Real_Data,1,sum)) - mean)/(mean^2) #Based on the alternative construction of the negative binomial distribution (see ?rnbinom)
    size = 20
  }
  
  if(otuLength == ""){
    otuLength <- length(hyp_vals) #This will be the number of OTU Counts that we can simulate
  }
  
  #hyp_vals is a vector whose elements are the hyperparameters of the dirichlet distribution
  
  
  #Generate the proportions for each OTU for each subject (CAN THIS BE DONE WITHOUT A LOOP?):
  eta_props <- c()
  for (i in 1:nsub){
    eta_props <- rbind(eta_props, MiSPU::rdirichlet(1, alpha = hyp_vals))
  }
  eta_props <- as.data.frame(eta_props)
  names(eta_props) <- names(hyp_vals)
  
  OTU_counts <- c()
  #Generate the OTU count for each subject and OTU:
  for (i in 1:nsub){
    OTU_counts <- rbind(OTU_counts, cbind(i, 1, t(rmultinom(1, size = rnbinom(1,mu = mean, size = size), prob = eta_props[i,]))))
  }
  OTU_counts <- as.data.frame(OTU_counts)
  
  names(OTU_counts) <- c("Subject", "Time", names(hyp_vals))
  
  #OTU_counts is a data frame that contains the first time point of our simulated data set.
  
  ##Record the average beta diversity between subjects to help track this mean change
  avg_beta_div <- data.frame(Time = NA, Avg_between_BD = NA, Avg_within_BD= NA)
  avg_beta_div <- rbind(avg_beta_div, cbind(Time = 1, Avg_between_BD= mean(vegdist(OTU_counts[,names(hyp_vals)], method = "bray")),
                                            Avg_within_BD = 0))
  
  
  
  #### Establish how microbiome will change ####
  #proportions below/above lambda will change:
  lambda <- lambda #threshold for small/large
  
  #Find starting mean relative abundance vector to identify otus into above/below threshold:
  rel_ab_start <- apply(eta_props, 2, mean)
  names(rel_ab_start) <- names(hyp_vals)
  
  small_otus <- names(rel_ab_start)[which(rel_ab_start < lambda)]
  large_otus <- names(rel_ab_start)[which(rel_ab_start >= lambda)]
  
  #Within the selected group for change, what proportion of OTUs should actually change:
  eta <- eta
  
  #Select the OTU to go to the mean and those to stay constant:
  mean_otus <- sample(large_otus, size = ceiling(eta*length(large_otus)))
  constant_otus <- names(hyp_vals)[!(names(hyp_vals) %in% mean_otus)]
  
  #Keep proportion of composition represented by changing group constant.
  # Relative abundance to be manipulated by large OTUs will always stay the same:
  rel_ab_available <- sum(rel_ab_start[mean_otus])
  
  
  #Establish goal proportions. This will determine the end proportions, on average, for all OTUs.
  # Future work should allow for different ways to identify the end proportions:
  
  goal_props <- hyp_vals #create a new object so we don't ruin the hyp_vals object
  names(goal_props) <- names(hyp_vals) #make sure the list is named
  
  goal_props[mean_otus] <- rel_ab_available/length(mean_otus)  #OTUs that will change will all go to the same proportion
  goal_props[constant_otus] <- rel_ab_start[constant_otus]  #These OTUs won't change from first proportion
  
  
  #Establish monotonic steps to go from start to end (linear if constant steps, exponential, logarithmic, etc)
  if(prop_cover == ""){
    prop_cover <- 1:(time_length-1)/(time_length-1) #Default is constant steps
  }
  
  #Determine how much each vector needs to adjust in total
  # to get to goal props at last time point:
  props_change <- (t(goal_props - t(eta_props)))
  
  
  #Keep zeros within a certain bound:
  center_zero <- round(mean(apply(OTU_counts,1,function(x) sum(x == 0))))
  spread_zero <- 4*sd(apply(OTU_counts,1,function(x) sum(x == 0)))
  
  low_limit <- round(center_zero - spread_zero)
  up_limit <- round(center_zero + spread_zero)
  
  
  ### Now we can generate the next time points! ###
  
  
  #Loop through time points, perturbing previous time point to get new time point values:
  
  individual_bd_matrix <- data.frame(Subject = NA, Time_Comp = NA, Beta_Diversity = NA) #Store beta diversity changes
  
  for(i in 2:time_length){
    within_bd_tracker <- c() #track beta diversity within subjects at all current and past generated times
    for(j in 1:nsub){
      
      #New proportions based upon desired change. Start at 
      eta_props_new <-eta_props[j,] + props_change[j,]*prop_cover[i-1]
      
      if(sum(eta_props_new) != 1){
        warning(paste("Relative Abundance vector for subject ", j, " at time ", i," does not sum to 1. It sums to ",
                      sum(eta_props_new),sep=""))
      }
      
      #Simulate counts for timepoint i
      tmp_vec <- c(j, i, rmultinom(1, size = rnbinom(1,mu = mean, size = size), prob = eta_props_new))
      names(tmp_vec) <- c("Subject", "Time", all_of(names(hyp_vals)))
      

      #How many otus were zero in previous time point?
      n_zeros <- (sum(OTU_counts[OTU_counts$Subject == j & OTU_counts$Time == i-1,names(hyp_vals)] == 0)) 
      
      if(n_zeros > (center_zero + spread_zero)){
        warning(paste("More than ", center_zero + spread_zero," OTUs were simulated to be zero at time point ", i-1,". Change the cutoff value
        to fix this warning",sep = ""))
      }
      
      if(n_zeros < (center_zero - spread_zero)){
        warning(paste("Less than ", center_zero - spread_zero," OTUs were simulated to be zero at time point ", i -1,". Change the cutoff value
                        to fix this warning.",sep = ""))
      }
      
      #How many zeros should be in current time point?
      n_zeros_cur <- round(qtruncnorm(runif(1,0,1),a = low_limit, b = up_limit, n_zeros, spread_zero/4))
      
      #Determine which OTUs will stay zero from the previous time point and adjust appropriately:
      prev_zero <- names(hyp_vals)[which(OTU_counts[OTU_counts$Subject == j & OTU_counts$Time == i-1,names(hyp_vals)] == 0)]
      
      not_zero <- c(sample(prev_zero[prev_zero %in% small_otus], size = round(non_zero_small_prob*length(prev_zero[prev_zero %in% small_otus]))),
                    sample(prev_zero[prev_zero %in% large_otus], size = round(non_zero_large_prob*length(prev_zero[prev_zero %in% large_otus]))))
      
      stay_zero <- prev_zero[!(prev_zero %in% not_zero)]
      
      #make sure all 'stay_zero' OTUs are zero (save current counts to redistribute)
      count_redist <- sum(tmp_vec[stay_zero])
      tmp_vec[stay_zero] <- 0
      
      #what OTUs are currently zero that do not belong to the 'not_zero' group
      stay_zero <- unique(c(stay_zero,names(hyp_vals)[tmp_vec[names(hyp_vals)]==0][!(names(hyp_vals)[tmp_vec[names(hyp_vals)]==0] %in% not_zero)]))
      
      #How many otus need to be added to this vector to reach the number of zeros?
      zeros_to_add <- n_zeros_cur - length(stay_zero)
      
      #Use sorted previous counts and make the 'zeros_to_add' smallest OTUs zero.
      
      #Take list of OTUs minus 'not_zero' otus and 'stay_zero' OTUS and select new OTUs to be zero:
      otus_consider <- names(hyp_vals)[!(names(hyp_vals) %in% c(stay_zero, not_zero))]
      if(length(otus_consider) >= zeros_to_add  & zeros_to_add > 0){
        stay_zero <- c(stay_zero,names(sort(OTU_counts[OTU_counts$Subject == j & OTU_counts$Time == i-1, all_of(otus_consider)]))[1:zeros_to_add])
      }
      
      if(length(otus_consider) < zeros_to_add){
        stay_zero <- c(stay_zero,names(sort(OTU_counts[OTU_counts$Subject == j & OTU_counts$Time == i-1, all_of(otus_consider)]))[1:length(otus_consider)])
      }
      
      if(zeros_to_add < 0){
        #remove OTUs that had the most counts at previous time point
        zeros_to_add <- -1*zeros_to_add
        not_zero <- c(not_zero,names(sort(OTU_counts[OTU_counts$Subject == j & OTU_counts$Time == i-1, all_of(stay_zero)],decreasing = T)[1:zeros_to_add]))
      }
      
      count_redist <- count_redist + sum(tmp_vec[stay_zero])
      tmp_vec[stay_zero] <- 0
      
      #Now redistribute the counts so that any OTU that is currently zero and should NOT be zero has at
      # least 1 count:
      no_longer_zero <- not_zero[which(tmp_vec[not_zero] == 0)]
      tmp_vec[no_longer_zero] <- 1
      
      count_redist <- count_redist - length(no_longer_zero)
      
      #Now redistribute the counts to all of the non-zero OTUs based on eta_props vector:
      if(count_redist <= 0){
        count_redist <- 40
      }
      
      tmp_additions <- (rmultinom(1, size = count_redist, prob = eta_props_new[names(hyp_vals)[!(names(hyp_vals) %in% stay_zero)]]))
      
      for(k in rownames(tmp_additions)){
        tmp_vec[k] <- tmp_vec[k] + tmp_additions[k,1]
      }
      
      OTU_counts <- rbind(OTU_counts, tmp_vec)
      
      within_bd_tracker <- c(within_bd_tracker, mean(vegdist(OTU_counts[OTU_counts$Subject == j,names(hyp_vals)],method = bd_metric)))
      individual_bd_matrix <- rbind(individual_bd_matrix, cbind(Subject = j, Time_Comp = paste(i-1," to ", i, sep = ""),
                                                                Beta_Diversity = vegdist(OTU_counts[OTU_counts$Subject == j & (OTU_counts$Time == i | OTU_counts$Time == i-1),names(hyp_vals)],method = bd_metric)))
      
    }
    
    tmp_OTU_mat <- OTU_counts[OTU_counts$Time == i,names(hyp_vals)]
    
    avg_beta_div <- rbind(avg_beta_div, cbind(Time = i, Avg_between_BD= mean(vegdist(tmp_OTU_mat, method = bd_metric)),
                                              Avg_within_BD = mean(within_bd_tracker)))
    
  } 
  
  out_return <- list(Data = OTU_counts, Perturbed_OTUs = mean_otus, Large_OTUs = large_otus, 
                     Small_OTUS = small_otus, Beta_Div = avg_beta_div, Individual_BD = individual_bd_matrix)
  
  return(out_return)	    
} 

#df needs to have rows as samples, and columns for otus, column for time called "Time", and subject ID
stacked_bar <- function(df, otus, horizontal = "Time", vertical = "Relative Abundance", col_vals = NULL){
  sim_1_props <- df
  
  sim_1_props[,otus] <- t(apply(df[,all_of(otus)], 1, function(x) x/sum(x)))
  sim_1_props_avg <- c()
  
  for(i in 1:length(unique(df[,horizontal]))){
    tmp_df <- sim_1_props[sim_1_props[,horizontal] == unique(df[,horizontal])[i],]
    sim_1_tmp <- cbind(X = i, t(apply(tmp_df[,otus], 2, mean)))
    sim_1_props_avg <- rbind(sim_1_props_avg,
                             sim_1_tmp)
  }
  
  sim_1_props_avg <- as.data.frame(sim_1_props_avg)
  colnames(sim_1_props_avg)[1] <- horizontal
  
  sim_1_long <- sim_1_props_avg %>% pivot_longer(cols = all_of(otus), names_to = "OTU", values_to = "Count")
  
  sim_1_long <- as.data.frame(sim_1_long)
  
  sim_1_long_last <- sim_1_long[sim_1_long[,horizontal] == unique(df[,horizontal])[i],]
  otus <- sim_1_long_last[order(sim_1_long_last[,"Count"]), "OTU"]
  
  sim_1_long$OTU <- factor(sim_1_long$OTU, levels = otus)
    
  g <- ggplot(sim_1_long , aes(x = sim_1_long[,horizontal], y = Count, fill = OTU)) + geom_bar(position = "fill", stat = "identity") +
    theme_bw() + theme(legend.position = "none") +
    labs(x = horizontal,
         y = vertical,
         title = paste("Average Microbiome composition at each ", horizontal,"\nby OTU",sep = ""))
  
  if(is.null(col_vals)){
    ct <- 12
    ot <- length(otus)
    col_list <- c(rep(brewer.pal(n = 12, name = "Paired"), ot/ct), brewer.pal(n = ot%%ct, name = "Paired"))
    g <- g + scale_fill_manual(values = col_list)
  }
  if(!is.null(col_vals)){
    ct <- length(col_vals)
    ot <- length(otus)
    col_list <- c(rep(col_vals, ot/ct), col_vals[1:(ot%%ct)])
    names(col_list) <- otus
    g <- g + scale_fill_manual(values = col_list)
  }
  
  return(g)
}