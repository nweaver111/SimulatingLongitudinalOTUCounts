################################################################################
##          Data Generation (Adapted for Server from Data Generation file)    ##
##                                                                            ##
##                Created: Nick Weaver, April 2021                            ##
##                                                                            ##
##                Update: July 3, 2021 to make sure everything is correct     ##
##                        -removed the seed selection line because unused     ##
##                        -now have 2 items in output, adjusted for this!     ##
##                        -now does all scenarios (including time trends!)    ##
##                                                                            ##
################################################################################

######## Prep Work ########

#Load WF 12 week OTU counts:
OTU_dat_12 <- read.csv("C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\WF_data.csv",
                       header = T, row.names = 1)

View(OTU_dat_12)

#Load dirichlet values:
vals.dat <- readRDS(file = "C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\hyp_vals.rds")

#Manipulate the dirmult information to prepare for simulation:
gamma <- vals.dat$pi
theta <- vals.dat$theta

hyp_vals <- gamma*((1-theta)/theta)

avg_change <- readRDS(file = "C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\avg_change.rds")

#Load the simulation code:
source("C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\Large Simulation Code.R")

########## OTU Simulation #############
##
## Goal:-Get simulated OTU counts needed to run all simulation frameworks
##      -Get data for all effect estimate combinations and each framework
##

###Identify Simulation Scenario to use
# 'scenario' from R code being run

#Possible scenario options:
# "CSKAT_Main"
# "CSKAT_Confounding"
# "CSKAT_30"
# "CSKAT_500"
# "CSKAT_Subject"
# "Nick_1_positive"
# "Nick_1_negative"
# "Nick_1_parabolic"
# "Nick_2"
# "Nick_3"
# "Nick_4"

### Simulate 100 samples at 7 time points for 210 OTUs using the correct scenario
seeeed <- seed_offset
set.seed(seeeed)

run <- as.numeric(rep_start):(as.numeric(rep_start) + as.numeric(replicates) - 1)
tracker <- 1

for(ru in run){
  
  lop <- ifelse(tracker <= 100, "1_100",
                ifelse(tracker <= 200, "101_200",
                       ifelse(tracker <= 300, "201_300",
                              ifelse(tracker <= 400, "301_400",
                                     ifelse(tracker <= 500, "401_500",
                                            ifelse(tracker <= 600, "501_600",
                                                   ifelse(tracker <= 700, "601_700",
                                                          ifelse(tracker <= 800, "701_800",
                                                                 ifelse(tracker <= 900, "801_900", "901_1000")
                                                          )
                                                   )
                                            )
                                     )
                              )
                       )
                )
  )
  
  tracker <- tracker + 1

#Create Method results object:
res <- data.frame(Run = NA, Method = NA, Statistic = NA, P_Value = NA)
scenario <- "CSKAT_Main"
a_set <- "Perturbed"

for(ai in 1:10){
  start_time <- Sys.time()
  sim_dat_info <- c()
  while(length(sim_dat_info)==0){
    tryCatch({
      sim_dat_info <- SimOTUCounts(Real_Data = OTU_dat_12, hyp_vals = hyp_vals, sim_scen="low", time_length = 7, nsub = 100,
                                   lambda = .01, #cutoff for small/large OTU designation
                                   eta = 0.5, #proportion of large OTUs that will be perturbed
                                   prop_cover = c(0.5, 0.75, 1, 0.75, 0.5, 0)) #controls change at each time point
    },error = function(e){
      
    })
  }
 
  end_time <- Sys.time()

  sim_dat <- sim_1 <- sim_dat_info$Data
  sim_1$zeros <- apply(sim_1[,3:197], 1, function(x) sum(x == 0))
  
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  ggplot(sim_dat_info$Beta_Div[-1,], aes(x = Time, y = Avg_within_BD)) + theme_bw() + 
    geom_point(aes(size=2)) + geom_point(aes(y = Avg_between_BD,color = "red", size = 2)) + 
    geom_line() + geom_line(aes(y = Avg_between_BD)) + 
    labs(y = "Average Beta Diversity", color = "Between Subject Beta Diversity")
  
  ggplot(sim_dat_info$Individual_BD[-1,], aes(x = Time_Comp, y = as.numeric(Beta_Diversity), group = Subject)) + 
    theme_bw() + geom_point(color = "gray") + geom_line(color = "gray") +
    labs(x = "Time Gap",
         y = "Beta Diversity",
         title = "Subject trends for Beta Diversity")
  
  ggplot(sim_dat_info$Beta_Div[-1,], aes(x = Time, y = as.numeric(Avg_within_BD))) + theme_bw() + 
    geom_point(aes(size=2)) + geom_point(aes(y = Avg_between_BD,color = "red", size = 2)) + 
    geom_line() + geom_line(aes(y = Avg_between_BD)) + 
    labs(y = "Average Beta Diversity", color = "Between Subject Beta Diversity")

  ggplot(sim_1, aes(x = as.factor(Time), y = zeros)) + theme_bw() +
    geom_boxplot() +
    labs(x = "Timepoint",
         y = "Number of OTUs\nwith no Counts",
         title = "Distribution of OTUs that are simulated\nto be zero at each Timepoint")
  ggsave(paste("C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\large_OTU_plots\\error\\Zeros_by_Timepoint_large_Perturbed_simulation_",ai,".png",sep = ""),
         width = 7, heigh = 7)
  
  otus <- names(sim_1)[3:197]
  
  
  stacked_bar(df = sim_1, otus = otus, col_vals = c("blue", "red", "green", "yellow",
                                                    "orange", "gray", "black",
                                                    "pink", "brown", "purple"))
  
  stacked_bar(df = sim_1, otus = otus)
  
  
  
  
  sim_1_props <- sim_1
  
  sim_1_props[,otus] <- t(apply(sim_1[,all_of(otus)], 1, function(x) x/sum(x)))
  
  sim_1_props_avg_1 <- sim_1_props %>% filter(Time == 1) 
  sim_1_props_avg_2 <- sim_1_props %>% filter(Time == 2)
  sim_1_props_avg_3 <- sim_1_props %>% filter(Time == 3)
  sim_1_props_avg_4 <- sim_1_props %>% filter(Time == 4)
  sim_1_props_avg_5 <- sim_1_props %>% filter(Time == 5)
  sim_1_props_avg_6 <- sim_1_props %>% filter(Time == 6)
  sim_1_props_avg_7 <- sim_1_props %>% filter(Time == 7)
  
  sim_1_props_avg<- cbind(Time = 1, t(apply(sim_1_props_avg_1[,otus], 2, mean)))
  sim_1_props_avg<- rbind(sim_1_props_avg,
                          cbind(Time = 2, t(apply(sim_1_props_avg_2[,otus], 2, mean))),
                          cbind(Time = 3, t(apply(sim_1_props_avg_3[,otus], 2, mean))),
                          cbind(Time = 4, t(apply(sim_1_props_avg_4[,otus], 2, mean))),
                          cbind(Time = 5, t(apply(sim_1_props_avg_5[,otus], 2, mean))),
                          cbind(Time = 6, t(apply(sim_1_props_avg_6[,otus], 2, mean))),
                          cbind(Time = 7, t(apply(sim_1_props_avg_7[,otus], 2, mean))))
  
  sim_1_props_avg <- as.data.frame(sim_1_props_avg)
  
  sim_1_long <- sim_1_props_avg %>% pivot_longer(cols = all_of(otus), names_to = "OTU", values_to = "Count")
  
  sim_1_long$Type <- "Large"
  sim_1_long[sim_1_long$OTU %in% sim_dat_info$Perturbed_OTUs, "Type"] <- "Perturbed"
  sim_1_long[sim_1_long$OTU %in% sim_dat_info$Small_OTUS, "Type"] <- "Small"
  
  ggplot(sim_1_long , aes(x = Time, y = Count, fill = OTU)) + geom_bar(position = "fill", stat = "identity") +
    theme_bw() + theme(legend.position = "none") +
    labs(x = "Timepoint",
         y = "Relative abundance",
         title = "Average Microbiome composition at each Timepoint\nby OTU")
  ggsave(paste("C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\large_OTU_plots\\error\\avg_composition_large_perturbed_by_OTU_sim_",ai,".png",sep=""),
         width = 7, heigh = 7)
  
  ggplot(sim_1_long , aes(x = Time, y = Count, fill = Type)) + geom_bar(position = "fill", stat = "identity") +
    theme_bw() +
    labs(x = "Timepoint",
         y = "Relative Abundance",
         title = "Average Microbiome composition at each Timepoint\nby OTU grouping")
  ggsave(paste("C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\large_OTU_plots\\error\\avg_composition_large_perturbed_by_grouping_sim_",ai,".png",sep=""),
         width = 7, heigh = 7)
  
  sim_1_long_perturbed <- sim_1_long %>% filter(OTU %in% c(sim_dat_info$Large_OTUs[1],all_of(sim_dat_info$Perturbed_OTUs)))
  
  ggplot(sim_1_long_perturbed , aes(x = Time, y = Count, fill = OTU)) + geom_bar(position = "fill", stat = "identity") +
    theme_bw() + theme(legend.position = "none") +
    labs(x = "Timepoint",
         y = "Relative Abundance",
         title = "Average abundance of perturbed OTUs at each Timepoint\nby OTU\n(Including 1 large OTU for reference)")
  ggsave(paste("C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\large_OTU_plots\\error\\avg_composition_perturbedOTUsOnly_large_perturbed_by_OTU_sim_",ai,".png",sep=""),
         width = 7, heigh = 7)
  
  ggplot(sim_1_long_perturbed , aes(x = Time, y = Count, fill = Type)) + geom_bar(position = "fill", stat = "identity") +
    theme_bw() +
    labs(x = "Timepoint",
         y = "Relative Abundance",
         title = "Average Abundance of perturbed OTUs at each Timepoint\nby OTU grouping\n(Including 1 large OTU for reference)")
  ggsave(paste("C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\large_OTU_plots\\error\\avg_composition_perturbedOTUsOnly_large_perturbed_by_grouping_sim_",ai,".png",sep=""),
         width = 7, heigh = 7)
  
  sim_1_long_large <- sim_1_long %>% filter(OTU %in% all_of(sim_dat_info$Large_OTUs))
  
  ggplot(sim_1_long_large , aes(x = Time, y = Count, fill = OTU)) + geom_bar(position = "fill", stat = "identity") +
    theme_bw() + theme(legend.position = "none") +
    labs(x = "Timepoint",
         y = "Relative Abundance",
         title = "Average Abundance of large OTUs at each Timepoint\nby OTU")
  ggsave(paste("C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\large_OTU_plots\\error\\avg_composition_largeOTUsOnly_large_perturbed_by_OTU_sim_",ai,".png",sep=""),
         width = 7, heigh = 7)
  
  
  sim_1_long_zero <- sim_1_long %>% filter(OTU %in% c(sim_dat_info$Perturbed_OTUs[1],all_of(sim_dat_info$Small_OTUS)))
  
  ggplot(sim_1_long_zero , aes(x = Time, y = Count, fill = OTU)) + geom_bar(position = "fill", stat = "identity") +
    theme_bw() + theme(legend.position = "none") +
    labs(x = "Timepoint",
         y = "Relative Abundance",
         title = "Average Abundance of small OTUs at each Timepoint\nby OTU grouping\n(Including 1 perturbed OTU for reference)")
  ggsave(paste("C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\large_OTU_plots\\error\\avg_composition_zeroOTUsOnly_large_perturbed_by_OTU_sim_",ai,".png",sep=""),
         width = 7, heigh = 7)
  
  ggplot(sim_1_long_zero , aes(x = Time, y = Count, fill = Type)) + geom_bar(position = "fill", stat = "identity") +
    theme_bw() +
    labs(x = "Timepoint",
         y = "Relative Abundance",
         title = "Average Abundance of zero OTUs at each Timepoint\nby OTU grouping\n(Including 1 perturbed OTU for reference)")
  ggsave(paste("C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\large_OTU_plots\\error\\avg_composition_zeroOTUsOnly_large_perturbed_by_grouping_sim_",ai,".png",sep=""),
         width = 7, heigh = 7)
  
  


  perturbed_otus <- sim_dat_info$Perturbed_OTUs
  sim_dat <- sim_dat_info$Data

  otus <- colnames(sim_dat)[3:length(sim_dat[1,])] #OTU names.


  ####Simulate the remaining covariates and/or noise parameters:

  #Simulate the 2 covariates:

  #covar 1 (subject dependent)
  for(l in 1:length(unique(sim_dat$Subject))){
    sim_dat[sim_dat$Subject == l, "Covar1"] <- rbinom(1, 1, prob = 0.5) #Same as CSKAT
  }

  #covar 2 (different for all observations)
  sim_dat[,"Covar2"] <- rnorm(length(sim_dat[,"Subject"])) #Same as CSKAT



  #Find active set by abundance:
  
  #Subset to 20% of most prevalent OTUs then randomly select half of them
  # - This will make the active set consist of 10% of the OTUs to match the random set.
  
  #Find relative abundance of each OTU in each site:
  props <- t(apply(sim_dat[,otus], 1, function(x) x/sum(x)))
  rownames(props) <- rownames(sim_dat)
  
  #Mean abundance of each OTU across all samples:
  mean_prop <- apply(props,2,mean)
  
  #Order by mean relative abundance:
  ord <- (order(mean_prop,decreasing = T))
  
  #subset to top 20%:
  nm <- .2*length(mean_prop)
  
  sum_OTU <- names(mean_prop[ord[1:nm]])
  
  #Now randomly select half of these OTUs for active cluster:
  act_OTU_clus <- sample(sum_OTU, size = round(.5*nm), replace = F) #List of OTU names in Active Set


  #Randomly select 10%
  act_OTU_rand <- sample(otus, size = .1*length(otus))


  #Randomly select same number of OTUs from perturbed list:
  if(length(perturbed_otus) >= length(act_OTU_clus)){
    act_OTU_pert <- sample(perturbed_otus, size = length(act_OTU_clus))
  }
  if(length(perturbed_otus) < length(act_OTU_clus)){
    act_OTU_pert <- perturbed_otus
  }

  #Store the vector for all three set types:
  for(i in 1:length(sim_dat$Subject)){
    sim_dat[i,"ClusteredEffect"] <- sum(sim_dat[i,act_OTU_clus])
    sim_dat[i,"RandMicrobiomeEffect"] <- sum(sim_dat[i,act_OTU_rand])
    sim_dat[i,"PertMicrobiomeEffect"] <- sum(sim_dat[i,act_OTU_pert])
  }


  #Scale the microbiome effect as described in CSKAT:
  sim_dat[,"ScaledClusteredEffect"] <- sapply(sim_dat$ClusteredEffect, function(x) (x-mean(sim_dat$ClusteredEffect))/sd(sim_dat$ClusteredEffect))
  sim_dat[,"ScaledRandMicrobiomeEffect"] <- sapply(sim_dat$RandMicrobiomeEffect, function(x) (x-mean(sim_dat$RandMicrobiomeEffect))/sd(sim_dat$RandMicrobiomeEffect))
  sim_dat[,"ScaledPertMicrobiomeEffect"] <- sapply(sim_dat$PertMicrobiomeEffect, function(x) (x-mean(sim_dat$PertMicrobiomeEffect))/sd(sim_dat$PertMicrobiomeEffect))


  #Now all X matrix data is ready to go, combine relevant information together:
  sim_dat_clustered <- sim_dat[,c("Subject", "Time", "Covar1", "Covar2", "ScaledClusteredEffect")]
  sim_dat_random <- sim_dat[,c("Subject", "Time", "Covar1", "Covar2", "ScaledRandMicrobiomeEffect")]
  sim_dat_perturbed <- sim_dat[,c("Subject", "Time", "Covar1", "Covar2", "ScaledPertMicrobiomeEffect")]

  ####Additional confounders ####

  #Simulate covar 2 so that it is a N(0,1) + the microbiome effect (scale(sum(activeset values)))
  sim_dat_clustered[,"Covar2_confounding"] <- sim_dat$ScaledClusteredEffect + rnorm(length(sim_dat[,"Subject"])) 
  sim_dat_random[,"Covar2_confounding"] <- sim_dat$ScaledRandMicrobiomeEffect + rnorm(length(sim_dat[,"Subject"]))
  sim_dat_perturbed[,"Covar2_confounding"] <- sim_dat$ScaledPertMicrobiomeEffect + rnorm(length(sim_dat[,"Subject"]))


  #CSKAT_Main scenario:
  if(scenario == "CSKAT_Main"){
    #Identify the active set:
    if(a_set == "Random"){
      tmp_dat <- sim_dat_random
      tmp_val <- "ScaledRandMicrobiomeEffect"
    }
  
    if(a_set == "Cluster"){
      tmp_dat <- sim_dat_clustered
      tmp_val <- "ScaledClusteredEffect"
    }
  
    if(a_set == "Perturbed"){
      tmp_dat <- sim_dat_perturbed
      tmp_val <- "ScaledPertMicrobiomeEffect"
    }
  
    j <- 0.5
    k <- 4
    i <- 1
  
    #tmp_nam <- paste("/home/math/weavenic/Dissertation/Simulations/Aim1/2021_08_10/Output/",scenario,"/",scenario,"_Simulations_RepStart_",rep_start,"_Beta_",k,"_randomSubject_",i,"_noise_",j,"_active_",a_set,"_seedOffset_",seed_offset,"_pertType_",pert_type,"/",lop,"/Run_",ru,sep="")
    #tmp_nam_act <- paste("/home/math/weavenic/Dissertation/Simulations/Aim1/2021_08_10/Active_sets/",scenario,"/",scenario,"_Simulations_RepStart_",rep_start,"_Beta_",k,"_randomSubject_",i,"_noise_",j,"_active_",a_set,"_seedOffset_",seed_offset,"_pertType_",pert_type,"/",lop,"/Run_",ru,sep="")
   
    #Noise:
    tmp_dat[,"Noise"] <- rnorm(length(tmp_dat[,"Subject"]), 0, j)
        
    #Subject Effect:
    for(l in 1:length(unique(tmp_dat$Subject))){
      tmp_dat[tmp_dat$Subject == l,"RandomEffect"] <- rnorm(1, 0, i)
    }
        
    #Microbiome Effect:
    tmp_dat[,"MicrobiomeEffect"] <- k*tmp_dat[,tmp_val]
        
    ###### Outcome ######
    tmp_dat[,"Outcome"] <- 0.5*(1 + tmp_dat[,"Covar1"]+tmp_dat[,"Covar2"]) +
    tmp_dat[,"MicrobiomeEffect"] + tmp_dat[,"RandomEffect"] + tmp_dat[,"Noise"]
        
    #Add OTU counts information to output data file.
    tmp_dat[,otus] <- sim_dat[,otus]
        
    #write.csv(tmp_dat, file = gzfile(paste(tmp_nam,".csv.gz",sep="")))

    pert_lists <- list(Random_Group = act_OTU_rand, Clustered_Group = act_OTU_clus)
    #saveRDS(pert_lists, file = paste(tmp_nam_act,"_ActiveSets.rds", sep=""))
  
  
    ###Testing###
  
    
    #Turn everything into proportions
    sum.sim <- apply(tmp_dat[,otus], 1, sum) #length should equal number of samples simulated
    
    
    props.sim <- tmp_dat[,otus]/sum.sim
    props.sim <- cbind(tmp_dat[,1:2], props.sim)
    
    
    ##Make some plots:
    
    #Try boxplot of zero counts:
    zeros.df.1 <- data.frame(Simulation = 1, Bin = "Zeros",Counts = apply(props.sim[props.sim$Time == 1,], 1, function(x) sum(x==0))/length(otus))
    zeros.df.2 <- data.frame(Simulation = 2, Bin = "Zeros",Counts = apply(props.sim[props.sim$Time == 2,], 1, function(x) sum(x==0))/length(otus))
    zeros.df.3 <- data.frame(Simulation = 3, Bin = "Zeros",Counts = apply(props.sim[props.sim$Time == 3,], 1, function(x) sum(x==0))/length(otus))
    zeros.df.4 <- data.frame(Simulation = 4, Bin = "Zeros",Counts = apply(props.sim[props.sim$Time == 4,], 1, function(x) sum(x==0))/length(otus))
    zeros.df.5 <- data.frame(Simulation = 5, Bin = "Zeros",Counts = apply(props.sim[props.sim$Time == 5,], 1, function(x) sum(x==0))/length(otus))
    zeros.df.6 <- data.frame(Simulation = 6, Bin = "Zeros",Counts = apply(props.sim[props.sim$Time == 6,], 1, function(x) sum(x==0))/length(otus))
    zeros.df.7 <- data.frame(Simulation = 7, Bin = "Zeros",Counts = apply(props.sim[props.sim$Time == 7,], 1, function(x) sum(x==0))/length(otus))
    
    
    zeros.df.1_2 <- data.frame(Simulation = 1, Bin = "<.001",Counts = apply(props.sim[props.sim$Time == 1,], 1, function(x) sum(0 < x & x < 0.001))/length(otus))
    zeros.df.2_2 <- data.frame(Simulation = 2, Bin = "<.001",Counts = apply(props.sim[props.sim$Time == 2,], 1, function(x) sum(0 < x & x < 0.001))/length(otus))
    zeros.df.3_2 <- data.frame(Simulation = 3, Bin = "<.001",Counts = apply(props.sim[props.sim$Time == 3,], 1, function(x) sum(0 < x & x < 0.001))/length(otus))
    zeros.df.4_2 <- data.frame(Simulation = 4, Bin = "<.001",Counts = apply(props.sim[props.sim$Time == 4,], 1, function(x) sum(0 < x & x < 0.001))/length(otus))
    zeros.df.5_2 <- data.frame(Simulation = 5, Bin = "<.001",Counts = apply(props.sim[props.sim$Time == 5,], 1, function(x) sum(0 < x & x < 0.001))/length(otus))
    zeros.df.6_2 <- data.frame(Simulation = 6, Bin = "<.001",Counts = apply(props.sim[props.sim$Time == 6,], 1, function(x) sum(0 < x & x < 0.001))/length(otus))
    zeros.df.7_2 <- data.frame(Simulation = 7, Bin = "<.001",Counts = apply(props.sim[props.sim$Time == 7,], 1, function(x) sum(0 < x & x < 0.001))/length(otus))
    
    
    zeros.df.1_3 <- data.frame(Simulation = 1, Bin = "<.005",Counts = apply(props.sim[props.sim$Time == 1,], 1, function(x) sum(0.001 <= x & x <0.005))/length(otus))
    zeros.df.2_3 <- data.frame(Simulation = 2, Bin = "<.005",Counts = apply(props.sim[props.sim$Time == 2,], 1, function(x) sum(0.001 <= x & x <0.005))/length(otus))
    zeros.df.3_3 <- data.frame(Simulation = 3, Bin = "<.005",Counts = apply(props.sim[props.sim$Time == 3,], 1, function(x) sum(0.001 <= x & x <0.005))/length(otus))
    zeros.df.4_3 <- data.frame(Simulation = 4, Bin = "<.005",Counts = apply(props.sim[props.sim$Time == 4,], 1, function(x) sum(0.001 <= x & x <0.005))/length(otus))
    zeros.df.5_3 <- data.frame(Simulation = 5, Bin = "<.005",Counts = apply(props.sim[props.sim$Time == 5,], 1, function(x) sum(0.001 <= x & x <0.005))/length(otus))
    zeros.df.6_3 <- data.frame(Simulation = 6, Bin = "<.005",Counts = apply(props.sim[props.sim$Time == 6,], 1, function(x) sum(0.001 <= x & x <0.005))/length(otus))
    zeros.df.7_3 <- data.frame(Simulation = 7, Bin = "<.005",Counts = apply(props.sim[props.sim$Time == 7,], 1, function(x) sum(0.001 <= x & x <0.005))/length(otus))
    
    
    
    zeros.df.1_4 <- data.frame(Simulation = 1, Bin = "<.01",Counts = apply(props.sim[props.sim$Time == 1,], 1, function(x) sum(0.005 <= x & x < 0.01))/length(otus))
    zeros.df.2_4 <- data.frame(Simulation = 2, Bin = "<.01",Counts = apply(props.sim[props.sim$Time == 2,], 1, function(x) sum(0.005 <= x & x < 0.01))/length(otus))
    zeros.df.3_4 <- data.frame(Simulation = 3, Bin = "<.01",Counts = apply(props.sim[props.sim$Time == 3,], 1, function(x) sum(0.005 <= x & x < 0.01))/length(otus))
    zeros.df.4_4 <- data.frame(Simulation = 4, Bin = "<.01",Counts = apply(props.sim[props.sim$Time == 4,], 1, function(x) sum(0.005 <= x & x < 0.01))/length(otus))
    zeros.df.5_4 <- data.frame(Simulation = 5, Bin = "<.01",Counts = apply(props.sim[props.sim$Time == 5,], 1, function(x) sum(0.005 <= x & x < 0.01))/length(otus))
    zeros.df.6_4 <- data.frame(Simulation = 6, Bin = "<.01",Counts = apply(props.sim[props.sim$Time == 6,], 1, function(x) sum(0.005 <= x & x < 0.01))/length(otus))
    zeros.df.7_4 <- data.frame(Simulation = 7, Bin = "<.01",Counts = apply(props.sim[props.sim$Time == 7,], 1, function(x) sum(0.005 <= x & x < 0.01))/length(otus))
    
    
    
    zeros.df.1_5 <- data.frame(Simulation = 1, Bin = "<.025",Counts = apply(props.sim[props.sim$Time == 1,], 1, function(x) sum(0.01 <= x & x <0.025))/length(otus))
    zeros.df.2_5 <- data.frame(Simulation = 2, Bin = "<.025",Counts = apply(props.sim[props.sim$Time == 2,], 1, function(x) sum(0.01 <= x & x <0.025))/length(otus))
    zeros.df.3_5 <- data.frame(Simulation = 3, Bin = "<.025",Counts = apply(props.sim[props.sim$Time == 3,], 1, function(x) sum(0.01 <= x & x <0.025))/length(otus))
    zeros.df.4_5 <- data.frame(Simulation = 4, Bin = "<.025",Counts = apply(props.sim[props.sim$Time == 4,], 1, function(x) sum(0.01 <= x & x <0.025))/length(otus))
    zeros.df.5_5 <- data.frame(Simulation = 5, Bin = "<.025",Counts = apply(props.sim[props.sim$Time == 5,], 1, function(x) sum(0.01 <= x & x <0.025))/length(otus))
    zeros.df.6_5 <- data.frame(Simulation = 6, Bin = "<.025",Counts = apply(props.sim[props.sim$Time == 6,], 1, function(x) sum(0.01 <= x & x <0.025))/length(otus))
    zeros.df.7_5 <- data.frame(Simulation = 7, Bin = "<.025",Counts = apply(props.sim[props.sim$Time == 7,], 1, function(x) sum(0.01 <= x & x <0.025))/length(otus))
    
    
    zeros.df.1_6 <- data.frame(Simulation = 1, Bin = "<.05",Counts = apply(props.sim[props.sim$Time == 1,], 1, function(x) sum(0.025 <= x & x <0.05))/length(otus))
    zeros.df.2_6 <- data.frame(Simulation = 2, Bin = "<.05",Counts = apply(props.sim[props.sim$Time == 2,], 1, function(x) sum(0.025 <= x & x <0.05))/length(otus))
    zeros.df.3_6 <- data.frame(Simulation = 3, Bin = "<.05",Counts = apply(props.sim[props.sim$Time == 3,], 1, function(x) sum(0.025 <= x & x <0.05))/length(otus))
    zeros.df.4_6 <- data.frame(Simulation = 4, Bin = "<.05",Counts = apply(props.sim[props.sim$Time == 4,], 1, function(x) sum(0.025 <= x & x <0.05))/length(otus))
    zeros.df.5_6 <- data.frame(Simulation = 5, Bin = "<.05",Counts = apply(props.sim[props.sim$Time == 5,], 1, function(x) sum(0.025 <= x & x <0.05))/length(otus))
    zeros.df.6_6 <- data.frame(Simulation = 6, Bin = "<.05",Counts = apply(props.sim[props.sim$Time == 6,], 1, function(x) sum(0.025 <= x & x <0.05))/length(otus))
    zeros.df.7_6 <- data.frame(Simulation = 7, Bin = "<.05",Counts = apply(props.sim[props.sim$Time == 7,], 1, function(x) sum(0.025 <= x & x <0.05))/length(otus))
    
    
    
    zeros.df.1_7 <- data.frame(Simulation = 1, Bin = "<.10",Counts = apply(props.sim[props.sim$Time == 1,], 1, function(x) sum(0.05 <= x & x <0.1))/length(otus))
    zeros.df.2_7 <- data.frame(Simulation = 2, Bin = "<.10",Counts = apply(props.sim[props.sim$Time == 2,], 1, function(x) sum(0.05 <= x & x <0.1))/length(otus))
    zeros.df.3_7 <- data.frame(Simulation = 3, Bin = "<.10",Counts = apply(props.sim[props.sim$Time == 3,], 1, function(x) sum(0.05 <= x & x <0.1))/length(otus))
    zeros.df.4_7 <- data.frame(Simulation = 4, Bin = "<.10",Counts = apply(props.sim[props.sim$Time == 4,], 1, function(x) sum(0.05 <= x & x <0.1))/length(otus))
    zeros.df.5_7 <- data.frame(Simulation = 5, Bin = "<.10",Counts = apply(props.sim[props.sim$Time == 5,], 1, function(x) sum(0.05 <= x & x <0.1))/length(otus))
    zeros.df.6_7 <- data.frame(Simulation = 6, Bin = "<.10",Counts = apply(props.sim[props.sim$Time == 6,], 1, function(x) sum(0.05 <= x & x <0.1))/length(otus))
    zeros.df.7_7 <- data.frame(Simulation = 7, Bin = "<.10",Counts = apply(props.sim[props.sim$Time == 7,], 1, function(x) sum(0.05 <= x & x <0.1))/length(otus))
    
    
    
    zeros.df.1_8 <- data.frame(Simulation = 1, Bin = "<1",Counts = apply(props.sim[props.sim$Time == 1,], 1, function(x) sum(0.1 <= x))/length(otus))
    zeros.df.2_8 <- data.frame(Simulation = 2, Bin = "<1",Counts = apply(props.sim[props.sim$Time == 2,], 1, function(x) sum(0.1 <= x))/length(otus))
    zeros.df.3_8 <- data.frame(Simulation = 3, Bin = "<1",Counts = apply(props.sim[props.sim$Time == 3,], 1, function(x) sum(0.1 <= x))/length(otus))
    zeros.df.4_8 <- data.frame(Simulation = 4, Bin = "<1",Counts = apply(props.sim[props.sim$Time == 4,], 1, function(x) sum(0.1 <= x))/length(otus))
    zeros.df.5_8 <- data.frame(Simulation = 5, Bin = "<1",Counts = apply(props.sim[props.sim$Time == 5,], 1, function(x) sum(0.1 <= x))/length(otus))
    zeros.df.6_8 <- data.frame(Simulation = 6, Bin = "<1",Counts = apply(props.sim[props.sim$Time == 6,], 1, function(x) sum(0.1 <= x))/length(otus))
    zeros.df.7_8 <- data.frame(Simulation = 7, Bin = "<1",Counts = apply(props.sim[props.sim$Time == 7,], 1, function(x) sum(0.1 <= x))/length(otus))
    
    
    zeros.df <- rbind(zeros.df.1, zeros.df.2, zeros.df.3, zeros.df.4, zeros.df.5, zeros.df.6, zeros.df.7,
                      zeros.df.1_2, zeros.df.2_2, zeros.df.3_2, zeros.df.4_2, zeros.df.5_2, zeros.df.6_2, zeros.df.7_2,
                      zeros.df.1_3, zeros.df.2_3, zeros.df.3_3, zeros.df.4_3, zeros.df.5_3, zeros.df.6_3, zeros.df.7_3,
                      zeros.df.1_4, zeros.df.2_4, zeros.df.3_4, zeros.df.4_4, zeros.df.5_4, zeros.df.6_4, zeros.df.7_4,
                      zeros.df.1_5, zeros.df.2_5, zeros.df.3_5, zeros.df.4_5, zeros.df.5_5, zeros.df.6_5, zeros.df.7_5,
                      zeros.df.1_6, zeros.df.2_6, zeros.df.3_6, zeros.df.4_6, zeros.df.5_6, zeros.df.6_6, zeros.df.7_6,
                      zeros.df.1_7, zeros.df.2_7, zeros.df.3_7, zeros.df.4_7, zeros.df.5_7, zeros.df.6_7, zeros.df.7_7,
                      zeros.df.1_8, zeros.df.2_8, zeros.df.3_8, zeros.df.4_8, zeros.df.5_8, zeros.df.6_8, zeros.df.7_8)
    
    
    
    
    ggplot(data = zeros.df, aes(y=Counts,fill=as.factor(Simulation), x=Bin)) +
      geom_boxplot(alpha = 0.3, width = .5) + 
      theme_bw() + 
      labs(x = "Proportion of single OTU in Sample", y ="Proportion of all OTUs within a Subject", title = "Simulation Comparisons\nDistribution of OTU proportions by Subject",fill="Timepoint") +
      scale_x_discrete(limits = c("Zeros", "<.001", "<.005", "<.01",
                                  "<.025", "<.05", "<.10", "<1"),
                       labels = c("Zeros" = "Zeros",
                                  "<.001"= "0 < x <.001", 
                                  "<.005"=".001 < x < .005", 
                                  "<.01"=".005 < x < .01", 
                                  "<.025"=".01 < x < .025", 
                                  "<.05"=".025 < x < .05", 
                                  "<.10"=".05 < x < .10", 
                                  "<1"=".10 < x < 1")) +
      theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, vjust = 0.5))
    ggsave(paste("C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\large_OTU_plots\\error\\large_perturbed_binned_counts_sim_",ai,".png",sep=""),
           width = 7, heigh = 7)  
    
    #Beta diveristy of OTUs?
    library(vegan)
    main <- vegdist(tmp_dat[order(tmp_dat$Subject),otus], method = "bray", diagonal = T, upper = T)

    ##Same thing, but with a subset to zoom in on heatmap!
    subset_main <- tmp_dat[tmp_dat$Subject == 6 | tmp_dat$Subject == 2 | tmp_dat$Subject == 68,]
    
    main_sub <- vegdist(subset_main[order(subset_main$Subject),otus], method = "bray", diagonal = T, upper = T)
    
    main_sub <- as.data.frame(as.matrix(main_sub))
    
    x <- colnames(main_sub)
    y <- colnames(main_sub)
    data <- expand.grid(X=x, Y=y)
    data$z <- c(main_sub[,1], main_sub[,2], main_sub[,3], main_sub[,4], main_sub[,5], main_sub[,6], main_sub[,7],
                main_sub[,8], main_sub[,9], main_sub[,10], main_sub[,11], main_sub[,12], main_sub[,13],main_sub[,14],
                main_sub[,15], main_sub[,16], main_sub[,17], main_sub[,18], main_sub[,19], main_sub[,20],main_sub[,21])
    
    # Heatmap 
    ggplot(data, aes(X, Y, fill= z)) + 
      geom_tile() +
      labs(x = "Observation",
           y = "Observation",
           fill = "Bray-Curtis Measure",
           title = "Bray-Curtis Dissimilarity Heatmap") +
      scale_fill_gradient(low="white", high="darkblue") + theme_bw() +
      theme(legend.position = "bottom")
    ggsave(paste("C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\large_OTU_plots\\error\\beta_div_heatmap_sim_",ai,".png",sep=""),
           width = 7, heigh = 7)
    
    
    #Hope to see the diagonal be more similar than other areas.
    
    ggplot(tmp_dat, aes(y = Outcome, x = ScaledPertMicrobiomeEffect)) + theme_bw() +
      geom_point()
    
    
  
    ####Methods####
    source("C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\SSKAT.R")
    library(CompQuadForm)
    library(vegan)
    library(ggplot2)
    library(lmerTest)
    library(stringr)
    library(metafor)
  
    ### Load Data ###
    dat <- tmp_dat
  
    dat <- dat[order(dat$Subject),]
    otus <- colnames(dat)[11:length(dat[1,])]
    dat$Subject <- as.factor(dat$Subject)
    dat$Covar1 <- as.factor(dat$Covar1)
    
    ### Use Relative Counts ###
    dat[,otus] <- (t(apply(dat[,otus], 1, function(x) x/sum(x))))
    #Filter so only use OTU with enough seen:
    otus <- names(which(apply(dat[,otus], 2, function(x) sum(x>0)>0.05) & 
                          apply(dat[,otus], 2, function(x) sum(x>.001)>0)))
    
    ##### CSKAT #####
  
    #Create Beta Diversity matrix
    K <- 1-as.matrix(vegdist(dat[,otus], method = "bray", diagonal = T, upper = T)) #Create a similarity measure
    aa <- cSKAT(Outcome ~ Covar1 + Covar2 + (1|Subject), data = dat, K = K)
  
    res <- rbind(res, cbind(Run = ai, Method = "CSKAT", Statistic = aa$Q.adj, P_Value = aa$p.value))
    
    
    
    ##### Beta Diversity as outcome #####
    
    #consecutive time points within subject only#
    bet_div <- c() #store beta diversity values here
  
    for(j in unique(dat$Subject)){
      sub_dat <- dat[dat$Subject==j,]
      bc <- as.matrix(vegdist(sub_dat[, otus]))
      tmp <- data.frame(Subject = as.vector(rep(j,6)), TimeDiff = as.vector(1:6), BetaDiversity = (bc[row(bc) == (col(bc) + 1)]),
                        Covar1 = as.vector(rep(unique(sub_dat[,"Covar1"]),6)))
      out_dif <- c()
      cov2_dif <- c()
      micro_dif <- c()
      for(k in 1:6){
        out_dif <- rbind(out_dif, sub_dat[sub_dat$Time == k+1,"Outcome"] - sub_dat[sub_dat$Time==k,"Outcome"])
        cov2_dif <- rbind(cov2_dif, sub_dat[sub_dat$Time == k+1,"Covar2"] - sub_dat[sub_dat$Time==k,"Covar2"])
        micro_dif <- rbind(micro_dif, sub_dat[sub_dat$Time == k+1,"ScaledPertMicrobiomeEffect"] - sub_dat[sub_dat$Time == k,"ScaledPertMicrobiomeEffect"])
      }
    
      bet_div <- rbind(bet_div,cbind(tmp, out_dif, cov2_dif, micro_dif))
    }
  
    ggplot(bet_div, aes(x = micro_dif, y = out_dif)) + theme_bw() + geom_point()
    ggplot(bet_div, aes(x = TimeDiff, y = BetaDiversity, group = Subject)) +
      theme_bw() + geom_point() + geom_line()
    ggplot(bet_div, aes(x = TimeDiff, y = out_dif, group = Subject)) +
      theme_bw() + geom_point() + geom_line()
    ggplot(bet_div, aes(x = micro_dif, y = BetaDiversity)) +theme_bw() +geom_point()
    ggplot(bet_div, aes(x = BetaDiversity, y = out_dif)) + theme_bw() +geom_point()
    ##Try a LMM for this data!
    tryCatch({
      lmod <- lmer(BetaDiversity ~ TimeDiff*out_dif + Covar1 + cov2_dif + (1|Subject), data = bet_div)
      bb <- summary(lmod)
    
      res <- rbind(res, cbind(Run = ai, Method = "LMM", Statistic = bb$coefficients["out_dif", "Estimate"], P_Value = bb$coefficients["out_dif", "Pr(>|t|)"]))
    },
    error=function(e){
    })
    
    
    ##Try quadratic LMM for this data!
    tryCatch({
      lmod_squared <- lmer(BetaDiversity ~ out_dif + I(out_dif^2) + Covar1 + cov2_dif + (1|Subject), data = bet_div)
      bb <- summary(lmod_squared)
    
      res <- rbind(res, cbind(Run = ai, Method = "LMM_Squared", Statistic = bb$coefficients["I(out_dif^2)", "Estimate"], P_Value = bb$coefficients["I(out_dif^2)", "Pr(>|t|)"]))
    },
    error=function(e){
    })
    
    
    ##Try two stage approach:
    tryCatch({
      first_stage <- lmList(BetaDiversity ~ cov2_dif + out_dif | Subject, data = bet_div)
    
      #Second stage for proof of concept:
      b <- lapply(first_stage, coef)
      V <- lapply(first_stage, vcov)
    
      estm <- rep(c("intercept", "covar2", "outcome"), length(b))
      subj <- rep(names(b), each=3)
    
      
      b <- unlist(b)
      V <- bldiag(V)
    
      res2 <- rma.mv(b ~ estm - 1, V, random = ~ estm | subj, struct="DIAG")
      cc <- summary(res2)
      res <- rbind(res, cbind(Run = ai, Method = "TwoStage", Statistic = cc$beta["estmoutcome",1], P_Value = cc$pval[3]))
    },
    error = function(e){
    })
    
    ##Try Permanova:
    perm = how(nperm = 1000)
    setBlocks(perm) <- with(dat, Subject) #Use this to permute within subjectID!
    t <- adonis2(dat[,otus] ~ Covar1 + Covar2 + Outcome, data = dat, permutations = perm, method = "bray")
  
    
    res <- rbind(res, cbind(Run = ai, Method = "PermANOVA", Statistic = t["Outcome", "F"], P_Value = t["Outcome", "Pr(>F)"]))
  

    print(ai)
  }
}
res <- res[-1,]

mod_perf <- "Power"


res$P_Value <- as.numeric(res$P_Value)
res$Statistic <- as.numeric(res$Statistic)


#Power/Type I error:
res[,mod_perf] <- ifelse(res$P_Value <0.05, 1, 0)

#CSKAT:
cskat_val <- sum(res[res$Method=="CSKAT", mod_perf])/sum(res$Method == "CSKAT")

#LMM:
lmm_val <- sum(res[res$Method == "LMM", mod_perf])/sum(res$Method == "LMM")

#LMM_Squared:
squared_val <- sum(res[res$Method == "LMM_Squared", mod_perf])/sum(res$Method == "LMM_Squared")

#TwoStage:
two_val <- sum(res[res$Method == "TwoStage", mod_perf])/sum(res$Method == "TwoStage")

#Permanova:
perm_val <- sum(res[res$Method == "PermANOVA", mod_perf])/sum(res$Method == "PermANOVA")

## Plots ##

ggplot(res, aes(x = Method, y = -log10(P_Value), fill = Method)) + geom_boxplot() + theme_bw() +
  labs(x = paste("Method (Type I Error)",sep=""),
       y = "-log10(P-Value)",
       title = "P-value comparison across methods") +
  geom_hline(yintercept = -log10(0.05),color = "gray", lty = 2) +
  scale_x_discrete(labels = c("CSKAT" = paste("CSKAT (", round(cskat_val,4),")",sep = ""),
                              "LMM" = paste("LMM (", round(lmm_val,4), ")",sep = ""),
                              "LMM_Squared" = paste("LMM_Squared (", round(squared_val,4), ")",sep = ""),
                              "TwoStage" = paste("Two Stage (",round(two_val,4),")",sep = ""),
                              "PermANOVA" = paste("PermANOVA (",round(perm_val, 4),")",sep = ""))) +
  theme(axis.text.x = element_text(angle = 90))
ggsave(filename = paste("C:\\Users\\nweav\\OneDrive\\Documents\\Aim 1\\Most_Recent_Server\\large_OTU_plots\\error\\Comparison_boxplot_large_OTU_power.jpg",sep=""),
       width = 7, height = 7)




