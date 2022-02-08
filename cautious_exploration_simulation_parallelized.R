require(ggplot2)
library(ggplot2)
library(plyr)
library(dplyr)
require(DEoptim)
####################################
########## SERVER SLURM ############
####################################
# grid of cells to be parallelized in slurm script
grid <- expand.grid(sigma=seq(1,20,1), lambda=seq(0,10,0.5))

# which cell to pick via command line
cell_id <- commandArgs(trailingOnly = TRUE)
cell <- grid[cell_id,]

####################################
############ SIMULATION ############
####################################
bandit_simulation <- function(threshold, sigma, lambda){
  threshold<-exp(threshold)
  ## setting parameters
  nparticipants <- 500  #number of participants per threshold per variance
  nblocks <- 20 #number of blocks
  ntrials <- 100 #number of trials per block
  ############################
  ##### the slot machine #####
  ############################
  generate_task <- function(sigma) {
    #input: sigma, the SD of the TRUE reward distribution
    #output: a function that acts as a two-armed bandit which:
    ##takes in the arm (1 or 2) and
    ##returns a reward according to the TRUE distribution
    #mu's are the means of the TRUE reward distribution (hard-coded, fixed but sampled)
    #sigma's are the standard deviations of the TRUE reward distribution
    mu <- rnorm(2, 0, sqrt(100)) #runif(2, 30, 70)
    sd <- c(sigma, sigma) #sample(100:100, 2, replace=FALSE)
    f <- function(arm) {
      reward = rnorm(1, mu[arm], sd[arm])
      return(reward)
    }
    return(f)
  }
  ##########################################################
  #### updating estimates using a bayesian mean tracker ####
  ##########################################################
  bmt <- function(arm, reward, prevPost) {
    #inputs: 
    ##arm, i.e. arm (1 or 2) of choice in this specific trial
    ##reward obtained in this specific trial
    ##prevPost, i.e. prosterior from the previous trial
    #output: 
    ##posterior, i.e. updated posterior (ESTIMATE of the reward distributions)
    ##i.e. updated 2*2 matrix for reward distributions, which looks like:
    ##         mu     var
    ## arm=1   mu[1]  var[1]
    ## arm=2   mu[2]  var[2]
    mu0  <- 0 #prior mean
    var0 <- 100 #prior variance
    vare <- 10  #error varriance
    if (is.null(prevPost)) { 
      #if no prior posterior, then it is the very first observation
      posterior <- data.frame(mu=rep(mu0,2), var=rep(var0,2))
    }
    else{ 
      #if previous posterior exists, update it
      posterior <- prevPost
    }
    #Kalman gain
    kGain <- posterior$var[arm] / (posterior$var[arm] + vare)
    #updating mean
    posterior$mu[arm]  <- posterior$mu[arm] + (kGain * (reward-posterior$mu[arm]))
    #updating variance
    posterior$var[arm] <- posterior$var[arm] * (1 - kGain)
    return(posterior)
  }
  
  ############################################################################################
  ##### accumulating evidence using a random walk model i.e. SINGLE accumulation process #####
  ############################################################################################
  pondering <- function (posterior, threshold) {
    #inputs:
    ##posterior, i.e. the updated posterior (ESTIMATE of reward distributions) from function "bmt()"
    ##threshold (hard-coded)
    #output:
    ##arm, i.e. the arm (1 or 2) chosen/deciced for this specific trial
    evidence <- 0
    while ((evidence <= threshold) & (evidence >= -threshold)) {
      evidence <- (evidence + rnorm(1, posterior$mu[1], sqrt(posterior$var[1])) -
                              rnorm(1, posterior$mu[2], sqrt(posterior$var[2])) )
    }
    arm <- ifelse(evidence>0, 1, 2)
    return(arm)
  }
  
  ##############################################################################################
  #### accumulating evidence using a race model i.e. two independent accumulation processes ####
  ##############################################################################################
  # pondering <- function (posterior, threshold) {
  #   #drift difussion model for the evidence accumulation process when making each decision on arm
  #   #inputs:
  #   ##posterior, i.e. the updated posterior (ESTIMATE of reward distributions) from function "bmt()"
  #   ##threshold (hard-coded)
  #   #output:
  #   ##arm, i.e. the arm (1 or 2) chosen/deciced for this specific trial
  #   evidence1 <- 0
  #   evidence2 <- 0
  #   while ((evidence1<=threshold)&(evidence2<=threshold)&(evidence1>=-threshold)&(evidence2>=-threshold)) {
  #     evidence1 <- evidence1 + rnorm(1, posterior$mu[1], sqrt(posterior$var[1]))
  #     evidence2 <- evidence2 + rnorm(1, posterior$mu[2], sqrt(posterior$var[2]))
  #     #print(paste("Accumulated evidence = ", evidence, "."))
  #   }
  #   arm <- ifelse(evidence1>evidence2, 1, 2)
  #   #cat("Arm of choice = ", arm, ".\n", sep = " ")
  #   return(arm)
  # }
  
  #########################
  ##### MAIN FUNCTION #####
  #########################
  reward_collect <- expand.grid(trial=1:ntrials, 
                                block=1:nblocks, 
                                participant=1:nparticipants, 
                                reward=0)
  #dcollect <- expand.grid(participant=1:nparticipants, mean=0)
  #learning_collect <- expand.grid(trial=1:ntrials, participant=1:nparticipants, reward_mean=0)
  ###looping over simulation (i.e. every row of dcollect = nparticipants * nthresholds)
  for (nsim in 1:nparticipants) { 
    #initialize matrices for each simulation
    #reward_matrix <- matrix(0, nrow=ntrials, ncol=nblocks)
    #looping over blocks withom each simulation
    for (nblock in 1:nblocks) {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ##wihtin each one task/block, what happens in the 1st trial?
      task <- generate_task(sigma)
      arm <- sample(1:2, 1) #for the 1st trial, chosing an arm completely randomly 
      reward <- task(arm) #reward from the 1st trial
      posterior <- bmt(arm, reward, prevPost = NULL) #updated posterior from the 1st trial
      
      ## every trial from 2nd trial onwards:
      for (ntrial in 2:ntrials) {
        armnew <- pondering(posterior, threshold)
        rewardnew <- task(armnew)
        posterior <- bmt(armnew, rewardnew, prevPost = posterior)
        
        # tracking rewards
        reward <- c(reward, rewardnew)
      } #finished all trials within the block
      
      reward_collect$reward[((nsim-1)*nblocks*ntrials+(nblock-1)*ntrials+1):((nsim-1)*nblocks*ntrials+nblock*ntrials)] <- reward
      #reward_matrix[, nblock] <- reward
    }  # finished all block within the participant / sim
    
    #calculating participant's mean (mean rewards averaging over all blocks in the simulation) & min (average rewards from the simulation's worst block) 
    #dcollect$mean[nsim] <- mean(apply(reward_matrix, 2, mean))
    #dcollect$min[nsim]  <- min(apply(reward_matrix, 2, mean))
    
    #storing average rewards per trial for plotting learning curve 
    #learning_collect$reward_mean[((nsim-1)*ntrials+1):(nsim*ntrials)] = apply(reward_matrix, 1, mean)
    #learning_collect$reward_min[((nsim-1)*ntrials+1):(nsim*ntrials)] = apply(reward_matrix, 1, min)
    
    #print(paste("Simulation No.", nsim))
  } # finished all participants / sim's
  
  # calculate utility
  reward_mean <- mean(reward_collect$reward)
  reward_sd <- sd(reward_collect$reward)
  utility <- reward_mean - lambda*reward_sd
  return(-utility)
} # end of bandit_simulation(threshold, sigma, lambda)

##### checking how utility outputs vary across repetitions
# utility <- rep(0,10)
# for(i in 1:10){
#   utility[i] <- bandit_simulation(0,10,0.1)
#   print(i)
# }
# utility

########################################################################
########## Global Optimization by Differential Evolution ###############
########################################################################

## DEoptim
threshold_lbound<- -5
threshold_ubound<-  6

## experiment matrix
#dcollect <- expand.grid(sigma=seq(1,15,2), lambda=seq(0,10,1), best_t=0)
#fit <- rep(0, nrow(dcollect))

fit <- DEoptim(bandit_simulation, lower=threshold_lbound, upper=threshold_ubound,
               sigma=cell$sigma, lambda=cell$lambda,
               DEoptim.control(itermax=50))
best_t<-fit$optim$bestmem
result_collect<-data.frame(cell$sigma, cell$lambda, best_t)
names(result_collect) <- c("sigma", "lambda", "best_t")

print(result_collect)
write.csv(result_collect, paste0("result_", result_collect$sigma, "_", result_collect$lambda, "_", cell_id,  ".csv"))

#ggplot(dcollect, aes(x=lambda, y=sigma, fill=best_t)) + geom_tile()

#### nSIMs = 1 ONLY for time estimation
# fit_test <- DEoptim(bandit_simulation, lower=threshold_lbound, upper=threshold_ubound, 
#                sigma=1, lambda=1, DEoptim.control(itermax=1))
# best_t_test <- fit_test$optim$bestmem

######## ERIC's code
# my_function_to_optimize<-function(threshold, sigma, lambda){
#   threshold<-exp(threshold)
#   #create bandit with sigma
#   #run 500 participants with threshold over your trials and blocks
#   #calulcate utility<-mean(rewards)+lambda*sd(rewards)
#   return(-utility)
# }