require(ggplot2)
library(ggplot2)
library(plyr)
library(dplyr)
library(tidyverse)
require(DEoptim)
####################################
########## SERVER SLURM ############
####################################
# grid of 40 cells to be parallelized in slurm script - for running on the mpi cluster
grid <- expand.grid(sigma=seq(1,15,2), lambda=seq(0,8,2))
# which cell to pick via command line
cell_id <- commandArgs(trailingOnly = TRUE)
cell <- grid[cell_id,]


####################################
############ SIMULATION ############
####################################
bandit_simulation <- function(threshold, sigma, lambda){
  # NOTE we use exponential scale for the thresholds
  threshold<-exp(threshold)
  # setting parameters
  nparticipants <- 500  #number of participants per condition
  nblocks <- 20 #number of blocks per participant
  ntrials <- 100 #number of trials per block
  
  ############################
  ##### the slot machine #####
  ############################
  generate_task <- function(sigma) {
    # input: 
    ## sigma = standard deviation of the TRUE underlying reward distributions (identical for arm 1 and arm 2)
    # output: a function f that acts as a 2-armed bandit which:
    ## takes in the arm of choice (1 or 2) and
    ## returns a reward according to the TRUE underlying reward distribution
    
    # means of the TRUE underlying reward distribution (sampled from a hard-coded Gaussian)
    mu <- rnorm(2, 0, sqrt(100)) #runif(2, 30, 70)
    # standard deviations of the TRUE underlying reward distribution (passed on from function input)
    sd <- c(sigma, sigma) #sample(100:100, 2, replace=FALSE)
    # function f to be returned
    f <- function(arm) {
      reward = rnorm(1, mu[arm], sd[arm])
      return(reward)
    }
    return(f)
  }
  
  ###########################################################
  #### updating posteriers using a bayesian mean tracker ####
  ###########################################################
  bmt <- function(arm, reward, prevPost) {
    # inputs: 
    ## arm = arm of choice (1 or 2) in a given trial
    ## reward = value of reward obtained in this trial
    ## prevPost = prosterior from the previous trial
    # output: 
    ## posterior = updated posterior (participant's ESTIMATE of the underlying reward distributions)
    ##             i.e. an updated 2*2 matrix/data frame that looks like:
    ##         mu     var
    ## arm=1   mu[1]  var[1]
    ## arm=2   mu[2]  var[2]
    
    # priors
    mu0  <- 0 #prior mean
    var0 <- 100 #prior variance
    vare <- 10  #error varriance - fixed
    # update posteriors
    if (is.null(prevPost)) { #if there's no prior posterior, then we have the very first observation
      posterior <- data.frame(mu=rep(mu0,2), var=rep(var0,2))
    }
    else{ #if previous posterior exists, pass it down so we can update it
      posterior <- prevPost
    }
    # calculating Kalman gain i.e. the learning rate
    kGain <- posterior$var[arm] / (posterior$var[arm] + vare)
    # updating mean of chosen arm
    posterior$mu[arm]  <- posterior$mu[arm] + (kGain * (reward-posterior$mu[arm]))
    # updating variance of chosen arm
    posterior$var[arm] <- posterior$var[arm] * (1 - kGain)
    return(posterior)
  }
  
  ###################################################################################################################
  ##### accumulating evidence for decision using a symmetric random walk model i.e. SINGLE accumulation process #####
  ###################################################################################################################
  pondering <- function (posterior, threshold) {
    #inputs:
    ## posterior = updated posterior (participant's ESTIMATE of the underlying reward distributions) from function bmt()
    ## threshold = reflecting how much pondering needs to be done before a decision can be reached (hard-coded)
    # output:
    ## arm, i.e. the arm (1 or 2) deciced in this trial
    
    # initial evidence is set to 0 because it is a SYMMETRIC random walk
    evidence <- 0
    # updating accumulated evidence by repeated internal sampling from posteriors of two arms
    while ((evidence <= threshold) & (evidence >= -threshold)) {
      evidence <- (evidence + rnorm(1, posterior$mu[1], sqrt(posterior$var[1])) -
                     rnorm(1, posterior$mu[2], sqrt(posterior$var[2])) )
    }
    # choosing arm 1 or 2 based on the sign of the accumulated evidence at the end of the random walk 
    arm <- ifelse(evidence>0, 1, 2)
    return(arm)
  }
  
  
  #########################
  ##### MAIN FUNCTION #####
  #########################
  # initiate a data frame to save reward values from every trial
  reward_collect <- expand.grid(trial=1:ntrials, 
                                block=1:nblocks, 
                                participant=1:nparticipants, 
                                reward=0)
  
  # SIMULATION looping over all participants
  for (nsim in 1:nparticipants) { 
    # looping over blocks within each participant
    for (nblock in 1:nblocks) {                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        ##wihtin each one task/block, what happens in the 1st trial?
      arm <- sample(1:2, 1) # for the 1st trial, chose arm 1 or 2 at 50/50
      task <- generate_task(sigma)
      reward <- task(arm) # obtain reward from the 1st trial
      posterior <- bmt(arm, reward, prevPost = NULL) # update posteriors at the end of the 1st trial
      
      # looping over trials from the 2nd trial onwards:
      for (ntrial in 2:ntrials) {
        armnew <- pondering(posterior, threshold)
        rewardnew <- task(armnew)
        posterior <- bmt(armnew, rewardnew, prevPost = posterior)
        # record reward values from each trial
        reward <- c(reward, rewardnew)
      }
      # collect all reward values from this block
      reward_collect$reward[ ((nsim-1)*nblocks*ntrials + (nblock-1)*ntrials + 1) : ((nsim-1)*nblocks*ntrials + nblock*ntrials) ] <- reward
    }
    #print(paste("Simulation No.", nsim))
  }
  
  # CODE FROM TOBI: Mean(SDs from each block)
  # compute reward mean and variance over trials within each block
  reward_summary <- reward_collect %>%
    group_by(block) %>%
    summarise(mean_over_trials = mean(reward),
              sd_over_trials   = sd(reward))
  # now we summarise the summary ;) by averaging over participants and blocks
  reward_mean <- mean(reward_summary$mean_over_trials)
  reward_sd <- mean(reward_summary$sd_over_trials)
  # OLD CODE
  #reward_mean <- mean(reward_collect$reward)
  #reward_sd <- sd(reward_collect$reward)
  
  utility <- reward_mean - lambda*reward_sd
  return(-utility)
} # the end of bandit_simulation(threshold, sigma, lambda)



########################################################################
########## Global Optimization by Differential Evolution ###############
########################################################################
# set upper and lower bound for the range of threshold to be searched through
threshold_lbound <- -5
threshold_ubound <-  6

# DEoptim()
fit <- DEoptim(bandit_simulation, lower=threshold_lbound, upper=threshold_ubound,
               sigma=cell$sigma, lambda=cell$lambda,
               DEoptim.control(itermax=50))
best_t <- fit$optim$bestmem
result_collect <- data.frame(cell$sigma, cell$lambda, best_t)
names(result_collect) <- c("sigma", "lambda", "best_t")

print(result_collect)
write.csv(result_collect, paste0("result_", result_collect$sigma, "_", result_collect$lambda, "_", cell_id,  ".csv"))














