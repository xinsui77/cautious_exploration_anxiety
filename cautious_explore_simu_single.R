require(ggplot2)
library(ggplot2)
library(plyr)
library(dplyr)

############################
##### the slot machine #####
############################
generate_param <- function() {
  # means of the TRUE underlying reward distribution
  mu <- rnorm(2, 0, sqrt(100)) #same as Gershman18 Experiment 1
  # standard deviations of the TRUE underlying reward distribution
  var <- runif(2, min=5, max=15) #same as Gershman18 Experiment 1
  param <- data.frame(mu, var)
  
  return(param)
}

generate_task <- function(param) {
  ### input: 
  # mu: 1*2 vector containing true means for arm 1 and 2
  # sd: 1*2 vector containing true standard deviations for arm 1 and 2
  ### output: a function f that acts as a 2-armed bandit which:
  # takes in the arm of choice (1 or 2) and
  # returns a reward according to the TRUE underlying reward distribution
  
  f <- function(arm) {
    reward = rnorm(1, param$mu[arm], sqrt(param$var[arm]))
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
  mu0  <- 0 #prior means
  var0 <- 100 #prior variances
  vare <- c(10,10)  #error varriances
  
  # update posteriors
  if (is.null(prevPost)) { #if there's no prior posterior, then we have the very first observation
    posterior <- data.frame(mu=rep(mu0,2), var=rep(var0,2))
  }
  else { #if previous posterior exists, pass it down so we can update it
    posterior <- prevPost
  }
  
  # calculating Kalman gain i.e. the learning rate
  kGain <- posterior$var[arm] / (posterior$var[arm] + vare[arm])
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
                            rnorm(1, posterior$mu[2], sqrt(posterior$var[2])))
  }
  # choosing arm 1 or 2 based on the sign of the accumulated evidence at the end of the random walk 
  arm <- ifelse(evidence>0, 1, 2)
  return(arm)
}





##################################################
################ MAIN FUNCTION ###################
##################################################

# thresholds
threshold_low <- 0
threshold_high <- 400
threshold_step <- 50
nthresholds = (threshold_high - threshold_low) /threshold_step + 1 

# experiment parameters
nparticipants <- 500  #number of participants per condition   ###Gershman18 = 500
nblocks <- 20 #number of blocks per participant
ntrials <- 100 #number of trials per block
nsims <- nparticipants * nthresholds #==nrow(dcollect)

# experiment grid, one row == one sim
dcollect <- expand.grid(participant=1:nparticipants, 
                        threshold=seq(threshold_low, threshold_high, threshold_step))

# initiate a data frame for saving true generative parameters for each block
block_collect <- expand.grid(block=1:nblocks,
                             nsim=1:nsims,
                             true_mu1=0,
                             true_mu2=0,
                             true_var1=0,
                             true_var2=0)
# initiate a data frame for saving arm choices, rewards, and posteriors from every trial
trial_collect <- expand.grid(trial=1:ntrials,
                             block=1:nblocks,
                             nsim=1:nsims,
                             arm=0,
                             reward=0,
                             mu1=0,
                             mu2=0,
                             var1=0,
                             var2=0)

# SIMULATION looping over all sims
for (nsim in 1:nsims) {
  # obtaining and saving sim-level parameters from dcollect
  threshold <- dcollect$threshold[nsim]
  
  # looping over blocks within each participant
  for (nblock in 1:nblocks) {
    # generate and save block-level parameters
    parameters <- generate_param()
    task <- generate_task(parameters)
    
    block_collect$true_mu1[(nsim-1)*nblocks+nblock] = parameters$mu[1]
    block_collect$true_mu2[(nsim-1)*nblocks+nblock] = parameters$mu[2]
    block_collect$true_var1[(nsim-1)*nblocks+nblock] = parameters$var[1]
    block_collect$true_var2[(nsim-1)*nblocks+nblock] = parameters$var[2]
    
    ## 1st trial
    arm <- sample(1:2, 1) # for the 1st trial, chose arm 1 or 2 at 50/50
    reward <- task(arm) # obtain reward from the 1st trial
    posterior <- bmt(arm, reward, prevPost = NULL) # update posteriors at the end of the 1st trial
    
    mu1 <- posterior$mu[1]
    mu2 <- posterior$mu[2]
    var1 <- posterior$var[1]
    var2 <- posterior$var[2]

    # from the 2nd trial onwards:
    for (ntrial in 2:ntrials) {
      armnew <- pondering(posterior, threshold)
      rewardnew <- task(armnew)
      posterior <- bmt(armnew, rewardnew, prevPost=posterior)
      
      mu1new <- posterior$mu[1]
      mu2new <- posterior$mu[2]
      var1new <- posterior$var[1]
      var2new <- posterior$var[2]
      
      # save trial-level values
      arm <- c(arm, armnew)
      reward <- c(reward, rewardnew)
      mu1 <- c(mu1, mu1new)
      mu2 <- c(mu2, mu2new)
      var1 <- c(var1, var1new)
      var2 <- c(var2, var2new)
      } # one block finished 
    
    # save all values from this block
    trial_collect$arm[((nsim-1)*nblocks*ntrials + (nblock-1)*ntrials + 1) : ((nsim-1)*nblocks*ntrials + nblock*ntrials) ] <- arm
    trial_collect$reward[((nsim-1)*nblocks*ntrials + (nblock-1)*ntrials + 1) : ((nsim-1)*nblocks*ntrials + nblock*ntrials) ] <- reward
    trial_collect$mu1[((nsim-1)*nblocks*ntrials + (nblock-1)*ntrials + 1) : ((nsim-1)*nblocks*ntrials + nblock*ntrials) ] <- mu1
    trial_collect$mu2[((nsim-1)*nblocks*ntrials + (nblock-1)*ntrials + 1) : ((nsim-1)*nblocks*ntrials + nblock*ntrials) ] <- mu2
    trial_collect$var1[((nsim-1)*nblocks*ntrials + (nblock-1)*ntrials + 1) : ((nsim-1)*nblocks*ntrials + nblock*ntrials) ] <- var1
    trial_collect$var2[((nsim-1)*nblocks*ntrials + (nblock-1)*ntrials + 1) : ((nsim-1)*nblocks*ntrials + nblock*ntrials) ] <- var2
    } # one sim finished
  print(paste("Sim No.", nsim, " out of ", nsims, ", Threshold = ", threshold))
}


write.csv(dcollect, "collect_sim.csv")
write.csv(block_collect, "collect_block.csv")
write.csv(trial_collect, "collect_trial.csv")