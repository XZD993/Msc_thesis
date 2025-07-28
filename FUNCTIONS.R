# Msc thesis:
# A group sequential three-arm two-stage trial design with rank-based endpoints
####### All necessary functions and packages for simulations ###################
################################################################################
#
#
#--------------------------- Package loading -----------------------------------
#
library(parallel)
library(tidyverse)
library(pbapply)  # faster parallel + nice progress bar
library(dplyr)
library(Rcpp)
library(microbenchmark)
#
#
# ----------------------------- Functions --------------------------------------
#
################################################################################
############################ Helper functions ##################################
################################################################################
#
### 1. Customized pairwise rank calculator (C++ version)
#
cppFunction('
NumericVector rank_func_cpp(NumericMatrix data) {
  int n = data.nrow();
  int p = data.ncol();
  NumericVector rr(n, 0.0);

  for (int i = 0; i < n - 1; ++i) {
    for (int j = i + 1; j < n; ++j) {
      double result = 0.5;

      for (int k = 0; k < p; ++k) {
        double Ai = data(i, k);
        double Bi = data(j, k);

        if (k == 0) { // survival time with 30-day cutoff
          if (Ai >= 30 && Bi >= 30) continue;
          else {
            if (Ai > Bi) { result = 1; break; }
            if (Ai < Bi) { result = 0; break; }
          }
        } else {
          if (Ai > Bi) { result = 1; break; }
          if (Ai < Bi) { result = 0; break; }
        }
      }

      rr[i] += result;
      rr[j] += (1.0 - result);
    }
  }

  for (int i = 0; i < n; ++i) {
    rr[i] += 1;
  }

  return rr;
}
')
#
#
#
#
################################################################################
################# Functions for simulating each design #########################
################################################################################
#
#
### ----------- 1. Data generator - Hierarchical endpoints
## The components of hierarchical endpoints:
# Component A. survival time: if 2) or 3), compare survival days, longer one is the winner. if 1), continue to compare the next component
# 1) both live longer than 30 days;
# 2) one live longer than 30 days and another shorter than 30 days;
# 3) both live shorter than 30 days.
#
# Component B. Kidney recovery: 0 - no recovery; 1 - recovery
dat_gen_multi <- function(nn=50, 
                          lambda=0.03, 
                          p = 0.3,
                          arm=1)
{
  # component 1: Survival data
  surv_time <- rexp(nn, rate = lambda)
  # component 2: Kidney recovery (Binary data)
  qq <- rbinom(nn, 1, p)
  
  data.frame(survtime = surv_time, recovery = qq, arm.inf = arm)
}
#
#
### ------- 2. WMW-Test statistics (Wij) calculator (for Hierarchical endpoint)
Wij_cal_multi <- function(dd.0=dat.arm.0,dd.1=dat.arm.1)
{
  # Note: dd.0 denotes control group and dd.1 denotes experimental group
  # pairwise rank
  dd.ij<-rbind(dd.0,dd.1)
  dd.ij$rank.all <- rank_func_cpp(as.matrix(dd.ij[,1:2]))
  arms<-unique(dd.ij$arm.inf)
  rankij_i <- dd.ij$rank.all[dd.ij$arm.inf == arms[1]] # pairwise rank of control group
  rankij_j <- dd.ij$rank.all[dd.ij$arm.inf == arms[2]] # pairwise rank of experimental group
  #
  # W.ij calculation
  #
  N <- nrow(dd.ij)
  ni <- length(rankij_i)
  nj <- length(rankij_j)
  var.ij <- sum((dd.ij$rank.all - (N+1)/2)^2)/(N-1)
  W.ij <- (mean(rankij_j) - mean(rankij_i)) * sqrt(ni*nj/N) / sqrt(var.ij)
  #
  # return(c(W.ij, mean(rankij_j), mean(rankij_i)))
  return(W.ij)
}
#
#
### ---------3. Simulation function for each stopping rules --------------------
#
################################################################################
############## 1) Simultaneous stopping rule - referred as MJW #################
# If one arm shows efficacy, the whole trial is stopped, and the arm is recommended as the success arm
#
sim_rfc_MJW <- function(i=1, n=c(150, 150),
                        lambda.all=c(0.03, 0.03, 0.03),  # c(treatment1, treatment2, control)
                        p.all=c(0.3, 0.3, 0.3),               # c(treatment1, treatment2, control)
                        boundaries=c(0, 2, 2))
{
  #  on.exit(browser())
  # sample size of each arm in stage1 
  n1 <- n[1]/3 
  # 
  # endpoint component 1
  lambda1 = lambda.all[1]
  lambda2 = lambda.all[2]
  lambda0 = lambda.all[3]
  # endpoint component 2
  p1 = p.all[1]
  p2 = p.all[2]
  p0 = p.all[3]
  #
  dat1_stage1<-dat_gen_multi(nn=n1, lambda=lambda1, p=p1, arm=1) # arm1 (treatment1)
  dat2_stage1<-dat_gen_multi(nn=n1, lambda=lambda2, p=p2, arm=2) # arm2 (treatment2)
  dat0_stage1<-dat_gen_multi(nn=n1, lambda=lambda0, p=p0, arm=0) # control arm
  
  #
  W01<-Wij_cal_multi(dat0_stage1,dat1_stage1)
  W02<-Wij_cal_multi(dat0_stage1,dat2_stage1)
  #
  c0=boundaries[1] # futility boundary
  c1=boundaries[2] # efficacy boundary at stage 1
  c2=boundaries[3] # efficacy boundary at stage 2
  #
  # -------------------------Decision in stage 1------------------------------ #
  #
  #                                     W01
  #                       futility(c0)        effective(c1)
  #                             |                 |   
  #                     0       |        1        |       4
  #                             |                 |
  #    futility(c0)  ---------------------------------------------
  #                             |                 |
  # W02                 2       |        3        |       4
  #                             |                 |
  #    effective(c1) ---------------------------------------------
  #                             |                 |
  #                     4       |        4        |       4
  #                             |                 |
  stop <- 0       # indicator for the above situations:
  # 0: stop trial due to futility
  # 1: continue with arm 1 and drop arm 2 due to futility
  # 2: continue with arm 2 and drop arm 1 due to futility
  # 3: continue with both arms
  # 4: stop trial if at least one arm shows efficacy
  #
  effective <- 0  # binary: 0 means accepting global null hypothesis, and 1 vice versa,
  effective_stage1 <- 0 # binary: 0 means no arm shows efficacy in stage 1, and 1 vice versa,
  effective_stage2 <- 0 # binary: 0 means no arm shows efficacy in stage 2, and 1 vice versa,
  winner <- 0 
  # winner == 0: no treatment arm is selected;
  # winner == 1: treatment arm 1 is selected;
  # winner == 2: treatment arm 2 is selected;
  # winner == 3: both arms are selected;
  
  if (W01 > c1 || W02 > c1) {
    stop <- 4
    effective_stage1 <- 1
    if (W01 > c1 && W02 > c1) {
      winner <- 3
    }
    else if (W01 > c1) {
      winner <- 1
    }
    else if (W02 > c1) {
      winner <- 2
    }
  } 
  else if (W01 < c0 && W02 < c0) {
    stop <- 0
  } 
  else if (W01 >= c0 && W01 <= c1 && W02 < c0) {
    stop <- 1
  } 
  else if (W02 >= c0 && W02 <= c1 && W01 < c0) {
    stop <- 2
  } 
  else if (W01 >= c0 && W01 <= c1 && W02 >= c0 && W02 <= c1) {
    stop <- 3
  }
  #
  #
  # -------------------------Decision in stage 2------------------------------ #
  #
  W_new <- -99   # W statistics if one arm is dropped at last stage
  if (stop %in% 1:2) {
    n2 <- n[2]/2
    arm_keep <- stop
    
    dat0_stage2 <- rbind(dat0_stage1, dat_gen_multi(n2, lambda0, p0, 0))
    dat_keep_stage2 <- rbind(get(paste0("dat", arm_keep, "_stage1")), 
                             dat_gen_multi(n2, lambda.all[arm_keep], p.all[arm_keep], arm_keep))
    
    W_new <- Wij_cal_multi(dat0_stage2, dat_keep_stage2)
    if (W_new > c2) {
      effective_stage2 <- 1
      winner <- arm_keep
    }
  } 
  else if (stop == 3) {
    n2 <- n[2]/3
    dat1_stage2 <- rbind(dat1_stage1, dat_gen_multi(n2, lambda1, p1, 1))
    dat2_stage2 <- rbind(dat2_stage1, dat_gen_multi(n2, lambda2, p2, 2))
    dat0_stage2 <- rbind(dat0_stage1, dat_gen_multi(n2, lambda0, p0, 0))
    
    W01 <- Wij_cal_multi(dat0_stage2, dat1_stage2)
    W02 <- Wij_cal_multi(dat0_stage2, dat2_stage2)
    
    if (W01 > c2 || W02 > c2) {
      effective_stage2 <- 1
      if (W01 > c1 && W02 > c1) {
        winner <- 3
      }
      else if (W01 > c1) {
        winner <- 1
      }
      else if (W02 > c1) {
        winner <- 2
      }
    }
  }
  # success of the trial, i.e. global test result
  # If efficacy was shown in at least one stage, the trial is successful
  if (any(c(effective_stage1, effective_stage2) == 1)) {
    effective <- 1
  }
  res<-c(Stop = stop,
         W.01 = W01, W.02 = W02, W.new = W_new, Arm.select = winner,
         Success.stage1 = effective_stage1, Success.stage2 = effective_stage2, Success=effective)
  return(res)
}
#start.time <- Sys.time()
#test1 <- sim_rfc_MJW()
#end.time <- Sys.time()
#run.time <- difftime(end.time, start.time, units = "mins")
#
#
#
################################################################################
########### 2) Simultaneous stopping rule plus drop-the-loser idea, ############
###################### referred as MJW-DTL #####################################
# Similar to simultaneous design, but if both test statistics are in continuation range, drop the treatment arm with smaller W and continue with others
#
sim_rfc_MJW_DTL <- function(i=1, n=c(150, 150),
                            lambda.all=c(0.03, 0.03, 0.03),  # c(treatment1, treatment2, control)
                            p.all=c(0.3, 0.3, 0.3),               # c(treatment1, treatment2, control)
                            boundaries=c(0, 2, 2))
{
  #  on.exit(browser())
  # sample size of each arm in stage1 
  n1 <- n[1]/3 
  # 
  # endpoint component 1
  lambda1 = lambda.all[1]
  lambda2 = lambda.all[2]
  lambda0 = lambda.all[3]
  # endpoint component 2
  p1 = p.all[1]
  p2 = p.all[2]
  p0 = p.all[3]
  #
  dat1_stage1<-dat_gen_multi(nn=n1, lambda=lambda1, p=p1, arm=1) # arm1 (treatment1)
  dat2_stage1<-dat_gen_multi(nn=n1, lambda=lambda2, p=p2, arm=2) # arm2 (treatment2)
  dat0_stage1<-dat_gen_multi(nn=n1, lambda=lambda0, p=p0, arm=0) # control arm
  
  #
  W01<-Wij_cal_multi(dat0_stage1,dat1_stage1)
  W02<-Wij_cal_multi(dat0_stage1,dat2_stage1)
  #
  c0=boundaries[1] # futility boundary
  c1=boundaries[2] # efficacy boundary at stage 1
  c2=boundaries[3] # efficacy boundary at stage 2
  #
  #
  # -------------------------Decision in stage 1------------------------------ #
  #
  #                                     W01
  #                       futility(c0)        effective(c1)
  #                             |                 |   
  #                     0       |        1        |       4
  #                             |                 |
  #    futility(c0)  ---------------------------------------------
  #                             |                 |
  # W02                 2       |     1/2/3       |       4
  #                             |                 |
  #    effective(c1) ---------------------------------------------
  #                             |                 |
  #                     4       |        4        |       4
  #                             |                 |
  stop <- 0       # indicator for the above situations:
  # 0: stop trial due to futility
  # 1: continue with arm 1 and drop arm 2 due to futility
  # 2: continue with arm 2 and drop arm 1 due to futility
  # 3: rare case. If two test statistics are the same, continue with both arm
  # 4: stop trial if at least one arm shows efficacy
  #
  effective <- 0  # binary: 0 means accepting global null hypothesis, and 1 vice versa,
  effective_stage1 <- 0 # binary: 0 means no arm shows efficacy in stage 1, and 1 vice versa,
  effective_stage2 <- 0 # binary: 0 means no arm shows efficacy in stage 2, and 1 vice versa,
  winner <- 0 
  # winner == 0: no treatment arm is selected;
  # winner == 1: treatment arm 1 is selected;
  # winner == 2: treatment arm 2 is selected;
  # winner == 3: both arms are selected;
  
  if (W01 > c1 || W02 > c1) {
    stop <- 4
    effective_stage1 <- 1
    if (W01 > c1 && W02 > c1) {
      winner <- 3
    }
    else if (W01 > c1) {
      winner <- 1
    }
    else if (W02 > c1) {
      winner <- 2
    }
  } 
  else if (W01 < c0 && W02 < c0) {
    stop <- 0
  } 
  else if (W01 >= c0 && W01 <= c1 && W02 < c0) {
    stop <- 1
  } 
  else if (W02 >= c0 && W02 <= c1 && W01 < c0) {
    stop <- 2
  } 
  else if (W01 >= c0 && W01 <= c1 && W02 >= c0 && W02 <= c1) {
    stop <- ifelse(W01 > W02, 1, ifelse(W01 < W02, 2, 3))
  }
  #
  #
  # -------------------------Decision in stage 2------------------------------ #
  #
  W_new <- -99   # W statistics if one arm is dropped at last stage
  if (stop %in% 1:2) {
    n2 <- n[2]/2
    arm_keep <- stop
    
    dat0_stage2 <- rbind(dat0_stage1, dat_gen_multi(n2, lambda0, p0, 0))
    dat_keep_stage2 <- rbind(get(paste0("dat", arm_keep, "_stage1")), 
                             dat_gen_multi(n2, lambda.all[arm_keep], p.all[arm_keep], arm_keep))
    
    W_new <- Wij_cal_multi(dat0_stage2, dat_keep_stage2)
    if (W_new > c2) {
      effective_stage2 <- 1
      winner <- arm_keep
    }
  } 
  else if (stop == 3) {
    n2 <- n[2]/3
    dat1_stage2 <- rbind(dat1_stage1, dat_gen_multi(n2, lambda1, p1, 1))
    dat2_stage2 <- rbind(dat2_stage1, dat_gen_multi(n2, lambda2, p2, 2))
    dat0_stage2 <- rbind(dat0_stage1, dat_gen_multi(n2, lambda0, p0, 0))
    
    W01 <- Wij_cal_multi(dat0_stage2, dat1_stage2)
    W02 <- Wij_cal_multi(dat0_stage2, dat2_stage2)
    
    if (W01 > c2 || W02 > c2) {
      effective_stage2 <- 1
      if (W01 > c1 && W02 > c1) {
        winner <- 3
      }
      else if (W01 > c1) {
        winner <- 1
      }
      else if (W02 > c1) {
        winner <- 2
      }
    }
  }
  # success of the trial, i.e. global test result
  # If efficacy was shown in at least one stage, the trial is successful
  if (any(c(effective_stage1, effective_stage2) == 1)) {
    effective <- 1
  }
  res<-c(Stop = stop,
         W.01 = W01, W.02 = W02, W.new = W_new, Arm.select = winner,
         Success.stage1 = effective_stage1, Success.stage2 = effective_stage2, Success=effective)
  return(res)
}
#
#test2 <- sim_rfc_MJW_DTL()
#
#
#
################################################################################
############ 3) Arm-specific stopping rule - referred as AS ####################
# drop treatment arm(s) due to both futility and efficacy, then continue with the others (if there is any) and control arm
#
sim_rfc_AS <- function(i=1, n=c(150, 150),
                       lambda.all=c(0.03, 0.03, 0.03),  # c(treatment1, treatment2, control)
                       p.all=c(0.3, 0.3, 0.3),               # c(treatment1, treatment2, control)
                       boundaries=c(0, 2, 2))
{
  #  on.exit(browser())
  # sample size of each arm in stage1 
  n1 <- n[1]/3 
  # 
  # endpoint component 1
  lambda1 = lambda.all[1]
  lambda2 = lambda.all[2]
  lambda0 = lambda.all[3]
  # endpoint component 2
  p1 = p.all[1]
  p2 = p.all[2]
  p0 = p.all[3]
  #
  dat1_stage1<-dat_gen_multi(nn=n1, lambda=lambda1, p=p1, arm=1) # arm1 (treatment1)
  dat2_stage1<-dat_gen_multi(nn=n1, lambda=lambda2, p=p2, arm=2) # arm2 (treatment2)
  dat0_stage1<-dat_gen_multi(nn=n1, lambda=lambda0, p=p0, arm=0) # control arm
  
  #
  W01<-Wij_cal_multi(dat0_stage1,dat1_stage1)
  W02<-Wij_cal_multi(dat0_stage1,dat2_stage1)
  #
  c0=boundaries[1] # futility boundary
  c1=boundaries[2] # efficacy boundary at stage 1
  c2=boundaries[3] # efficacy boundary at stage 2
  #
  # -------------------------Decision in stage 1-------------------------------#
  # 
  #                                     W01
  #                       futility(c0)        effective(c1)
  #                             |                 |   
  #                     0       |        1        |       4
  #                             |                 |
  #    futility(c0)  --------------------------------------------
  #                             |                 |
  # W02                 2       |        3        |       2
  #                             |                 |
  #    effective(c1) --------------------------------------------
  #                             |                 |
  #                     4       |        1        |       4
  #                             |                 |
  stop <- 0       # indicator for the above situations:
  # 0: stop trial due to futility shown from all treatments
  # 1: continue with arm 1 and drop arm 2 due to futility or efficacy
  # 2: continue with arm 2 and drop arm 1 due to futility or efficacy
  # 3: continue with both arm
  # 4: stop trial because all treatments are dropped due to either futility or efficacy
  #
  effective <- 0
  effective_stage1 <- 0
  effective_stage2 <- 0
  winner <- 0
  # winner == 0: no treatment arm is selected;
  # winner == 1: treatment arm 1 is selected;
  # winner == 2: treatment arm 2 is selected;
  # winner == 3: both arms are selected;
  
  if (W01 > c1 || W02 > c1) {
    effective_stage1 <- 1
    if ((W01 > c1) && (W02 >= c0 && W02 <= c1)){
      stop <- 2
      winner <- 1
    }
    else if ((W02 > c1) && (W01 >= c0 && W01 <= c1)){
      stop <- 1
      winner <- 2
    }
    else {
      stop <- 4
      winner <- ifelse(all(c(W01, W02) > c1), 3, 
                       ifelse(W01 > c1, 1, 2))
    }
  } 
  else if (W01 < c0 && W02 < c0) {
    stop <- 0
  } 
  else if (W01 >= c0 && W01 <= c1 && W02 < c0) {
    stop <- 1
  } 
  else if (W02 >= c0 && W02 <= c1 && W01 < c0) {
    stop <- 2
  } 
  else if (W01 >= c0 && W01 <= c1 && W02 >= c0 && W02 <= c1) {
    stop <- 3
  }
  #
  #
  # -------------------------Decision in stage 2------------------------------ #
  #
  W_new <- -99
  if (stop %in% 1:2) {
    n2 <- n[2]/2
    arm_keep <- stop
    
    dat0_stage2 <- rbind(dat0_stage1, dat_gen_multi(n2, lambda0, p0, 0))
    dat_keep_stage2 <- rbind(get(paste0("dat", arm_keep, "_stage1")), 
                             dat_gen_multi(n2, lambda.all[arm_keep], p.all[arm_keep], arm_keep))
    
    W_new <- Wij_cal_multi(dat0_stage2, dat_keep_stage2)
    if (W_new > c2) {
      effective_stage2 <- 1
      winner <- ifelse(winner != 0, 3, arm_keep)
    }
  } 
  else if (stop == 3) {
    n2 <- n[2]/3
    dat1_stage2 <- rbind(dat1_stage1, dat_gen_multi(n2, lambda1, p1, 1))
    dat2_stage2 <- rbind(dat2_stage1, dat_gen_multi(n2, lambda2, p2, 2))
    dat0_stage2 <- rbind(dat0_stage1, dat_gen_multi(n2, lambda0, p0, 0))
    
    W01 <- Wij_cal_multi(dat0_stage2, dat1_stage2)
    W02 <- Wij_cal_multi(dat0_stage2, dat2_stage2)
    
    if (W01 > c2 || W02 > c2) {
      effective_stage2 <- 1
      if (W01 > c1 && W02 > c1) {
        winner <- 3
      }
      else if (W01 > c1) {
        winner <- 1
      }
      else if (W02 > c1) {
        winner <- 2
      }
    }
  }
  # success of the trial, i.e. global test result
  # If efficacy was shown in at least one stage, the trial is successful
  if (any(c(effective_stage1, effective_stage2) == 1)) {
    effective <- 1
  }
  res<-c(Stop = stop,
         W.01 = W01, W.02 = W02, W.new = W_new, Arm.select = winner,
         Success.stage1 = effective_stage1, Success.stage2 = effective_stage2, Success=effective)
  return(res)
}

#test3 <- sim_rfc_AS()
#
#
#
#
################################################################################
########### 4) Arm-specific stopping rule plus drop-the-loser idea #############
############################ referred as AS-DTL ################################
# similar to AS design, but drop the treatment arm with smaller test statistics when all of them are in continuation range
#
sim_rfc_AS_DTL <- function(i=1, n=c(150, 150),
                           lambda.all=c(0.03, 0.03, 0.03),  # c(treatment1, treatment2, control)
                           p.all=c(0.3, 0.3, 0.3),               # c(treatment1, treatment2, control)
                           boundaries=c(0, 2, 2))
{
  #  on.exit(browser())
  # sample size of each arm in stage1 
  n1 <- n[1]/3 
  # 
  # endpoint component 1
  lambda1 = lambda.all[1]
  lambda2 = lambda.all[2]
  lambda0 = lambda.all[3]
  # endpoint component 2
  p1 = p.all[1]
  p2 = p.all[2]
  p0 = p.all[3]
  #
  dat1_stage1<-dat_gen_multi(nn=n1, lambda=lambda1, p=p1, arm=1) # arm1 (treatment1)
  dat2_stage1<-dat_gen_multi(nn=n1, lambda=lambda2, p=p2, arm=2) # arm2 (treatment2)
  dat0_stage1<-dat_gen_multi(nn=n1, lambda=lambda0, p=p0, arm=0) # control arm
  
  #
  W01<-Wij_cal_multi(dat0_stage1,dat1_stage1)
  W02<-Wij_cal_multi(dat0_stage1,dat2_stage1)
  #
  c0=boundaries[1] # futility boundary
  c1=boundaries[2] # efficacy boundary at stage 1
  c2=boundaries[3] # efficacy boundary at stage 2
  #
  # -------------------------Decision in stage 1-------------------------------#
  # 
  #                                     W01
  #                       futility(c0)        effective(c1)
  #                             |                 |   
  #                     0       |        1        |       4
  #                             |                 |
  #    futility(c0)  --------------------------------------------
  #                             |                 |
  # W02                 2       |     1/2/3       |       2
  #                             |                 |
  #    effective(c1) --------------------------------------------
  #                             |                 |
  #                     4       |        1        |       4
  #                             |                 |
  stop <- 0       # indicator for the above situations:
  # 0: stop trial due to futility shown from all treatments
  # 1: continue with arm 1, and drop arm 2 due to 1) futility or 2) efficacy or 3) inferior test statistics
  # 2: continue with arm 2, and drop arm 2 due to 1) futility or 2) efficacy or 3) inferior test statistics
  # 3: continue with both arm because they have the same values of test statistics
  # 4: stop trial because all treatments are dropped due to either futility or efficacy
  #
  effective <- 0
  effective_stage1 <- 0
  effective_stage2 <- 0
  winner <- 0
  # winner == 0: no treatment arm is selected;
  # winner == 1: treatment arm 1 is selected;
  # winner == 2: treatment arm 2 is selected;
  # winner == 3: both arms are selected;
  
  if (W01 > c1 || W02 > c1) {
    effective_stage1 <- 1
    if ((W01 > c1) && (W02 >= c0 && W02 <= c1)){
      stop <- 2
      winner <- 1
    }
    else if ((W02 > c1) && (W01 >= c0 && W01 <= c1)){
      stop <- 1
      winner <- 2
    }
    else {
      stop <- 4
      winner <- ifelse(all(c(W01, W02) > c1), 3, 
                       ifelse(W01 > c1, 1, 2))
    }
  } 
  else if (W01 < c0 && W02 < c0) {
    stop <- 0
  } 
  else if (W01 >= c0 && W01 <= c1 && W02 < c0) {
    stop <- 1
  } 
  else if (W02 >= c0 && W02 <= c1 && W01 < c0) {
    stop <- 2
  } 
  else if (W01 >= c0 && W01 <= c1 && W02 >= c0 && W02 <= c1) {
    stop <- ifelse(W01 > W02, 1, ifelse(W01 < W02, 2, 3))
  }
  #
  # -------------------------Decision in stage 2------------------------------ #
  #
  W_new <- -99
  if (stop %in% 1:2) {
    n2 <- n[2]/2
    arm_keep <- stop
    
    dat0_stage2 <- rbind(dat0_stage1, dat_gen_multi(n2, lambda0, p0, 0))
    dat_keep_stage2 <- rbind(get(paste0("dat", arm_keep, "_stage1")), 
                             dat_gen_multi(n2, lambda.all[arm_keep], p.all[arm_keep], arm_keep))
    
    W_new <- Wij_cal_multi(dat0_stage2, dat_keep_stage2)
    if (W_new > c2) {
      effective_stage2 <- 1
      winner <- ifelse(winner != 0, 3, arm_keep)
    }
  } 
  else if (stop == 3) {
    n2 <- n[2]/3
    dat1_stage2 <- rbind(dat1_stage1, dat_gen_multi(n2, lambda1, p1, 1))
    dat2_stage2 <- rbind(dat2_stage1, dat_gen_multi(n2, lambda2, p2, 2))
    dat0_stage2 <- rbind(dat0_stage1, dat_gen_multi(n2, lambda0, p0, 0))
    
    W01 <- Wij_cal_multi(dat0_stage2, dat1_stage2)
    W02 <- Wij_cal_multi(dat0_stage2, dat2_stage2)
    
    if (W01 > c2 || W02 > c2) {
      effective_stage2 <- 1
      if (W01 > c1 && W02 > c1) {
        winner <- 3
      }
      else if (W01 > c1) {
        winner <- 1
      }
      else if (W02 > c1) {
        winner <- 2
      }
    }
  }
  # success of the trial, i.e. global test result
  # If efficacy was shown in at least one stage, the trial is successful
  if (any(c(effective_stage1, effective_stage2) == 1)) {
    effective <- 1
  }
  res<-c(Stop = stop,
         W.01 = W01, W.02 = W02, W.new = W_new, Arm.select = winner,
         Success.stage1 = effective_stage1, Success.stage2 = effective_stage2, Success=effective)
  return(res)
}

#test4 <- sim_rfc_AS_DTL()
#
#
#
################################################################################
########### Functions for computing operating characteristics ##################
################################################################################
### 1. Type I error calculator
prob_cal <- function(n, 
                     lambda.all=c(0.06, 0.06, 0.06),  # c(treatment1, treatment2, control)
                     p.all=c(0.2, 0.2, 0.2),               # c(treatment1, treatment2, control)
                     boundaries,
                     clu,
                     n_sim,
                     sim_rfc) {
  # generate simulation result for each boundary and catch error
  results. <- tryCatch({
    pboptions(type = "none") # avoid duplicate progress bars inside
    parLapply(clu, 1:n_sim, function(x) sim_rfc(
      n = n,
      lambda.all = lambda.all,
      p.all = p.all,
      boundaries = boundaries
    ))
  }, error = function(e) {
    message("Error during parallel simulation: ", e)
    return(NULL)
  })
  
  # process results
  if (is.null(results.)) {
    return(rep(NA, 2))
  }
  
  results_df. <- as.data.frame(do.call(rbind, results.))
  
  n_total <- nrow(results_df.)
  n_win_stage1 <- nrow(results_df.[results_df.$Success.stage1 == 1,])
  n_win_stage2 <- nrow(results_df.[!results_df.$Success.stage1 == 1 & results_df.$Success.stage2 == 1,])
  
  p1 <- ifelse(n_total > 0, n_win_stage1 / n_total, NA) # alpha*_1
  p2 <- ifelse(n_total > 0, n_win_stage2 / n_total, NA) # alpha*_2 - alpha*_1
  
  # clean up in case of memory overload
  rm(list = grep("\\.$", ls(all.names = T), value = T))
  gc()
  
  c(p1, p2)
}
#
#
### 2. Power and other selection probability calculator
prob_cal_power_others <- function(n, 
                                  lambda.all=c(0.03, 0.06, 0.06),  # c(treatment1, treatment2, control)
                                  p.all=c(0.5, 0.2, 0.2),               # c(treatment1, treatment2, control)
                                  boundaries,
                                  clu,
                                  n_sim,
                                  sim_rfc) {
  # generate simulation result for each boundary and catch error
  results. <- tryCatch({
    pboptions(type = "none") # avoid duplicate progress bars inside
    parLapply(clu, 1:n_sim, function(x) sim_rfc(
      n = n,
      lambda.all = lambda.all,
      p.all = p.all,
      boundaries = boundaries
    ))
  }, error = function(e) {
    message("Error during parallel simulation: ", e)
    return(NULL)
  })
  
  # process results
  if (is.null(results.)) {
    return(rep(NA, 2))
  }
  
  results_df. <- as.data.frame(do.call(rbind, results.))
  
  n_total <- nrow(results_df.)
  n_win_total <- nrow(results_df.[results_df.$Success == 1,])
  n_select_arm1 <- nrow(results_df.[results_df.$Arm.select == 1,]) 
  n_select_arm2 <- nrow(results_df.[results_df.$Arm.select == 2,])
  n_select_both <- nrow(results_df.[results_df.$Arm.select == 3,])
  
  # calculate different probabilities
  p_power <- ifelse(n_total > 0, n_win_total / n_total, NA)
  p_select_arm1 <- ifelse(n_total > 0, n_select_arm1 / n_total, NA)
  p_select_arm2 <- ifelse(n_total > 0, n_select_arm2 / n_total, NA)
  p_select_both <- ifelse(n_total > 0, n_select_both / n_total, NA)
  
  # clean up in case of memory overload
  rm(list = grep("\\.$", ls(all.names = T), value = T))
  gc()
  
  c(Power = p_power, Select_arm1 = p_select_arm1, Select_arm2 = p_select_arm2, Select_both = p_select_both)
}
#
#
#
### 3. Boundaries finder - based on "Nelder-Mead" method
findBoundary_alphas_optim <- function(n.seq = c(150, 150), 
                                      c.start = c(2.5, 2), 
                                      lambda = c(0.06, 0.06, 0.06),
                                      p = c(0.2, 0.2, 0.2),
                                      alphas_fun = 1, 
                                      phi = NULL, 
                                      alpha = 0.05, 
                                      side = 1,
                                      clu, 
                                      n_sim,
                                      sim_rfc) {
  # Time fractions for spending function
  t_seq <- cumsum(n.seq) / sum(n.seq)
  # Alpha-spending function
  alpha_seq <- alpha_spend(iuse = alphas_fun, alpha = alpha, side = side, t = t_seq)[["pd"]]
  alpha1 <- alpha_seq[1]
  alpha2 <- alpha_seq[2]
  
  # Define constrained loss function
  loss_fn <- function(par) {
    boundaries <- c(0, par[1], par[2])
    p <- prob_cal(n = n.seq,
                  lambda.all = lambda,
                  p.all = p,
                  boundaries = boundaries,
                  clu = clu,
                  n_sim = n_sim,
                  sim_rfc = sim_rfc)
    
    if (any(is.na(p))) return(1e6)  # Penalize invalid simulation results
    
    # minimize the sum of square difference
    base_loss <- sum((p - c(alpha1, alpha2))^2)
    
    return(base_loss)
  }
  
  # Run optimization - "Nelder-Mead" method
  opt_res <- optim(par = c.start, 
                   fn = loss_fn, 
                   method = "Nelder-Mead",
                   control = list(maxit = 5000, abstol = 1e-7))
  
  # Final evaluation
  final_c1 <- opt_res$par[1]
  final_c2 <- opt_res$par[2]
  # put the optimized c1 and c2 back to function to calculate actual type I error
  final_p  <- prob_cal(n = n.seq,
                       lambda.all = lambda,
                       p.all = p,
                       boundaries = c(0, final_c1, final_c2),
                       clu = clu,
                       n_sim = n_sim,
                       sim_rfc = sim_rfc)
  
  # Return results
  data.frame(
    c1 = final_c1, 
    c2 = final_c2, 
    alpha1_target = alpha1, 
    p1 = final_p[1],
    alpha2_target = alpha2, 
    p2 = final_p[2],
    total_error = final_p[1] + final_p[2],
    converged = opt_res$convergence
  )
}
#
#
# 4. Power and selection probabily calculator wrapped with different boundaries
# Note: to maintain 5% alpha level, different sample sizes bring different boundaries values
power_wrapper_with_boundary <- function(n_total,
                                        ratio_stage1,
                                        lambda.all,
                                        p.all,
                                        alphas_fun,
                                        alpha,
                                        side,
                                        c.start,
                                        clu,
                                        n_sim,
                                        sim_rfc) {
  n1 <- floor(n_total * ratio_stage1)
  n2 <- n_total - n1
  n.seq <- c(n1, n2)
  
  # 1. Optimize boundaries for this n.seq to maintain 5% alpha level
  boundary_res <- findBoundary_alphas_optim(n.seq = n.seq,
                                            c.start = c.start,
                                            lambda = rep(0.06, 3),
                                            p = rep(0.2, 3),
                                            alphas_fun = alphas_fun,
                                            alpha = alpha,
                                            side = side,
                                            clu = clu,
                                            n_sim = n_sim,
                                            sim_rfc = sim_rfc)
  
  #if (boundary_res$converged != 0) {
  #  warning("Boundary optimization failed to converge.")
  #  return(c(Power = NA, Select_arm1 = NA, Select_both = NA, c1 = NA, c2 = NA))
  #}
  
  boundaries <- c(0, boundary_res$c1, boundary_res$c2)
  
  # 2. Evaluate power with true effect size
  power_result <- prob_cal_power_others(n = n.seq,
                                        lambda.all = lambda.all,
                                        p.all = p.all,
                                        boundaries = boundaries,
                                        clu = clu,
                                        n_sim = n_sim,
                                        sim_rfc = sim_rfc)
  
  c(Power = power_result[["Power"]],
    Select_arm1 = power_result[["Select_arm1"]],
    Select_arm2 = power_result[["Select_arm2"]],
    Select_both = power_result[["Select_both"]],
    c1 = boundary_res$c1,
    c2 = boundary_res$c2,
    Error = boundary_res$total_error)
}
#
#
#
#
# 5. Wrap boundaries, error, power and sample size finder all together  -------
Threearm_twostage_design <- function(
                                     power_target = 0.8,
                                     n_min = 250,
                                     n_max = 500,
                                     tol = 0.01,
                                     max_iter = 20,
                                     verbose = TRUE,
                                     ratio_stage1 = 1/2,
                                     lambda.all = c(0.04, 0.06, 0.06),
                                     p.all = c(0.4, 0.2, 0.2),
                                     alphas_fun = 1,
                                     alpha = 0.05,
                                     side = 1,
                                     c.start = c(2.8, 2),
                                     clu,
                                     n_sim = 10000,
                                     sim_rfc     # define stopping rules
  ) {
    iter <- 1
    best_result <- NULL
    
    while (n_max - n_min > 2 && iter <= max_iter) {
      n_mid <- floor((n_min + n_max) / 2)
      
      res <- power_wrapper_with_boundary(n_total = n_mid,
                                         ratio_stage1 = ratio_stage1,
                                         lambda.all = lambda.all,
                                         p.all = p.all,
                                         alphas_fun = alphas_fun,
                                         alpha = alpha,
                                         side = side,
                                         c.start = c.start,
                                         clu = clu,
                                         n_sim = n_sim,
                                         sim_rfc = sim_rfc)
      
      if (is.na(res[["Power"]])) {
        stop("Simulation or boundary optimization returned NA.")
      }
      # process tracking
      if (verbose) {
        cat(sprintf("n = %d → Power = %.3f, c1 = %.2f, c2 = %.2f\n",
                    n_mid, res[["Power"]], res["c1"], res["c2"]))
      }
      
      best_result <- res
      
      # Use tol to accept near match
      if (abs(res[["Power"]] - power_target) <= tol) {
        break
      }
      
      if (res[["Power"]] < power_target) {
        n_min <- n_mid + 1
      } else {
        n_max <- n_mid
      }
      
      iter <- iter + 1
    }
    
    return(list(n_total = n_mid,
                power = best_result[["Power"]],
                error = best_result[["Error"]],
                select_arm1 = best_result[["Select_arm1"]],
                select_arm2 = best_result[["Select_arm2"]],
                select_both = best_result[["Select_both"]],
                c1 = best_result[["c1"]],
                c2 = best_result[["c2"]]))
  }
#
#
#