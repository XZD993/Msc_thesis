################################################################################
#FWER, Power, and Other operating characteristics in standard two-arm group sequential trial ##
############# For calculation in Section 6.2, 7.3, and 8 of thesis #############
#
#### package loading
library(mvtnorm)
source("FUNCTIONS.R")
#
#
lambda.all <- c(0.04, 0.06, 0.06)  # c(treatment1, treatment2, control)
p.all <- c(0.4, 0.2, 0.2)               # c(treatment1, treatment2, control)
#
##### MJW
#
N_MJW <- 344
#
### Endpoints
#
# endpoint component A
lambda1 = lambda.all[1]
lambda2 = lambda.all[2]
lambda0 = lambda.all[3]
# endpoint component B
p1 = p.all[1]
p2 = p.all[2]
p0 = p.all[3]
#
### Data generating
# stage 1
dat1_stage1<-dat_gen_multi(nn=N_MJW/6, lambda=lambda1, p=p1, arm=1) # arm1 (treatment1)
dat2_stage1<-dat_gen_multi(nn=N_MJW/6, lambda=lambda2, p=p2, arm=2) # arm2 (treatment2)
dat0_stage1<-dat_gen_multi(nn=N_MJW/6, lambda=lambda0, p=p0, arm=0) # control arm
# stage 2
dat1_stage2 <- rbind(dat1_stage1, 
                     dat_gen_multi(N_MJW/6, lambda=lambda1, p=p1, arm=1))
dat2_stage2 <- rbind(dat2_stage1, 
                     dat_gen_multi(N_MJW/6, lambda=lambda2, p=p2, arm=2))
dat0_stage2 <- rbind(dat0_stage1, 
                     dat_gen_multi(N_MJW/6, lambda=lambda0, p=p0, arm=0))
#
### Test statistics W01_1, W02_1, W01_2, W02_2
W01_1<-Wij_cal_multi(dat0_stage1,dat1_stage1)
W02_1<-Wij_cal_multi(dat0_stage1,dat2_stage1)
W01_2<-Wij_cal_multi(dat0_stage2,dat1_stage2)
W02_2<-Wij_cal_multi(dat0_stage2,dat2_stage2)

