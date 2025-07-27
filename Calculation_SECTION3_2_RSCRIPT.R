################################################################################
### The FWER and power calculation in standard two-arm group sequential trail ##
############# For calculation in Section 3.2 of thesis #########################
#
#### package loading
library(mvtnorm)
#
#
#### -------------------  FWER calculation
## Type I error at stage 1
FWER_1 <- pnorm(-1.96)
#
## Type I error at stage 2
vec_mean_H0 <- c(0, 0)  # mean vector under null hypothesis
cov <- 1/sqrt(2) # covariance between Z1 and Z2
mat_cov <- matrix(c(1, cov,
                    cov, 1),
                  nrow = 2, ncol = 2)
FWER_2.1 <- 
  pmvnorm(lower = c(-Inf, -Inf), upper = c(1.96, Inf), mean = vec_mean_H0, sigma = mat_cov)-
  pmvnorm(lower = c(-Inf, -Inf), upper = c(0, Inf), mean = vec_mean_H0, sigma = mat_cov)-
  pmvnorm(lower = c(-Inf, -Inf), upper = c(1.96, 1.96), mean = vec_mean_H0, sigma = mat_cov)+
  pmvnorm(lower = c(-Inf, -Inf), upper = c(0, 1.96), mean = vec_mean_H0, sigma = mat_cov)
# It can also be calculated as below:
FWER_2.2 <- 
  pmvnorm(lower = c(0, -Inf), upper = c(1.96, Inf), mean = vec_mean_H0, sigma = mat_cov)-
  pmvnorm(lower = c(0, -Inf), upper = c(1.96, 1.96), mean = vec_mean_H0, sigma = mat_cov)
#
## FWER
FWER <- FWER_1 + FWER_2.2
#
#
#### ------------------- Power calculation
## Power at stage 1
Power_1 <- 1- pnorm(1.96, 2, 1)
#
## Power at stage 2
vec_mean_H1 <- c(2, 2*sqrt(2))  # mean vector under alternative hypothesis
cov <- 1/sqrt(2) # covariance between Z1 and Z2
mat_cov <- matrix(c(1, cov,
                    cov, 1),
                  nrow = 2, ncol = 2)
Power_2.1 <- 
  pmvnorm(lower = c(-Inf, -Inf), upper = c(1.96, Inf), mean = vec_mean_H1, sigma = mat_cov)-
  pmvnorm(lower = c(-Inf, -Inf), upper = c(0, Inf), mean = vec_mean_H1, sigma = mat_cov)-
  pmvnorm(lower = c(-Inf, -Inf), upper = c(1.96, 1.96), mean = vec_mean_H1, sigma = mat_cov)+
  pmvnorm(lower = c(-Inf, -Inf), upper = c(0, 1.96), mean = vec_mean_H1, sigma = mat_cov)
# It can also be calculated as below:
Power_2.2 <- 
  pmvnorm(lower = c(0, -Inf), upper = c(1.96, Inf), mean = vec_mean_H1, sigma = mat_cov)-
  pmvnorm(lower = c(0, -Inf), upper = c(1.96, 1.96), mean = vec_mean_H1, sigma = mat_cov)
#
## FWER
Power <- Power_1 + Power_2.2

