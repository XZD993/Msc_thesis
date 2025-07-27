################################################################################
#FWER, Power, andOther operating characteristics in standard two-arm group sequential trial ##
############# For calculation in Section 6.2, 7.3, and 8 of thesis #############
#
#### package loading
library(mvtnorm)
#
#### global variable
## W_arm1_stage1, W_arm2_stage1, W_arm1_stage2, W_arm2_stage2
# covariance matrix
Sigma <- matrix(c(
  1, 0.5, 1/sqrt(2), 1/(2*sqrt(2)),
  0.5, 1, 1/(2*sqrt(2)), 1/sqrt(2),
  1/sqrt(2), 1/(2*sqrt(2)), 1, 0.5,
  1/(2*sqrt(2)), 1/sqrt(2), 0.5, 1
), nrow = 4, byrow = TRUE)  
#
#
################## Section 6.2 ----- FWER calculation ##########################
# mean vector under null hypothesis
vec_mean_H0 <- c(rep(0, 4))  

## Type I error at stage 1
FWER_1 <- 1 - pmvnorm(
  upper = c(2.87, 2.87),
  mean = vec_mean_H0[1:2],
  sigma = Sigma[1:2, 1:2]
)
#
## Type I error at stage 2
FWER_2 <- pmvnorm(
  lower = c(0, 0),
  upper = c(2.87, 2.87),
  mean = vec_mean_H0[1:2],
  sigma = Sigma[1:2, 1:2]
) - pmvnorm(
  lower = c(0, 0, -Inf, -Inf),
  upper = c(2.87, 2.87, 1.95, 1.95),
  mean = vec_mean_H0,
  sigma = Sigma
)
#
FWER_3 <- pmvnorm(
  lower = c(-Inf, 0),
  upper = c(0, 2.87),
  mean = vec_mean_H0[1:2],
  sigma = Sigma[1:2, 1:2]
) - pmvnorm(
  lower = c(-Inf, 0, -Inf, -Inf),
  upper = c(0, 2.87, 1.95, 1.95),
  mean = vec_mean_H0,
  sigma = Sigma
)
#
FWER_4 <- pmvnorm(
  lower = c(0, -Inf),
  upper = c(2.87, 0),
  mean = vec_mean_H0[1:2],
  sigma = Sigma[1:2, 1:2]
) - pmvnorm(
  lower = c(0, -Inf, -Inf, -Inf),
  upper = c(2.87, 0, 1.95, 1.95),
  mean = vec_mean_H0,
  sigma = Sigma
)
#
FWER <- FWER_1 + FWER_2 + FWER_3 + FWER_4
#
#
############ Section 7.3 -------- Power calculation ############################
#
### Disjunctive power
# mean vector under global alternative hypothesis
vec_mean_dis <- c(2, 2, 2*sqrt(2), 2*sqrt(2)) 
## Power at stage 1
dis_power_1 <- 1 - pmvnorm(
  upper = c(2.87, 2.87),
  mean = vec_mean_dis[1:2],
  sigma = Sigma[1:2, 1:2]
)
#
## Power at stage 2
dis_power_2 <- pmvnorm(
  lower = c(0, 0),
  upper = c(2.87, 2.87),
  mean = vec_mean_dis[1:2],
  sigma = Sigma[1:2, 1:2]
) - pmvnorm(
  lower = c(0, 0, -Inf, -Inf),
  upper = c(2.87, 2.87, 1.95, 1.95),
  mean = vec_mean_dis,
  sigma = Sigma
)
#
dis_power_3 <- pmvnorm(
  lower = c(-Inf, 0),
  upper = c(0, 2.87),
  mean = vec_mean_dis[1:2],
  sigma = Sigma[1:2, 1:2]
) - pmvnorm(
  lower = c(-Inf, 0, -Inf, -Inf),
  upper = c(0, 2.87, 1.95, 1.95),
  mean = vec_mean_dis,
  sigma = Sigma
  )
#
dis_power_4 <- pmvnorm(
  lower = c(0, -Inf),
  upper = c(2.87, 0),
  mean = vec_mean_dis[1:2],
  sigma = Sigma[1:2, 1:2]
) - pmvnorm(
  lower = c(0, -Inf, -Inf, -Inf),
  upper = c(2.87, 0, 1.95, 1.95),
  mean = vec_mean_dis,
  sigma = Sigma
  )
#
## Disjunctive power
dis_power <- dis_power_1 + dis_power_2 + dis_power_3 + dis_power_4
#
#
#
#
### LFC power
# mean vector under global alternative hypothesis
vec_mean_LFC <- c(2, 0, 2*sqrt(2), 0) 
## Power at stage 1
LFC_power_1 <- 1 - pmvnorm(
  upper = c(2.87, 2.87),
  mean = vec_mean_LFC[1:2],
  sigma = Sigma[1:2, 1:2]
)
#
## Power at stage 2
LFC_power_2 <- pmvnorm(
  lower = c(0, 0),
  upper = c(2.87, 2.87),
  mean = vec_mean_LFC[1:2],
  sigma = Sigma[1:2, 1:2]
) - pmvnorm(
  lower = c(0, 0, -Inf, -Inf),
  upper = c(2.87, 2.87, 1.95, 1.95),
  mean = vec_mean_LFC,
  sigma = Sigma
)
#
LFC_power_3 <- pmvnorm(
  lower = c(-Inf, 0),
  upper = c(0, 2.87),
  mean = vec_mean_LFC[1:2],
  sigma = Sigma[1:2, 1:2]
) - pmvnorm(
  lower = c(-Inf, 0, -Inf, -Inf),
  upper = c(0, 2.87, 1.95, 1.95),
  mean = vec_mean_LFC,
  sigma = Sigma
)
#
LFC_power_4 <- pmvnorm(
  lower = c(0, -Inf),
  upper = c(2.87, 0),
  mean = vec_mean_LFC[1:2],
  sigma = Sigma[1:2, 1:2]
) - pmvnorm(
  lower = c(0, -Inf, -Inf, -Inf),
  upper = c(2.87, 0, 1.95, 1.95),
  mean = vec_mean_LFC,
  sigma = Sigma
)
#
## Disjunctive power
LFC_power <- LFC_power_1 + LFC_power_2+ LFC_power_3 + LFC_power_4
#
#
#
############## Section 8 --------- Other operating characteristics #############
#
#
#### Section 8.1
## The probability of both arms being recommended if both are genuinely effective - p1
# mean vector under global alternative hypothesis
vec_mean_dis <- c(2, 2, 2*sqrt(2), 2*sqrt(2)) 
#
# Term 1: P(W_{01}^1 > 2.87, W_{02}^1 > 2.87)
p1.term1 <- pmvnorm(
  lower = c(2.87, 2.87),
  upper = rep(Inf, 2),
  mean = vec_mean_dis[1:2],
  sigma = Sigma[1:2, 1:2]
)
# Term 2: P(0 ≤ W_{01}^1 ≤ 2.87, 0 ≤ W_{02}^1 ≤ 2.87, W_{01}^2 > 1.95, W_{02}^2 > 1.95)
p1.term2 <- pmvnorm(
  lower = c(0, 0, 1.95, 1.95),
  upper = c(2.87, 2.87, Inf, Inf),
  mean = vec_mean_dis,
  sigma = Sigma
)
# p1
p1 <- p1.term1 + p1.term2
#
#
#### Section 8.2
### The probability of one specific experimental arm being recommended - p2
#
# mean vector under global alternative hypothesis
vec_mean_LFC <- c(2, 0, 2*sqrt(2), 0) 
#
# Term 1: P(W_{01}^1 > 2.87, W_{02}^1 > 2.87)
p2.term1 <- pmvnorm(
  lower = c(2.87, -Inf),
  upper = c(Inf, 2.87),
  mean = vec_mean_LFC[1:2],
  sigma = Sigma[1:2, 1:2]
)
# Term 2: P(0 ≤ W_{01}^1 ≤ 2.87, 0 ≤ W_{02}^1 ≤ 2.87, W_{01}^2 > 1.95, W_{02}^2 > 1.95)
p2.term2 <- pmvnorm(
  lower = c(0, -Inf, 1.95, -Inf),
  upper = c(2.87, 2.87, Inf, 1.95),
  mean = vec_mean_LFC,
  sigma = Sigma
)
# p1
p2 <- p2.term1 + p2.term2
#
#
#### Section 8.3
#### The conditional probability - p3
#
### Input
## Observed value
w_a <- c(2, 2)
#
## mean vectors under global alternative hypothesis
mean_a <- c(2, 2)
mean_b <- c(2*sqrt(2), 2*sqrt(2))
#
## covariance matrices:
Sigma_aa <- matrix(c(1, 0.5,
                     0.5, 1),
                   byrow = T, nrow = 2)
Sigma_ab <- matrix(c(1/sqrt(2), 1/(2*sqrt(2)),
                     1/(2*sqrt(2)), 1/sqrt(2)),
                   byrow = T, nrow = 2)
Sigma_ba <- Sigma_ab
Sigma_bb <- Sigma_aa
#
### Calculate conditional distribution
#
## Conditional mean vector
mean_cond <- mean_b + Sigma_ba %*% solve(Sigma_aa) %*% (w_a - mean_a)
#
## Conditional covariance matrix
Sigma_cond <- Sigma_bb - Sigma_ba %*% solve(Sigma_aa) %*% Sigma_ab
#
### Calculate the probability
# p3
p3 <- 1 - pmvnorm(
  lower = c(-Inf, -Inf),
  upper = c(1.95, 1.95),
  mean = as.vector(mean_cond[,1]),
  sigma = Sigma_cond
)
#
#