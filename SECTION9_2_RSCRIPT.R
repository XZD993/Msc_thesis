# Msc thesis:
# A group sequential three-arm two-stage trial design with rank-based endpoints
########################## Code for Section 9.2 ################################
################################################################################
#
#
#
# ----------------------- Global settings --------------------------------------
source("FUNCTIONS.R")

Sys.setenv(OPENBLAS_NUM_THREADS = 4)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")

n_cores <- 48
#
#
#
# ---------------------- Boundaries finding for 4 designs ----------------------
#
#
# ------- MJW
# Set up the cluster
start_time1 <- Sys.time()

clu <- makeCluster(n_cores)
on.exit(stopCluster(clu)) # auto shutdown cluster if error

clusterEvalQ(clu, source("FUNCTIONS.R"))


result_MJW <- findBoundary_alphas_optim(n.seq = c(150, 150), c.start = c(2.8, 2), alphas_fun = 1,
                                        clu = clu, n_sim = 100000, sim_rfc = sim_rfc_MJW)

stopCluster(clu)

end_time1 <- Sys.time()
cat("\nRun completed in", round(difftime(end_time1, start_time1, units = "mins"), 2), "minutes.\n")

print(result_MJW)
#
#
#
# ------- MJW-DTL
# Set up the cluster
start_time2 <- Sys.time()

clu <- makeCluster(n_cores)
on.exit(stopCluster(clu)) # auto shutdown cluster if error

clusterEvalQ(clu, source("FUNCTIONS.R"))

result_MJW_DTL <- findBoundary_alphas_optim(n.seq = c(150, 150), c.start = c(2.8, 2), alphas_fun = 1,
                                            clu = clu, n_sim = 100000, sim_rfc = sim_rfc_MJW_DTL)

stopCluster(clu)

end_time2 <- Sys.time()
cat("\nRun completed in", round(difftime(end_time2, start_time2, units = "mins"), 2), "minutes.\n")

print(result_MJW_DTL)
#
#
#
# ------- AS
# Set up the cluster
start_time3 <- Sys.time()

clu <- makeCluster(n_cores)
on.exit(stopCluster(clu)) # auto shutdown cluster if error

clusterEvalQ(clu, source("FUNCTIONS.R"))

result_AS <- findBoundary_alphas_optim(n.seq = c(150, 150), c.start = c(2.8, 2), alphas_fun = 1,
                                       clu = clu, n_sim = 100000, sim_rfc = sim_rfc_AS)

stopCluster(clu)

end_time3 <- Sys.time()
cat("\nRun completed in", round(difftime(end_time3, start_time3, units = "mins"), 2), "minutes.\n")

print(result_AS)
#
#
#
# ------- AS-DTL
# Set up the cluster
start_time4 <- Sys.time()

clu <- makeCluster(n_cores)
on.exit(stopCluster(clu)) # auto shutdown cluster if error

clusterEvalQ(clu, source("FUNCTIONS.R"))


result_AS_DTL <- findBoundary_alphas_optim(n.seq = c(150, 150), c.start = c(2.8, 2), alphas_fun = 1,
                                           clu = clu, n_sim = 100000, sim_rfc = sim_rfc_AS_DTL)

stopCluster(clu)

end_time4 <- Sys.time()
cat("\nRun completed in", round(difftime(end_time4, start_time4, units = "mins"), 2), "minutes.\n")

print(result_AS_DTL)
#
#
# Aggregate results
res_bound_find_4rules <- rbind(result_MJW, result_MJW_DTL, result_AS, result_AS_DTL)
res_bound_find_4rules$design <- c("MJW", "MJW-DTL", "AS", "AS-DTL")
res_bound_find_4rules$run.time <- c(round(difftime(end_time1, start_time1, units = "mins"), 2),
                                    round(difftime(end_time2, start_time2, units = "mins"), 2),
                                    round(difftime(end_time3, start_time3, units = "mins"), 2),
                                    round(difftime(end_time4, start_time4, units = "mins"), 2))
#
#
saveRDS(res_bound_find_4rules, file = "TABLE1_DATA.rds")
#
#
#
#
# ------------------------------------------------------------------------------
### Appendix A.2 
#The impact of starting point on the performance of Nelder-Mead method
## use MJW-DTL as an example
# n_sim == 100000, maxit==5000, abstol == 1e-7，sweeping
ini_simplex <- list(c(1, 1),
                    c(1.5, 1),
                    c(1.5, 1.5),
                    c(2, 1.5),
                    c(2, 2),
                    c(2.5, 2),
                    c(2.5, 2.5),
                    c(3, 2.5),
                    c(3, 3),
                    c(3.5, 3),
                    c(3.5, 3.5),
                    c(4, 3.5),
                    c(4, 4),
                    c(4.5, 4),
                    c(4.5, 4.5),
                    c(5, 4.5),
                    c(5, 5))
res_ini_simplex <- data.frame()

clu <- makeCluster(n_cores)
on.exit(stopCluster(clu)) # auto shutdown cluster if error
clusterEvalQ(clu, source("FUNCTIONS.R"))

for (i in 1:length(ini_simplex)) {
  start_time <- Sys.time()
  
  res_i <- findBoundary_alphas_optim(n.seq = c(150, 150), c.start = ini_simplex[[i]], alphas_fun = 1,
                                     clu = clu, n_sim = 100000, sim_rfc = sim_rfc_MJW_DTL)
  
  end_time <- Sys.time()
  
  res_i$run_time <- round(difftime(end_time, start_time, units = "mins"), 2)
  
  res_ini_simplex <- rbind(res_ini_simplex, res_i)
}


ini_simplex_dat <- as.data.frame(do.call(rbind, ini_simplex))
res_ini_simplex_sweeping <- cbind(ini_simplex_dat, res_ini_simplex)
names(res_ini_simplex_sweeping)[1] <- "c1_initial"
names(res_ini_simplex_sweeping)[2] <- "c2_initial"
saveRDS(res_ini_simplex_sweeping, file = "TABLE_MULTI_MULTI_ini_simplex_sweep.rds")
save.image(file = "MULTI_ini_simplex_sweep_20072025.RData")
