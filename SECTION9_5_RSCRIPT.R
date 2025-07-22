# Msc thesis:
# A group sequential three-arm two-stage trial design with rank-based endpoints
###################### Code for section 9.5 ####################################
################################################################################
#
#
# ----------------------- Global settings --------------------------------------
source("FUNCTIONS.R")

Sys.setenv(OPENBLAS_NUM_THREADS = 4)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")
#
n_cores <- 48
#
#
#
# ------- Trial characteristics calculation for different stopping rules -------
#
### General setting:
# 1) alpha: 5% (one-sided)
# 2) target power: 80%
# 3) LFC hypothesis: 
# lambda.all = c(0.04, 0.06, 0.06), 
# which means 16.5% patient live longer than 30 days in control group and 30% in arm 1, and
# p.all = c(0.4, 0.2, 0.2),
# which means 20% kidney recovery in control group and 40% in arm 1.
#
#
### --------------- MJW stopping rule
#
start_time_MJW <- Sys.time()
clu <- makeCluster(n_cores)
on.exit(stopCluster(clu)) # auto shutdown cluster if error
clusterEvalQ(clu, source("FUNCTIONS.R"))
#
MJW <- Threearm_twostage_design(clu = clu, sim_rfc = sim_rfc_MJW, n_sim = 100000)
#
#stopCluster(clu)
end_time_MJW <- Sys.time()
cat("\nRun completed in", round(difftime(end_time_MJW, start_time_MJW, units = "mins"), 2), "minutes.\n")
print(MJW)
#
#
## --------------- MJW-DTL stopping rule
#
start_time_MJW_DTL <- Sys.time()
clu <- makeCluster(n_cores)
on.exit(stopCluster(clu)) # auto shutdown cluster if error
clusterEvalQ(clu, source("Rcpp_MULTI_FUNCTIONS_ALPHAandPOWER.R"))
#
MJW_DTL <- Threearm_twostage_design(clu = clu, sim_rfc = sim_rfc_MJW_DTL, n_sim = 100000)
#
#stopCluster(clu)
end_time_MJW_DTL <- Sys.time()
cat("\nRun completed in", round(difftime(end_time_MJW_DTL, start_time_MJW_DTL, units = "mins"), 2), "minutes.\n")
print(MJW_DTL)
#
#
## --------------- AS stopping rule
#
start_time_AS <- Sys.time()
clu <- makeCluster(n_cores)
on.exit(stopCluster(clu)) # auto shutdown cluster if error
clusterEvalQ(clu, source("Rcpp_MULTI_FUNCTIONS_ALPHAandPOWER.R"))
#
AS <- Threearm_twostage_design(clu = clu, sim_rfc = sim_rfc_AS, n_sim = 100000)
#
stopCluster(clu)
end_time_AS <- Sys.time()
cat("\nRun completed in", round(difftime(end_time_AS, start_time_AS, units = "mins"), 2), "minutes.\n")
print(AS)
#
#
## --------------- AS-DTL stopping rule
#
start_time_AS_DTL <- Sys.time()
clu <- makeCluster(n_cores)
on.exit(stopCluster(clu)) # auto shutdown cluster if error
clusterEvalQ(clu, source("Rcpp_MULTI_FUNCTIONS_ALPHAandPOWER.R"))
#
AS_DTL <- Threearm_twostage_design(clu = clu, sim_rfc = sim_rfc_AS_DTL, n_sim = 100000)
#
stopCluster(clu)
end_time_AS_DTL <- Sys.time()
cat("\nRun completed in", round(difftime(end_time_AS_DTL, start_time_AS_DTL, units = "mins"), 2), "minutes.\n")
print(AS_DTL)
#
#
# ----------------------- Process the output -----------------------------------
res_table_MJW <- do.call(cbind, MJW)
res_table_MJW_DTL <- do.call(cbind, MJW_DTL)
res_table_AS <- do.call(cbind, AS)
res_table_AS_DTL <- do.call(cbind, AS_DTL)
res_table <- rbind(res_table_MJW, res_table_MJW_DTL, res_table_AS, res_table_AS_DTL)
res_table <- as.data.frame(res_table)
res_table$run_time <- c(round(difftime(end_time_MJW, start_time_MJW, units = "mins"), 2),
                        round(difftime(end_time_MJW_DTL, start_time_MJW_DTL, units = "mins"), 2),
                        round(difftime(end_time_AS, start_time_AS, units = "mins"), 2),
                        round(difftime(end_time_AS_DTL, start_time_AS_DTL, units = "mins"), 2))
rownames(res_table) <- c("MJW", "MJW_DTL", "AS", "AS_DTL")
saveRDS(res_table, file = "TABLE4_DATA.rds")

#
#