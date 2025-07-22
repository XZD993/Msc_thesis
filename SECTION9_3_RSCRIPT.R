# Msc thesis:
# A group sequential three-arm two-stage trial design with rank-based endpoints
###################### Code for section 9.3 ####################################
################################################################################
#
#
#
#
# ----------------------- Global settings --------------------------------------
source("FUNCTIONS.R")
#
# Boundaries values for all four stopping rules
Bounds <- readRDS("SECTION9_2_RSCRIPT.rds")
# Note: If no available boundaries, please use "MULTI_Section9_6.R" first to get corresponding boundaries

Sys.setenv(OPENBLAS_NUM_THREADS = 4)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")
#
n_cores <- 48
#
#
#
# ------------ Empirical power for 4 designs under different hypotheses --------
#### ------- 1. MJW
#
### boundaries values
c1 <- round(Bounds$c1[Bounds$design == "MJW"],2)
c2 <- round(Bounds$c2[Bounds$design == "MJW"],2)
#
### Simulation
#
start_time <- Sys.time()
clu <- makeCluster(n_cores)
on.exit(stopCluster(clu)) # auto shutdown cluster if error
clusterEvalQ(clu, source("FUNCTIONS.R"))
clusterExport(clu, varlist = c("c1", "c2"))

## 1. lambda.all=c(0.03, 0.06, 0.06), p.all=c(0.5, 0.2, 0.2) - LFC hypothesis
p_MJW_1 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.03, 0.06, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.5, 0.2, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW
)
#
## 2. lambda.all=c(0.03, 0.05, 0.06), p.all=c(0.5, 0.3, 0.2)
p_MJW_2 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.03, 0.05, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.5, 0.3, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW
)
#
## 3. lambda.all=c(0.03, 0.04, 0.06), p.all=c(0.5, 0.4, 0.2)
p_MJW_3 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.03, 0.04, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.5, 0.4, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW
)
#
## 4. lambda.all=c(0.03, 0.03, 0.06), p.all=c(0.5, 0.5, 0.2) - global alternative hypothesis
p_MJW_4 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.03, 0.03, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.5, 0.5, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW
)
#
## 5. lambda.all=c(0.04, 0.06, 0.06), p.all=c(0.4, 0.2, 0.2)
p_MJW_5 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.04, 0.06, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.4, 0.2, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW
)
#
## 6. lambda.all=c(0.04, 0.05, 0.06), p.all=c(0.4, 0.3, 0.2)
p_MJW_6 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.04, 0.05, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.4, 0.3, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW
)
#
## 7. lambda.all=c(0.04, 0.04, 0.06), p.all=c(0.4, 0.4, 0.2)
p_MJW_7 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.04, 0.04, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.4, 0.4, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW
)
#
## 8. lambda.all=c(0.05, 0.06, 0.06), p.all=c(0.3, 0.2, 0.2)
p_MJW_8 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.05, 0.06, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.3, 0.2, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW
)
#
## 9. lambda.all=c(0.05, 0.05, 0.06), p.all=c(0.3, 0.3, 0.2)
p_MJW_9 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.05, 0.05, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.3, 0.3, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW
)
stopCluster(clu)

# aggregate data
table_MJW <- do.call(rbind, mget(paste0("p_MJW_", 1:9)))

end_time <- Sys.time()
cat("\nRun completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes.\n")
print(table_MJW)
#
#
#### ------- 2. MJW-DTL
### boundaries values
c1 <- round(Bounds$c1[Bounds$design == "MJW-DTL"],2)
c2 <- round(Bounds$c2[Bounds$design == "MJW-DTL"],2)
#
### Simulation
#
start_time <- Sys.time()
clusterExport(clu, varlist = c("c1", "c2"))
#
## 1. lambda.all=c(0.03, 0.06, 0.06), p.all=c(0.5, 0.2, 0.2) - LFC hypothesis
p_MJW_DTL_1 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.03, 0.06, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.5, 0.2, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW_DTL
)
#
## 2. lambda.all=c(0.03, 0.05, 0.06), p.all=c(0.5, 0.3, 0.2)
p_MJW_DTL_2 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.03, 0.05, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.5, 0.3, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW_DTL
)
#
## 3. lambda.all=c(0.03, 0.04, 0.06), p.all=c(0.5, 0.4, 0.2)
p_MJW_DTL_3 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.03, 0.04, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.5, 0.4, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW_DTL
)
#
## 4. lambda.all=c(0.03, 0.03, 0.06), p.all=c(0.5, 0.5, 0.2) - global alternative hypothesis
p_MJW_DTL_4 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.03, 0.03, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.5, 0.5, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW_DTL
)
#
## 5. lambda.all=c(0.04, 0.06, 0.06), p.all=c(0.4, 0.2, 0.2)
p_MJW_DTL_5 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.04, 0.06, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.4, 0.2, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW_DTL
)
#
## 6. lambda.all=c(0.04, 0.05, 0.06), p.all=c(0.4, 0.3, 0.2)
p_MJW_DTL_6 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.04, 0.05, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.4, 0.3, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW_DTL
)
#
## 7. lambda.all=c(0.04, 0.04, 0.06), p.all=c(0.4, 0.4, 0.2)
p_MJW_DTL_7 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.04, 0.04, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.4, 0.4, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW_DTL
)
#
## 8. lambda.all=c(0.05, 0.06, 0.06), p.all=c(0.3, 0.2, 0.2)
p_MJW_DTL_8 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.05, 0.06, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.3, 0.2, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW_DTL
)
#
## 9. lambda.all=c(0.05, 0.05, 0.06), p.all=c(0.3, 0.3, 0.2)
p_MJW_DTL_9 <- prob_cal_power(n = c(150, 150),
                          lambda.all=c(0.05, 0.05, 0.06),  # c(treatment1, treatment2, control)
                          p.all=c(0.3, 0.3, 0.2),               # c(treatment1, treatment2, control)
                          boundaries = c(0, c1, c2),
                          clu = clu,
                          n_sim = 100000,
                          sim_rfc = sim_rfc_MJW_DTL
)
stopCluster(clu)

# aggregate data
table_MJW_DTL <- do.call(rbind, mget(paste0("p_MJW_DTL_", 1:9)))

end_time <- Sys.time()
cat("\nRun completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes.\n")
print(table_MJW_DTL)
#
#
#### ------- 3. AS
### boundaries values
c1 <- round(Bounds$c1[Bounds$design == "AS"],2)
c2 <- round(Bounds$c2[Bounds$design == "AS"],2)
#
### Simulation
#
start_time <- Sys.time()
clusterExport(clu, varlist = c("c1", "c2"))
#
## 1. lambda.all=c(0.03, 0.06, 0.06), p.all=c(0.5, 0.2, 0.2) - LFC hypothesis
p_AS_1 <- prob_cal_power(n = c(150, 150),
                              lambda.all=c(0.03, 0.06, 0.06),  # c(treatment1, treatment2, control)
                              p.all=c(0.5, 0.2, 0.2),               # c(treatment1, treatment2, control)
                              boundaries = c(0, c1, c2),
                              clu = clu,
                              n_sim = 100000,
                              sim_rfc = sim_rfc_AS
)
#
## 2. lambda.all=c(0.03, 0.05, 0.06), p.all=c(0.5, 0.3, 0.2)
p_AS_2 <- prob_cal_power(n = c(150, 150),
                              lambda.all=c(0.03, 0.05, 0.06),  # c(treatment1, treatment2, control)
                              p.all=c(0.5, 0.3, 0.2),               # c(treatment1, treatment2, control)
                              boundaries = c(0, c1, c2),
                              clu = clu,
                              n_sim = 100000,
                              sim_rfc = sim_rfc_AS
)
#
## 3. lambda.all=c(0.03, 0.04, 0.06), p.all=c(0.5, 0.4, 0.2)
p_AS_3 <- prob_cal_power(n = c(150, 150),
                              lambda.all=c(0.03, 0.04, 0.06),  # c(treatment1, treatment2, control)
                              p.all=c(0.5, 0.4, 0.2),               # c(treatment1, treatment2, control)
                              boundaries = c(0, c1, c2),
                              clu = clu,
                              n_sim = 100000,
                              sim_rfc = sim_rfc_AS
)
#
## 4. lambda.all=c(0.03, 0.03, 0.06), p.all=c(0.5, 0.5, 0.2) - global alternative hypothesis
p_AS_4 <- prob_cal_power(n = c(150, 150),
                              lambda.all=c(0.03, 0.03, 0.06),  # c(treatment1, treatment2, control)
                              p.all=c(0.5, 0.5, 0.2),               # c(treatment1, treatment2, control)
                              boundaries = c(0, c1, c2),
                              clu = clu,
                              n_sim = 100000,
                              sim_rfc = sim_rfc_AS
)
#
## 5. lambda.all=c(0.04, 0.06, 0.06), p.all=c(0.4, 0.2, 0.2)
p_AS_5 <- prob_cal_power(n = c(150, 150),
                              lambda.all=c(0.04, 0.06, 0.06),  # c(treatment1, treatment2, control)
                              p.all=c(0.4, 0.2, 0.2),               # c(treatment1, treatment2, control)
                              boundaries = c(0, c1, c2),
                              clu = clu,
                              n_sim = 100000,
                              sim_rfc = sim_rfc_AS
)
#
## 6. lambda.all=c(0.04, 0.05, 0.06), p.all=c(0.4, 0.3, 0.2)
p_AS_6 <- prob_cal_power(n = c(150, 150),
                              lambda.all=c(0.04, 0.05, 0.06),  # c(treatment1, treatment2, control)
                              p.all=c(0.4, 0.3, 0.2),               # c(treatment1, treatment2, control)
                              boundaries = c(0, c1, c2),
                              clu = clu,
                              n_sim = 100000,
                              sim_rfc = sim_rfc_AS
)
#
## 7. lambda.all=c(0.04, 0.04, 0.06), p.all=c(0.4, 0.4, 0.2)
p_AS_7 <- prob_cal_power(n = c(150, 150),
                              lambda.all=c(0.04, 0.04, 0.06),  # c(treatment1, treatment2, control)
                              p.all=c(0.4, 0.4, 0.2),               # c(treatment1, treatment2, control)
                              boundaries = c(0, c1, c2),
                              clu = clu,
                              n_sim = 100000,
                              sim_rfc = sim_rfc_AS
)
#
## 8. lambda.all=c(0.05, 0.06, 0.06), p.all=c(0.3, 0.2, 0.2)
p_AS_8 <- prob_cal_power(n = c(150, 150),
                              lambda.all=c(0.05, 0.06, 0.06),  # c(treatment1, treatment2, control)
                              p.all=c(0.3, 0.2, 0.2),               # c(treatment1, treatment2, control)
                              boundaries = c(0, c1, c2),
                              clu = clu,
                              n_sim = 100000,
                              sim_rfc = sim_rfc_AS
)
#
## 9. lambda.all=c(0.05, 0.05, 0.06), p.all=c(0.3, 0.3, 0.2)
p_AS_9 <- prob_cal_power(n = c(150, 150),
                              lambda.all=c(0.05, 0.05, 0.06),  # c(treatment1, treatment2, control)
                              p.all=c(0.3, 0.3, 0.2),               # c(treatment1, treatment2, control)
                              boundaries = c(0, c1, c2),
                              clu = clu,
                              n_sim = 100000,
                              sim_rfc = sim_rfc_AS
)


# aggregate data
table_AS <- do.call(rbind, mget(paste0("p_AS_", 1:9)))

end_time <- Sys.time()
cat("\nRun completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes.\n")
print(table_AS)
#
#
#### ------- 4. AS-DTL
### boundaries values
c1 <- round(Bounds$c1[Bounds$design == "AS-DTL"],2)
c2 <- round(Bounds$c2[Bounds$design == "AS-DTL"],2)
#
### Simulation
#
start_time <- Sys.time()
clusterExport(clu, varlist = c("c1", "c2"))
#
## 1. lambda.all=c(0.03, 0.06, 0.06), p.all=c(0.5, 0.2, 0.2) - LFC hypothesis
p_AS_DTL_1 <- prob_cal_power(n = c(150, 150),
                         lambda.all=c(0.03, 0.06, 0.06),  # c(treatment1, treatment2, control)
                         p.all=c(0.5, 0.2, 0.2),               # c(treatment1, treatment2, control)
                         boundaries = c(0, c1, c2),
                         clu = clu,
                         n_sim = 100000,
                         sim_rfc = sim_rfc_AS_DTL
)
#
## 2. lambda.all=c(0.03, 0.05, 0.06), p.all=c(0.5, 0.3, 0.2)
p_AS_DTL_2 <- prob_cal_power(n = c(150, 150),
                         lambda.all=c(0.03, 0.05, 0.06),  # c(treatment1, treatment2, control)
                         p.all=c(0.5, 0.3, 0.2),               # c(treatment1, treatment2, control)
                         boundaries = c(0, c1, c2),
                         clu = clu,
                         n_sim = 100000,
                         sim_rfc = sim_rfc_AS_DTL
)
#
## 3. lambda.all=c(0.03, 0.04, 0.06), p.all=c(0.5, 0.4, 0.2)
p_AS_DTL_3 <- prob_cal_power(n = c(150, 150),
                         lambda.all=c(0.03, 0.04, 0.06),  # c(treatment1, treatment2, control)
                         p.all=c(0.5, 0.4, 0.2),               # c(treatment1, treatment2, control)
                         boundaries = c(0, c1, c2),
                         clu = clu,
                         n_sim = 100000,
                         sim_rfc = sim_rfc_AS_DTL
)
#
## 4. lambda.all=c(0.03, 0.03, 0.06), p.all=c(0.5, 0.5, 0.2) - global alternative hypothesis
p_AS_DTL_4 <- prob_cal_power(n = c(150, 150),
                         lambda.all=c(0.03, 0.03, 0.06),  # c(treatment1, treatment2, control)
                         p.all=c(0.5, 0.5, 0.2),               # c(treatment1, treatment2, control)
                         boundaries = c(0, c1, c2),
                         clu = clu,
                         n_sim = 100000,
                         sim_rfc = sim_rfc_AS_DTL
)
#
## 5. lambda.all=c(0.04, 0.06, 0.06), p.all=c(0.4, 0.2, 0.2)
p_AS_DTL_5 <- prob_cal_power(n = c(150, 150),
                         lambda.all=c(0.04, 0.06, 0.06),  # c(treatment1, treatment2, control)
                         p.all=c(0.4, 0.2, 0.2),               # c(treatment1, treatment2, control)
                         boundaries = c(0, c1, c2),
                         clu = clu,
                         n_sim = 100000,
                         sim_rfc = sim_rfc_AS_DTL
)
#
## 6. lambda.all=c(0.04, 0.05, 0.06), p.all=c(0.4, 0.3, 0.2)
p_AS_DTL_6 <- prob_cal_power(n = c(150, 150),
                         lambda.all=c(0.04, 0.05, 0.06),  # c(treatment1, treatment2, control)
                         p.all=c(0.4, 0.3, 0.2),               # c(treatment1, treatment2, control)
                         boundaries = c(0, c1, c2),
                         clu = clu,
                         n_sim = 100000,
                         sim_rfc = sim_rfc_AS_DTL
)
#
## 7. lambda.all=c(0.04, 0.04, 0.06), p.all=c(0.4, 0.4, 0.2)
p_AS_DTL_7 <- prob_cal_power(n = c(150, 150),
                         lambda.all=c(0.04, 0.04, 0.06),  # c(treatment1, treatment2, control)
                         p.all=c(0.4, 0.4, 0.2),               # c(treatment1, treatment2, control)
                         boundaries = c(0, c1, c2),
                         clu = clu,
                         n_sim = 100000,
                         sim_rfc = sim_rfc_AS_DTL
)
#
## 8. lambda.all=c(0.05, 0.06, 0.06), p.all=c(0.3, 0.2, 0.2)
p_AS_DTL_8 <- prob_cal_power(n = c(150, 150),
                         lambda.all=c(0.05, 0.06, 0.06),  # c(treatment1, treatment2, control)
                         p.all=c(0.3, 0.2, 0.2),               # c(treatment1, treatment2, control)
                         boundaries = c(0, c1, c2),
                         clu = clu,
                         n_sim = 100000,
                         sim_rfc = sim_rfc_AS_DTL
)
#
## 9. lambda.all=c(0.05, 0.05, 0.06), p.all=c(0.3, 0.3, 0.2)
p_AS_DTL_9 <- prob_cal_power(n = c(150, 150),
                         lambda.all=c(0.05, 0.05, 0.06),  # c(treatment1, treatment2, control)
                         p.all=c(0.3, 0.3, 0.2),               # c(treatment1, treatment2, control)
                         boundaries = c(0, c1, c2),
                         clu = clu,
                         n_sim = 100000,
                         sim_rfc = sim_rfc_AS_DTL
)

# aggregate data
table_AS_DTL <- do.call(rbind, mget(paste0("p_AS_DTL_", 1:9)))

end_time <- Sys.time()
cat("\nRun completed in", round(difftime(end_time, start_time, units = "mins"), 2), "minutes.\n")
print(table_AS_DTL)
#
#
table_probs <- cbind(table_MJW, table_MJW_DTL, table_AS, table_AS_DTL)
saveRDS(table_probs, file = "TABLE2_DATA.rds")

