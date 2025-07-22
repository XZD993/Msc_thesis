# Msc thesis:
# A group sequential three-arm two-stage trial design with rank-based endpoints
###################### Code for section 9.4 ####################################
################################################################################
#
#
source("FUNCTIONS.R")
#
#
# ----------------------- Global settings --------------------------------------
#
Sys.setenv(OPENBLAS_NUM_THREADS = 4)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")
#
n_cores <- 48
#
# ------------------ helper function -------------------------------------------
#
#
findBoundary_alphas_optim_futility <- function(n.seq = c(150, 150), 
                                               c.futility,
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
    boundaries <- c(c.futility, par[1], par[2])
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
                       boundaries = c(c.futility, final_c1, final_c2),
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
power_wrapper_futilities <- function(n_total,
                                     ratio_stage1,
                                     lambda.all,
                                     p.all,
                                     alphas_fun,
                                     alpha,
                                     side,
                                     c.start,
                                     c.futility,
                                     clu,
                                     n_sim,
                                     sim_rfc) {
  n1 <- floor(n_total * ratio_stage1)
  n2 <- n_total - n1
  n.seq <- c(n1, n2)
  
  # 1. Optimize boundaries for this n.seq to maintain 5% alpha level
  boundary_res <- findBoundary_alphas_optim_futility(n.seq = n.seq,
                                                     c.futility = c.futility,
                                                     c.start = c.start,
                                                     lambda = rep(0.06, 3),
                                                     p = rep(0.2, 3),
                                                     alphas_fun = alphas_fun,
                                                     alpha = alpha,
                                                     side = side,
                                                     clu = clu,
                                                     n_sim = n_sim,
                                                     sim_rfc = sim_rfc)
  gc()
  #if (boundary_res$converged != 0) {
  #  warning("Boundary optimization failed to converge.")
  #  return(c(Power = NA, Select_arm1 = NA, Select_both = NA, c1 = NA, c2 = NA))
  #}
  
  boundaries <- c(c.futility, boundary_res$c1, boundary_res$c2)
  
  # 2. Evaluate power with true effect size
  power_result <- prob_cal_power(n = n.seq,
                                 lambda.all = lambda.all,
                                 p.all = p.all,
                                 boundaries = boundaries,
                                 clu = clu,
                                 n_sim = n_sim,
                                 sim_rfc = sim_rfc)
  gc()
  
  c(Power = power_result[["Power"]],
    #Select_arm1 = power_result["Select_arm1"],
    #Select_both = power_result["Select_both"],
    c0 = c.futility,
    c1 = boundary_res$c1,
    c2 = boundary_res$c2)
}
#
#
# -------------------- Main sample size calculator -----------------------------
find_sample_size_futilities <- function(power_target = 0.8,
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
                                        c.futility = 0,
                                        clu,
                                        n_sim = 100000,
                                        sim_rfc
) {
  iter <- 1
  best_result <- NULL
  
  while (n_max - n_min > 2 && iter <= max_iter) {
    n_mid <- floor((n_min + n_max) / 2)
    
    res <- power_wrapper_futilities(n_total = n_mid,
                                    ratio_stage1 = ratio_stage1,
                                    lambda.all = lambda.all,
                                    p.all = p.all,
                                    alphas_fun = alphas_fun,
                                    alpha = alpha,
                                    side = side,
                                    c.start = c.start,
                                    c.futility = c.futility,
                                    clu = clu,
                                    n_sim = n_sim,
                                    sim_rfc = sim_rfc)
    
    if (is.na(res[["Power"]])) {
      stop("Simulation or boundary optimization returned NA.")
    }
    
    if (verbose) {
      cat(sprintf("n = %d → Power = %.3f, c0 = %.2f, c1 = %.2f, c2 = %.2f\n",
                  n_mid, res[["Power"]], res["c0"], res["c1"], res["c2"]))
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
    # clean up
    rm(res)
    gc()
    
    iter <- iter + 1
  }
  
  return(list(n_total = n_mid,
              power = best_result[["Power"]],
              #select_arm1 = final["Select_arm1"],
              #select_both = final["Select_both"],
              c0 = best_result[["c0"]],
              c1 = best_result[["c1"]],
              c2 = best_result[["c2"]],
              within_tol = abs(best_result[["Power"]] - power_target) <= tol))
}
#
#
#
# ----- sample size for different futilities in different stopping rules  ------
#
### General setting:
# 1) alpha: 5% (one-sided)
# 2) target power: 80%
# 3) LFC hypothesis: 
# lambda.all = c(0.04, 0.06, 0.06), p.all = c(0.4, 0.2, 0.2),
# 4) Futility boundaries: -2, -1.5, -1, -0.5, 0, 0.5
#
#
## cluster setting
start_time_ss_fut <- Sys.time()
clu <- makeCluster(n_cores)
on.exit(stopCluster(clu)) # auto shutdown cluster if error
clusterEvalQ(clu, source("FUNCTIONS.R"))
#
#
### -------------------  MJW stopping rule  ------------------------------------
## c0 = -2
ss_fut_MJW_minus2 <- find_sample_size_futilities(c.futility = -2, clu = clu, sim_rfc = sim_rfc_MJW)
#
## c0 = -1.5
ss_fut_MJW_minus1.5 <- find_sample_size_futilities(c.futility = -1.5, clu = clu, sim_rfc = sim_rfc_MJW)
#
## c0 = -1
ss_fut_MJW_minus1 <- find_sample_size_futilities(c.futility = -1, clu = clu, sim_rfc = sim_rfc_MJW)
#
## c0 = -0.5
ss_fut_MJW_minus0.5 <- find_sample_size_futilities(c.futility = -0.5, clu = clu, sim_rfc = sim_rfc_MJW)
#
## c0 = 0
ss_fut_MJW_0 <- find_sample_size_futilities(c.futility = 0, clu = clu, sim_rfc = sim_rfc_MJW)
#
## c0 = 0.5
ss_fut_MJW_0.5 <- find_sample_size_futilities(c.futility = 0.5, clu = clu, sim_rfc = sim_rfc_MJW)
#
## c0 = 1
ss_fut_MJW_1 <- find_sample_size_futilities(c.futility = 1, clu = clu, sim_rfc = sim_rfc_MJW)
#
#
### ------------------- MJW-DTL stopping rule  ---------------------------------
## c0 = -2
ss_fut_MJW_DTL_minus2 <- find_sample_size_futilities(c.futility = -2, clu = clu, sim_rfc = sim_rfc_MJW_DTL)
#
## c0 = -1.5
ss_fut_MJW_DTL_minus1.5 <- find_sample_size_futilities(c.futility = -1.5, clu = clu, sim_rfc = sim_rfc_MJW_DTL)
#
## c0 = -1
ss_fut_MJW_DTL_minus1 <- find_sample_size_futilities(c.futility = -1, clu = clu, sim_rfc = sim_rfc_MJW_DTL)
#
## c0 = -0.5
ss_fut_MJW_DTL_minus0.5 <- find_sample_size_futilities(c.futility = -0.5, clu = clu, sim_rfc = sim_rfc_MJW_DTL)
#
## c0 = 0
ss_fut_MJW_DTL_0 <- find_sample_size_futilities(c.futility = 0, clu = clu, sim_rfc = sim_rfc_MJW_DTL)
#
## c0 = 0.5
ss_fut_MJW_DTL_0.5 <- find_sample_size_futilities(c.futility = 0.5, clu = clu, sim_rfc = sim_rfc_MJW_DTL)
#
## c0 = 1
ss_fut_MJW_DTL_1 <- find_sample_size_futilities(c.futility = 1, clu = clu, sim_rfc = sim_rfc_MJW_DTL)
#
#
### ------------------------ AS stopping rule  ---------------------------------
## c0 = -2
ss_fut_AS_minus2 <- find_sample_size_futilities(c.futility = -2, clu = clu, sim_rfc = sim_rfc_AS)
#
## c0 = -1.5
ss_fut_AS_minus1.5 <- find_sample_size_futilities(c.futility = -1.5, clu = clu, sim_rfc = sim_rfc_AS)
#
## c0 = -1
ss_fut_AS_minus1 <- find_sample_size_futilities(c.futility = -1, clu = clu, sim_rfc = sim_rfc_AS)
#
## c0 = -0.5
ss_fut_AS_minus0.5 <- find_sample_size_futilities(c.futility = -0.5, clu = clu, sim_rfc = sim_rfc_AS)
#
## c0 = 0
ss_fut_AS_0 <- find_sample_size_futilities(c.futility = 0, clu = clu, sim_rfc = sim_rfc_AS)
#
## c0 = 0.5
ss_fut_AS_0.5 <- find_sample_size_futilities(c.futility = 0.5, clu = clu, sim_rfc = sim_rfc_AS)
#
## c0 = 1
ss_fut_AS_1 <- find_sample_size_futilities(c.futility = 1, clu = clu, sim_rfc = sim_rfc_AS)
#
#
### ------------------------ AS-DTL stopping rule  -----------------------------
## c0 = -2
ss_fut_AS_DTL_minus2 <- find_sample_size_futilities(c.futility = -2, clu = clu, sim_rfc = sim_rfc_AS_DTL)
#
## c0 = -1.5
ss_fut_AS_DTL_minus1.5 <- find_sample_size_futilities(c.futility = -1.5, clu = clu, sim_rfc = sim_rfc_AS_DTL)
#
## c0 = -1
ss_fut_AS_DTL_minus1 <- find_sample_size_futilities(c.futility = -1, clu = clu, sim_rfc = sim_rfc_AS_DTL)
#
## c0 = -0.5
ss_fut_AS_DTL_minus0.5 <- find_sample_size_futilities(c.futility = -0.5, clu = clu, sim_rfc = sim_rfc_AS_DTL)
#
## c0 = 0
ss_fut_AS_DTL_0 <- find_sample_size_futilities(c.futility = 0, clu = clu, sim_rfc = sim_rfc_AS_DTL)
#
## c0 = 0.5
ss_fut_AS_DTL_0.5 <- find_sample_size_futilities(c.futility = 0.5, clu = clu, sim_rfc = sim_rfc_AS_DTL)
#
## c0 = 1
ss_fut_AS_DTL_1 <- find_sample_size_futilities(c.futility = 1, clu = clu, sim_rfc = sim_rfc_AS_DTL)
#
#
#
#stopCluster(clu)
end_time_ss_fut <- Sys.time()
cat("\nRun completed in", round(difftime(end_time_ss_fut, start_time_ss_fut, units = "mins"), 2), "minutes.\n")
#
#
### ------------------------ Process outputs -----------------------------------
## --------- MJW
# Named list
ss_fut_list <- list(
  minus2 = ss_fut_MJW_minus2,
  minus1.5 = ss_fut_MJW_minus1.5,
  minus1 = ss_fut_MJW_minus1,
  minus0.5 = ss_fut_MJW_minus0.5,
  zero = ss_fut_MJW_0,
  plus0.5 = ss_fut_MJW_0.5,
  plus1 = ss_fut_MJW_1
)
# Apply cbind to each dataset in the list
ss_fut_tables <- lapply(ss_fut_list, function(x) do.call(cbind, x))
# Combine all into one table
ss_fut_table_MJW <- do.call(rbind, ss_fut_tables)
rownames(ss_fut_table_MJW) <- names(ss_fut_list)
#
#
## ----------- MJW-DTL
# Named list
ss_fut_list <- list(
  minus2 = ss_fut_MJW_DTL_minus2,
  minus1.5 = ss_fut_MJW_DTL_minus1.5,
  minus1 = ss_fut_MJW_DTL_minus1,
  minus0.5 = ss_fut_MJW_DTL_minus0.5,
  zero = ss_fut_MJW_DTL_0,
  plus0.5 = ss_fut_MJW_DTL_0.5,
  plus1 = ss_fut_MJW_DTL_1
)
# Apply cbind to each dataset in the list
ss_fut_tables <- lapply(ss_fut_list, function(x) do.call(cbind, x))
# Combine all into one table
ss_fut_table_MJW_DTL <- do.call(rbind, ss_fut_tables)
rownames(ss_fut_table_MJW_DTL) <- names(ss_fut_list)
#
#
### ------------- AS
# Named list
ss_fut_list <- list(
  minus2 = ss_fut_AS_minus2,
  minus1.5 = ss_fut_AS_minus1.5,
  minus1 = ss_fut_AS_minus1,
  minus0.5 = ss_fut_AS_minus0.5,
  zero = ss_fut_AS_0,
  plus0.5 = ss_fut_AS_0.5,
  plus1 = ss_fut_AS_1
)
# Apply cbind to each dataset in the list
ss_fut_tables <- lapply(ss_fut_list, function(x) do.call(cbind, x))
# Combine all into one table
ss_fut_table_AS <- do.call(rbind, ss_fut_tables)
rownames(ss_fut_table_AS) <- names(ss_fut_list)
#
#
### ------------- AS-DTL
# Named list
ss_fut_list <- list(
  minus2 = ss_fut_AS_DTL_minus2,
  minus1.5 = ss_fut_AS_DTL_minus1.5,
  minus1 = ss_fut_AS_DTL_minus1,
  minus0.5 = ss_fut_AS_DTL_minus0.5,
  zero = ss_fut_AS_DTL_0,
  plus0.5 = ss_fut_AS_DTL_0.5,
  plus1 = ss_fut_AS_DTL_1
)
# Apply cbind to each dataset in the list
ss_fut_tables <- lapply(ss_fut_list, function(x) do.call(cbind, x))
# Combine all into one table
ss_fut_table_AS_DTL <- do.call(rbind, ss_fut_tables)
rownames(ss_fut_table_AS_DTL) <- names(ss_fut_list)
#
#
### Save results
ss_fut_table <- cbind(as.data.frame(ss_fut_table_MJW)[,1:2],
                      as.data.frame(ss_fut_table_MJW_DTL)[,1:2],
                      as.data.frame(ss_fut_table_AS)[,1:2],
                      as.data.frame(ss_fut_table_AS_DTL)[,1:2])
saveRDS(ss_fut_table, file = "TABLE3_DATA.rds")



