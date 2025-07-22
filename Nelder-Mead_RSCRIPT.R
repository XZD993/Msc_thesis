###################### Simulation and Data processing ##########################
#
Sys.setenv(OPENBLAS_NUM_THREADS = 4)
Sys.getenv("OPENBLAS_NUM_THREADS")
Sys.setenv(OMP_NUM_THREADS = 1)
Sys.getenv("OMP_NUM_THREADS")
source("FUNCTIONS.R")
# -----------------------Package loading----------------------------------------
#
library(parallel)
library(tidyverse)
library(dplyr)
library(gtsummary)
library(flextable)
library(plotly)
library(GGally)
#
# ----------------------- Global settings --------------------------------------
## Boundaries setting
#
# create all boundaries combinations
bo_vals <- seq(1, 3, by = 0.01)

boundary_grid <- list()
for (c1 in bo_vals) {
  for (c2 in bo_vals) {
    boundary_grid[[length(boundary_grid)+1]] <- c(0, c1, c2)
  }
}
#
## parameter setting of simulation
n_sim <- 10000
param_combos <- expand.grid(sim_id = 1:n_sim,
                            boundary_id = seq_along(boundary_grid))
n_cores <- 48
#
#
# ---------------------- Helper Simulation Function ----------------------------
## The simulation function for each batch. 
# It is like the bullet in a gun and we need a trigger to shoot it.
sim_bullet <- function(batch_id, param_combos_batch, lambda.all, hypothesis) {
  cat("Running batch", batch_id, "with", nrow(param_combos_batch), "simulations...\n")
  
  clu <- makeCluster(n_cores)
  clusterExport(clu, c("sim_rfc_MJW", "dat_gen_uni", "Wij_cal_uni",
                       "boundary_grid"))
  
  param_sweep <- function(row_id){
    combo <- param_combos_batch[row_id,]
    sim_res <- sim_rfc_MJW(i = combo$sim_id,
                           n=c(150, 150),
                           lambda.all = lambda.all,
                           boundaries = boundary_grid[[combo$boundary_id]])
    c(sim_res, boundary_id = combo$boundary_id)
  }
  
  #start.time <- Sys.time()
  results <- parLapply(clu, 1:nrow(param_combos_batch), param_sweep)
  stopCluster(clu)
  #end.time <- Sys.time()
  
  results_df <- as.data.frame(do.call(rbind, results))
  results_df$boundary <- sapply(results_df$boundary_id, function(x) {
    paste(boundary_grid[[x]], collapse = ",")
  })
  #results_df$run_time <- as.numeric(end.time - start.time, units = "secs")
  
  saveRDS(results_df, paste0(gsub("-", "", Sys.Date()), "_sim_", hypothesis, "_batch_", batch_id, ".rds"))
  
  # remove big objects
  rm(results_df, results) 
  # clear memory after each batch
  gc()
}
#
#
#
## Main simulation function, i.e. simulation trigger
# run simulation in batches
sim_trigger <- function(hypothesis = "null", lambda.all = c(1, 1, 1), n_batches = 80) {
  total_rows <- nrow(param_combos)
  batch_size <- ceiling(total_rows / n_batches)
  
  for (b in 1:n_batches) {
    start_id <- (b - 1) * batch_size + 1
    end_id <- min(b * batch_size, total_rows)
    param_combos_batch <- param_combos[start_id:end_id, ]
    
    sim_bullet(batch_id = b,
               param_combos_batch = param_combos_batch,
               lambda.all = lambda.all,
               hypothesis = hypothesis)
  }
}
#
# ---------------------------------- Simulation --------------------------
### Null hypothesis
## 1. boundaries: c1: 1.00-3.00, c2: 1.00-3.00 
start.time_null <- Sys.time()

sim_trigger(hypothesis = "null", lambda.all = c(1, 1, 1), n_batches = 80)

end.time_null <- Sys.time()
run.time_null <- as.numeric(end.time_null - start.time_null, units = "mins")
#
#
#
# -------------------- Results processing --------------------------------------
### 
#
### 1. Null hypothesis - 20250418 simulation null, boundaries: c1: 1.00-3.00, c2: 1.00-3.00
## data loading
# Due to huge amount of data. split the data into 4 parts
#
# part 1: 1-20
null_files. <- sprintf("/nfsmb/koll/zhidong.xie/MSC thesis/code/20250619_sim_null_batch_%d.rds", 1:20)
res_df_null. <- do.call(rbind, lapply(null_files., readRDS))
#
# get the result table
res_table_null_1 <- res_df_null. %>% 
  group_by(boundary) %>% 
  summarise(n.total = n(),
            p.stage1 = sum(Success.stage1 == 1)/n.total,
            p.stage2 = sum(Success.stage2 == 1)/n.total,
            p.global = sum(Success == 1)/n.total) %>% 
  ungroup() %>% 
  select(p.stage1, p.stage2, p.global, boundary)
#
# clean up in case of memory overload
rm(list = grep("\\.$", ls(all.names = T), value = T))
gc()
#
#
# part 2: 21-40
null_files. <- sprintf("/nfsmb/koll/zhidong.xie/MSC thesis/code/20250619_sim_null_batch_%d.rds", 21:40)
res_df_null. <- do.call(rbind, lapply(null_files., readRDS))
#
# get the result table
res_table_null_2 <- res_df_null. %>% 
  group_by(boundary) %>% 
  summarise(n.total = n(),
            p.stage1 = sum(Success.stage1 == 1)/n.total,
            p.stage2 = sum(Success.stage2 == 1)/n.total,
            p.global = sum(Success == 1)/n.total) %>% 
  ungroup() %>% 
  select(p.stage1, p.stage2, p.global, boundary)
#
# clean up in case of memory overload
rm(list = grep("\\.$", ls(all.names = T), value = T))
gc()
#
#
# part 3: 41-60
null_files. <- sprintf("/nfsmb/koll/zhidong.xie/MSC thesis/code/20250619_sim_null_batch_%d.rds", 41:60)
res_df_null. <- do.call(rbind, lapply(null_files., readRDS))
#
# get the result table
res_table_null_3 <- res_df_null. %>% 
  group_by(boundary) %>% 
  summarise(n.total = n(),
            p.stage1 = sum(Success.stage1 == 1)/n.total,
            p.stage2 = sum(Success.stage2 == 1)/n.total,
            p.global = sum(Success == 1)/n.total) %>% 
  ungroup() %>% 
  select(p.stage1, p.stage2, p.global, boundary)
#
# clean up in case of memory overload
rm(list = grep("\\.$", ls(all.names = T), value = T))
gc()
#
#
# part 4: 61-80
null_files. <- sprintf("/nfsmb/koll/zhidong.xie/MSC thesis/code/20250619_sim_null_batch_%d.rds", 61:80)
res_df_null. <- do.call(rbind, lapply(null_files., readRDS))
#
# get the result table
res_table_null_4 <- res_df_null. %>% 
  group_by(boundary) %>% 
  summarise(n.total = n(),
            p.stage1 = sum(Success.stage1 == 1)/n.total,
            p.stage2 = sum(Success.stage2 == 1)/n.total,
            p.global = sum(Success == 1)/n.total) %>% 
  ungroup() %>% 
  select(p.stage1, p.stage2, p.global, boundary)
#
# clean up in case of memory overload
rm(list = grep("\\.$", ls(all.names = T), value = T))
gc()
#
#
### aggregate dataset
#
#
res_table_null <- res_table_null_1 %>% 
  bind_rows(res_table_null_2, res_table_null_3, res_table_null_4) %>% 
  separate(col = boundary, into = c("c0", "c1", "c2"), sep = ",", convert = T)
#
# save result
saveRDS(res_table_null, "Data_Nelder-Mead")
#
#
#
#
# ------------------------ Figure 6 & 7 ----------------------------------------
##### alpha
## data loading
dat <- readRDS("Data_Nelder-Mead")
#
#
#
### Figure 6
# 3D plot
dat$sig <- ifelse(dat$p.global < 0.05, "< 5%", ">= 5%")

plot_ly(dat, x = ~c1, y = ~c2, z = ~p.global*100, 
        type = 'scatter3d', mode = 'markers', color = ~sig, colors = c("#89ABE3", "#EA738D"), 
        size = 1
        # marker = list(color = ifelse(res_table_null$sig == 1, "red", "blue"), size = 2)
) %>%
  layout(title = list(text = "Type I error under different boundary values"),
         scene = list(xaxis = list(title = 'c1'),
                      yaxis = list(title = 'c2'),
                      zaxis = list(title = 'alpha (%)')))
#
#
#
### Figure 7

#
# Loss function - differentiated by the value of loss function
## threshold 0.001
alphas <- alpha_spend(iuse = 1, alpha = 0.05, side = 1, t = c(0.5, 1))$pd
dat$alpha1 <- alphas[1]
dat$alpha2 <- alphas[2]
dat$loss <- (dat$p.stage1 - dat$alpha1)^2 + 
  (dat$p.stage2 - dat$alpha2)^2
dat$sig <- ifelse(dat$loss < 0.001, "loss <= 0.001", "loss > 0.001")

plot_ly(dat, x = ~c1, y = ~c2, z = ~loss, 
        type = 'scatter3d', mode = 'markers', color = ~sig, colors = c("#1a80bb", "#ea801c"),
        size = 1
        # marker = list(color = ifelse(res_table_null$sig == 1, "red", "blue"), size = 2)
) %>%
  layout(title = list(text = "Loss function values with boundaries from 1 to 3"),
         scene = list(xaxis = list(title = 'c1'),
                      yaxis = list(title = 'c2'),
                      zaxis = list(title = 'loss')))
#
#
## threshold 0.0001
alphas <- alpha_spend(iuse = 1, alpha = 0.05, side = 1, t = c(0.5, 1))$pd
dat$alpha1 <- alphas[1]
dat$alpha2 <- alphas[2]
dat$loss <- (dat$p.stage1 - dat$alpha1)^2 + 
  (dat$p.stage2 - dat$alpha2)^2
dat$sig <- ifelse(dat$loss < 0.0001, "loss <= 0.0001", "loss > 0.0001")

plot_ly(dat, x = ~c1, y = ~c2, z = ~loss, 
        type = 'scatter3d', mode = 'markers', color = ~sig, colors = c("#1a80bb", "#ea801c"),
        size = 1
        # marker = list(color = ifelse(res_table_null$sig == 1, "red", "blue"), size = 2)
) %>%
  layout(title = list(text = "Loss function values with boundaries from 1 to 3"),
         scene = list(xaxis = list(title = 'c1'),
                      yaxis = list(title = 'c2'),
                      zaxis = list(title = 'loss')))
#
#
## threshold 0.00001
alphas <- alpha_spend(iuse = 1, alpha = 0.05, side = 1, t = c(0.5, 1))$pd
dat$alpha1 <- alphas[1]
dat$alpha2 <- alphas[2]
dat$loss <- (dat$p.stage1 - dat$alpha1)^2 + 
  (dat$p.stage2 - dat$alpha2)^2
dat$sig <- ifelse(dat$loss < 0.00001, "loss <= 0.00001", "loss > 0.00001")

plot_ly(dat, x = ~c1, y = ~c2, z = ~loss, 
        type = 'scatter3d', mode = 'markers', color = ~sig, colors = c("#1a80bb", "#ea801c"),
        size = 1
        # marker = list(color = ifelse(res_table_null$sig == 1, "red", "blue"), size = 2)
) %>%
  layout(title = list(text = "Loss function values with boundaries from 1 to 3"),
         scene = list(xaxis = list(title = 'c1'),
                      yaxis = list(title = 'c2'),
                      zaxis = list(title = 'loss')))