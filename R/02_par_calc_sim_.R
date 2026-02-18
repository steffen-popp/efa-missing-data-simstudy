# 1. load simulated data and functions -----------------------------------------
rm(list = ls())
source("./R/header.R")
source("./R/function_covariances.R")
load("./R/data_and_results/1000_sim_data.RData")

# 2. specify current setting ---------------------------------------------------
# possible values: "mcar_small", "mcar_large", "mar_small", "mar_large"
current_set <- "mar_large" 

# parameters for covariance function (m and ci)
setting_config <- list(mcar_small = list(m = 10, ci = FALSE),
                       mcar_large = list(m = 30, ci = FALSE),
                       mar_small  = list(m = 10, ci = TRUE),
                       mar_large  = list(m = 30, ci = TRUE))

# 3. parallel calculation ------------------------------------------------------
plan(multisession, workers = parallel::detectCores() - 1)

handlers(handler_pbcol(adjust = 1.0,
                       complete = function(s) cli::bg_red(cli::col_black(s)),
                       incomplete = function(s) cli::bg_cyan(cli::col_black(s))
))

# start calculation
results_this_setting <- with_progress({
  p <- progressor(steps = length(sim_data))
  
  # map over 1000 runs in list
  future_map(sim_data, function(run_data) {
    p()
    
    # extract die parameter for actual setting
    conf <- setting_config[[current_set]]
    
    # call compute_covariances function with the appropriate parameters
    result <- compute_covariances(X_comp = run_data$true_data,
                                  X_miss = run_data$settings[[current_set]],
                                  mi_m   = conf$m,
                                  run_ci = conf$ci)
    
    # add true covariance-matrix to result object
    result$sigma_true <- run_data$sigma_true
    return(result)
    
  }, .options = furrr_options(seed = TRUE))
})

# 4. save results --------------------------------------------------------------
saveRDS(results_this_setting, 
        file = paste0("./R/data_and_results/results_", current_set, ".rds"))

# 5. set plan to sequential ----------------------------------------------------
plan(sequential)