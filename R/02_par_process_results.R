load("./R/data_and_results/1000_sim_data.RData")

# 1. reformat simulation data --------------------------------------------------
data <- lapply(sim_data, function(iteration) {
  list(X = iteration$true_data,
       X_mcar_small = iteration$settings$mcar_small,
       X_mcar_large = iteration$settings$mcar_large,
       X_mar_small  = iteration$settings$mar_small,
       X_mar_large  = iteration$settings$mar_large)
})


# 2. load and pre-process results ----------------------------------------------
mcar_small <- readRDS("./R/data_and_results/results_mcar_small.rds") %>% 
  list_transpose()
mcar_large <- readRDS("./R/data_and_results/results_mcar_large.rds") %>% 
  list_transpose()
mar_small <- readRDS("./R/data_and_results/results_mar_small.rds") %>% 
  list_transpose()
mar_large <- readRDS("./R/data_and_results/results_mar_large.rds") %>% 
  list_transpose()

covariances_by_setting <- list(mcar_small = mcar_small$covariance,
                               mcar_large = mcar_large$covariance,
                               mar_small = mar_small$covariance,
                               mar_large = mar_large$covariance)


# 3. pre-process data for confidence intervals ---------------------------------
conf_int <- list(
  mar_small = lapply(mar_small$imputed_data, function(x) x$mi$var_explained),
  mar_large = lapply(mar_large$imputed_data, function(x) x$mi$var_explained))

