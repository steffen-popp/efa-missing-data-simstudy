# 1. load data, results, and functions -----------------------------------------
rm(list = ls())
source("./R/header.R")
source("./R/function_mahalanobis_measures.R")
source("./R/02_par_process_results.R")

# initiate variables for processing
settings <- names(covariances_by_setting)
methods <- c("listwise", "pairwise", "ipca_reg", "mi")



# 2 distance measures for covariance matrices ----------------------------------
## 2.1 mahalanobis-based distance ----------------------------------------------
# function to calculate the mahalanobis-based distance
calculate_mahal_by_setting <- function(cov_list, setting_name, N, methods){
  # imap_dfr iterates over cov_list and binds the results into a single df
  # .x refers to the content of the list element, .y refers to the simrun
  imap_dfr(cov_list, ~{
    S_true <- .x$true
    
    # iterate over each method specified in the 'methods' vector
    map_dfr(methods, function(method){
      S_curr <- .x[[method]]
      # check if matrix has valid entries (no NAs or Infs)
      # if NAs exist, its not positiv definit
      if (any(is.na(S_curr)) || any(is.infinite(S_curr))) {
        is_pd <- FALSE
      } else {
        # calculate eigenvalues (use try catch to handle potential errors)
        ev <- tryCatch(eigen(S_curr, only.values = TRUE)$values, 
                       error = function(e) return(-1))
        is_pd <- all(ev > 0) # check if all eigenvalues are positive
      }
      
      # calculate distance measure only if matrix is positiv definit
      dist_val <- if(is_pd) {
        calculate_mahalanobis_distance(S_comp = S_true, S_miss = S_curr, N = N)
      } else {
        NA
      }
      
      # return table with results
      tibble(run = .y,
             method = method,
             value = dist_val,
             setting = setting_name)
    })
  })
}

# calculation of all distances
mahala_nobis_distances_all <- imap_dfr(covariances_by_setting, ~{
  calculate_mahal_by_setting(.x, .y, N = 100, methods = methods)
})



## 2.2 Frobenius-Norm ----------------------------------------------------------
calc_frob <- function(S1, S2) {
  sqrt(sum((S1 - S2)^2))
}

# calculation for all settings
# imap_dfr is used to iterate over settings and combine results into a single df
frob_distances_all <- imap_dfr(covariances_by_setting, 
                               function(cov_list, setting_name) {
  # loop through each individual simulation run within the current setting
  imap_dfr(cov_list, function(run_data, run_id) {
    S_true <- run_data$true # Extract the truth covariance matrix
    
    # iterate through the list of methods to compare against the truth
    map_dfr(methods, function(m) {
      S_miss <- run_data[[m]] # extract the estimated matrix for method 'm'
      
      # check for NAs, if any exist, calculation isn't possible
      dist <- if(any(is.na(S_miss))) NA else calc_frob(S_true, S_miss)
      
      # return table with results
      tibble(run = run_id,
             method = m,
             value = dist,
             setting = setting_name)
    })
  }) 
})


# 4. measures about eigenvalues ------------------------------------------------
# create container list
eigenvalues_all <- list()

for(i in settings){
  # extract data for the current simulation setting
  cov_data <- covariances_by_setting[[i]]
  
  ## i) calculate eigenvalues --------------------------------------------------
  eigenvalues <- imap(cov_data, ~ {
    tibble(simrun = .y, 
           !!!map(.x, function(x){ # !!! unpacks the list of matrices
             tryCatch(eigen(x, symmetric = TRUE)$values,
                      error = function(e) rep(NA, ncol(x)))
           }))}) %>% 
    imap(~ {.x %>% mutate(k = row_number(), .before = 1)}) # adds rownumber (k)
  
  ## ia) create a lookup table to track which method/run resulted in a pd matrix
  is_pd_df <- imap_dfr(eigenvalues, ~ {
    # check columns (excluding 'k' and 'simrun') for positive eigenvalues
    map_dfc(.x[-c(1,2)], function(ev_vector) {
      all(ev_vector > 0, na.rm = TRUE)}) %>% 
      mutate(simrun = .y)}) %>% 
    # reshape to long format
    pivot_longer(-simrun, names_to = "method", values_to = "is_pd")
  
  ## ii) proportion of non-pd matrices -----------------------------------------
  # check for negative eigenvalues in covariance matrices 
  prop_pd_matrices <- imap_dfr(eigenvalues, ~ {
    tibble(simrun = .y, !!!map_lgl(.x[-c(1,2)], ~ any(.x < 0)))}) %>% 
    select(-simrun) %>% 
    colMeans() 
  
  ## iii) Differences between prop. of explained variances ---------------------
  # Calculate cumulative proportion of explained variance for each method
  eigen_shares <- eigenvalues %>% 
    bind_rows() %>% 
    group_by(simrun) %>% 
    mutate(across(c(true, listwise, pairwise, ipca_reg, mi),
                  ~ cumsum(.x) / sum(.x)))
  
  # calculate differences to true proportions
  eigen_share_diff <- eigen_shares %>% 
    mutate(across(c(mi, ipca_reg, pairwise, listwise), ~ .x - true))
  
  # reshape for plotting to long format
  eigen_share_diff_long <- eigen_share_diff %>%
    pivot_longer(cols = c(mi, ipca_reg, pairwise, listwise),
                 names_to = "method",
                 values_to = "gamma_diff") %>% 
    left_join(is_pd_df, by = c("simrun", "method")) %>% # to filter pd matrices
    filter(is_pd == TRUE) %>% 
    mutate(method = recode(method, mi = "MI", ipca_reg = "IPCA-reg",
                           pairwise = "Pairwise", listwise = "Listwise"),
           method = factor(method, levels = c("Listwise", "Pairwise", 
                                              "IPCA-reg", "MI")))
  
  
  ## iv) save data -------------------------------------------------------------
  eigenvalues_all[[i]] <- list(
    eigenvalue = eigenvalues,
    eigenvalue_share = eigen_shares,
    eigenvalue_share_diff = eigen_share_diff,
    eigenvalue_share_diff_long = eigen_share_diff_long,
    proportion_pd_matrices = prop_pd_matrices)

}



# 3. things about confidence intervals and coverage ----------------------------
## 3.1 calculate confidence interval -------------------------------------------
df_all_confint <- map_df(conf_int, ~bind_rows(.x, .id = "run"), 
                         .id = "condition")


# 2. calculate averaged ci's
conf_int_averaged <- df_all_confint %>%
  group_by(condition, n_pc) %>%
  summarise(
    var_explained = mean(var_explained, na.rm = TRUE),
    ci_boot_lower = mean(ci_boot_lower, na.rm = TRUE),
    ci_boot_upper = mean(ci_boot_upper, na.rm = TRUE),
    ci_fieller_lower = mean(ci_fieller_lower, na.rm = TRUE),
    ci_fieller_upper = mean(ci_fieller_upper, na.rm = TRUE),
    width_boot = mean(ci_boot_upper - ci_boot_lower, na.rm = TRUE),
    width_fieller = mean(ci_fieller_upper - ci_fieller_lower, na.rm = TRUE),
    .groups = "drop"
  )

conf_int_averaged %>% filter(condition == "mar_small")
conf_int_averaged %>% filter(condition == "mar_large")


# 3.2 calculate coverage -------------------------------------------------------
true_shares_df <- map_dfr(eigenvalues_all, ~{.x$eigenvalue_share %>% 
    select(simrun, k, true)}, .id = "condition")

coverage_data <- df_all_confint %>%
  mutate(run = as.integer(run)) %>% 
  # joinen by condition (=setting), run and n_pc (=k)
  left_join(true_shares_df, by = c("condition", "run" = "simrun", "n_pc" = "k"))

coverage <- coverage_data %>%
  group_by(condition, n_pc) %>%
  summarise(coverage_boot = mean(true >= ci_boot_lower & 
                                   true <= ci_boot_upper, 
                                 na.rm = TRUE),
            coverage_fieller = mean(true >= ci_fieller_lower & 
                                      true <= ci_fieller_upper, 
                                    na.rm = TRUE),
            rmse = sqrt(mean((var_explained - true)^2, na.rm = TRUE)),
            rrmse = rmse / mean(true, na.rm = TRUE),
            .groups = "drop")

coverage %>% filter(condition == "mar_small")
coverage %>% filter(condition == "mar_large")


# local save -------------------------------------------------------------------
# save(eigenvalues_all, 
#      file = "./R/data_and_results/par_eigenvalues.RData")
# save(mahala_nobis_distances_all, 
#      file = "./R/data_and_results/par_mahala_nobis.RData")
# save(frob_distances_all, 
#      file = "./R/data_and_results/par_frob_distances.RData")
# save(coverage, file = "./R/data_and_results/par_coverage.RData")
# save(conf_int_averaged, file = "./R/data_and_results/par_averaged_ci.RData")