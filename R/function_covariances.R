# function to calculate the covariance matrix ----------------------------------
compute_covariances <- function(X_comp, X_miss, mi_m = 10, run_ci= FALSE) {
  
  # calculate true covariance matrix
  cov_true <- cov(X_comp)

  # listwise deletion covariance
  cov_listwise <- tryCatch(cov(X_miss, use = "complete.obs"),
                           error = function(e) matrix(NA, 
                                                      ncol = ncol(X_miss), 
                                                      nrow = ncol(X_miss)))

  # pairwise deletion covariance
  cov_pairwise <- cov(X_miss, use = "pairwise.complete.obs")

  # IPCA-reg imputed covariance matrices
  ipca_reg <- missMDA::imputePCA(X_miss,
                                 method = "Regularized",
                                 coeff.ridge = 1, # regularization parameter
                                 maxiter = 1000) # maximum nr of iterations
  cov_ipca_reg <- ipca_reg$completeObs %>% cov()

  
  # Multiple Imputation covariance matrix
  cov_mi <- mifa::mifa(as.data.frame(X_miss),
                       m = mi_m, # for small amount of missing 10, for large 30
                       maxit = 5,
                       n_boot = 500, 
                       ci = if(run_ci) "both" else FALSE,
                       printFlag = FALSE
                       )
  # uses pmm as default imputation method

  # return as list
  list(covariance = list(true = cov_true,
                         listwise = cov_listwise,
                         pairwise = cov_pairwise,
                         ipca_reg = cov_ipca_reg,
                         mi = cov_mi$cov_combined),
       imputed_data = list(ipca_reg = ipca_reg,
                           mi = cov_mi))

}
