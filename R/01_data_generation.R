rm(list = ls())
source("./R/header.R")

# func - generate_simulation_data ----------------------------------------------
generate_simulation_data <- function(n = 100, p = 15) {
  # 1. generate complete
  lambda <- c(50.0, 48.0, 45.0, 25.0, 20.0, 10.0, 5.0,
              5.0, 1.0, 1.0, 0.5, 0.5, 0.5, 0.1, 0.1)
  
  sigma <- eiginv(evals = lambda, symmetric = TRUE)
  R <- chol(sigma)
  Z <- matrix(rnorm(n * p), nrow = n, ncol = p)
  X <- Z %*% R
  
  # 2. create missings
  # MCAR small
  X_mcar_small <- X
  X_mcar_small[matrix(runif(n*p), n, p) < 0.05] <- NA
  
  # MCAR large
  X_mcar_large <- X
  X_mcar_large[matrix(runif(n*p), n, p) < 0.30] <- NA
  
  # MAR small
  X_mar_small <- mice::ampute(data = X, 
                              prop = 0.05, 
                              mech = "MAR", 
                              bycases = FALSE)$amp
  
  # MAR large
  pattern_mat <- mice::ampute(data = X, prop = 0.05, mech = "MAR", 
                              bycases = FALSE, run = FALSE)$patterns
  
  for (j in 1:p) { 
    pattern_mat[j, sample(1:p, sample(5:10, 1))] <- 0 
  }
  
  X_mar_large <- mice::ampute(data = X, 
                              prop = 0.30, 
                              mech = "MAR", 
                              bycases = FALSE, 
                              patterns = pattern_mat)$amp
  
  list(true_data = X,
       sigma_true = sigma,
       settings = list(mcar_small = X_mcar_small,
                       mcar_large = X_mcar_large,      
                       mar_small  = X_mar_small,      
                       mar_large  = X_mar_large))
}


# data generation --------------------------------------------------------------
sim_data <- map(1:1000, ~generate_simulation_data(n = 100, p = 15))


# local save of results --------------------------------------------------------
save(sim_data, file = "./R/data_and_results/1000_sim_data.RData")
