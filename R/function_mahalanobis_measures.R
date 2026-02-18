# function to calculate mahalanobis distance -----------------------------------
calculate_mahalanobis_distance <- function(S_comp, S_miss, N) {

  # 1. calculate den difference and generate vech
  delta <- vech(S_miss - S_comp)
  
  # 2. calculate kronecker product
  S_kron_S <- S_comp %x% S_comp
  
  # 3. create elimination matrix H and calculate Variance of vech(S_comp)
  H <- elimination.matrix(n = ncol(S_comp))
  Var_vech_Scomp <- 2 * (N - 1) * (H %*% S_kron_S %*% t(H))
  
  # 4. calculate inverse of Var_vech_Scomp
  # catch errors, e.g. if Var_vech_Scomp is singular
  Var_inv <- tryCatch(expr = solve(Var_vech_Scomp),
                      error = function(e) {
    warning("Varianzmatrix ist singulär oder fehlerhaft. Kann Mahalanobis-Distanz nicht berechnen.")
    return(NULL)
  })
  
  if (is.null(Var_inv)) {
    return(NA)
    }
  
  # 5. finally calculate the distance
  distance <- sqrt(t(delta) %*% Var_inv %*% delta)
  
  # 6. return the result
  as.numeric(distance)
}