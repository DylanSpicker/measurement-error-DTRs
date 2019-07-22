RC <- function(W, Z=NULL, return_var=FALSE) {
  # Compute the Regression Calibration imputed estimates of 'W'
  # Parameters
  #    W: Observed Data Matrix w/ Error + Replicates
  #       > Rows represent Observations; Columns Replicates
  #       > NA are permitted for samples w/o that replicate
  #       > K's are inferred from the matrix
  #    Z: Observed Covariate Matrix w/o Error:
  #       > NULL Z refers to the case wherein no other covariates
  #           are considered
  #       > Otherwise, rows represent observations
  #         each column is a covariate
  #       > Assumed that all samples are present (i.e no NA)
  #   return_var: Boolean
  #       > If true (default is false) then the estimated variance parameters
  #         are returned
  # Returns:
  #     X.imp: The imputed X values; will be Nx1 vector (return_var = false)
  #     ret_obj: Object with $X.imp, $Sigma.XX, $Sigma.UU, $Mu.W
  #              and $Sigma.XZ, $Sigma.ZZ, $Z.bar if Z is not NULL
  X.imp <- NULL
  
  n <- dim(W)[1] # Extract the Sample Size
  W.bar <- rowMeans(W, na.rm=T) # Compute the Mean; Remove NA's
  K <- apply(W, 1, function(x) sum(!is.na(x))) # Compute the Replicate Counts
  
  Mu.W <- sum(K*W.bar,na.rm=T)/sum(K) # Mean over W
  v <- sum(K) - sum(K^2)/sum(K)
  
  # Compute Sigma.UU
  # TODO: This can be more efficient, but for now it's fine
  Sigma.UU <- 0
  for(ii in 1:n) {
    Sigma.UU <- Sigma.UU + sum((W[ii,] - W.bar[ii])^2, na.rm=T)
  }
  Sigma.UU <- Sigma.UU/sum(K-1)
  
  Sigma.XX <- (sum(K*(W.bar-Mu.W)*(W.bar-Mu.W)) - (n-1)*Sigma.UU)/v
  
  if (! is.null(Z)) {
    # We do have Error-Free Covariates
    Z.bar <- colMeans(Z) # Means of Error-Free Covariates
    Sigma.ZZ <- cov(Z) # Covariance Matrix for Z
    
    # Initialize Cross Matrix
    n.Z <- dim(Z)[2]
    
    Sigma.XZ <- matrix(nrow=1, ncol=n.Z)
    
    # Loop over the Error Free Covariates
    for (ii in 1:n.Z) {
      Sigma.XZ[,ii] <- sum(K*(W.bar - Mu.W)%*%t(Z[,ii] - Z.bar[ii]))/v
    }
    
    for (ii in 1:n) {
      inv_mat <- rbind(cbind(as.matrix(Sigma.XX + Sigma.UU/K[ii]), Sigma.XZ),cbind(t(Sigma.XZ), Sigma.ZZ))
      if (is.singular.matrix(inv_mat)) {
        X.imp[ii] <- Mu.W + cbind(as.matrix(Sigma.XX), Sigma.XZ)%*%
                          pseudoinverse(inv_mat)%*%
                          t(cbind(W.bar[ii] - Mu.W, t(Z[ii,] - Z.bar)))  
      } else {
        X.imp[ii] <- Mu.W + cbind(as.matrix(Sigma.XX), Sigma.XZ)%*%
                          solve(inv_mat)%*%
                          t(cbind(W.bar[ii] - Mu.W, t(Z[ii,] - Z.bar)))
      }
      
    }
  } else {
    # No need for Matrix Computation, everything is a scalar
    for (ii in 1:n) {
      X.imp[ii] <- Mu.W + Sigma.XX*(W.bar[ii] - Mu.W)/(Sigma.XX + Sigma.UU/K[ii])  
    }
  }
  
  if (return_var) {
    ret_obj = list()
    ret_obj$X.imp <- X.imp
    ret_obj$Sigma.XX <- Sigma.XX
    ret_obj$Sigma.UU <- Sigma.UU
    if (! is.null(Z)){ 
      ret_obj$Sigma.XZ <- Sigma.XZ
      ret_obj$Z.bar <- Z.bar
      ret_obj$Sigma.ZZ <- Sigma.ZZ
    }
    ret_obj$Mu.W <- Mu.W
    return (ret_obj)
  }
  return (X.imp)
  
}