set.seed(314159265)
source("rCalibration.R") # Include file for Regression Calibration

# Libraries Used 
library(xtable)     # For outputting tabular results for Latex
library(foreach)    # For Parallelization of Simulations
library(doParallel) # To go along w/ foreach
library(matrixcalc) 
library(corpcor)

# User Defined Functions
expit <- function(w){return(1/(1+exp(-w)))}

## Simulation Controls:
##    Set the following variables to 'true' or 'false' based on which settings
##    you are interested in running.
sensitivity1 <- F # Sensitivity Analysis for the alpha parameters (1 stage)
sensitivity2 <- F # Sensitivity Analysis for heteroskedastic errors (pseudo outcomes)
sensitivity3 <- F # Sensitivity Analysis for the Ratio of Variances in Replicates

# Sensitivity Analysis #1:
# This is run to test the sensitivity of the fit, in a single stage, to the 
# alpha parameters (specifically the intercept parameter)
if (sensitivity1) {
  n <- 1000
  replicates <- 500
  param_estimates <- 2
  
  scenarios <- list(list(a10=-3,a11=-3),list(a10=-1,a11=-3),list(a10=0,a11=-3),list(a10=1,a11=-3),list(a10=3,a11=-3),
                    list(a10=-3,a11=-1),list(a10=-1,a11=-1),list(a10=0,a11=-1),list(a10=1,a11=-1),list(a10=3,a11=-1),
                    list(a10=-3,a11=0),list(a10=-1,a11=0),list(a10=0,a11=0),list(a10=1,a11=0),list(a10=3,a11=0),
                    list(a10=-3,a11=1),list(a10=-1,a11=1),list(a10=0,a11=1),list(a10=1,a11=1),list(a10=3,a11=1),
                    list(a10=-3,a11=3),list(a10=-1,a11=3),list(a10=0,a11=3),list(a10=1,a11=3),list(a10=3,a11=3))
  
  registerDoParallel(cores=4)
  
  l1 <- foreach(s = scenarios) %dopar% {
    results <- matrix(nrow=replicates, ncol=param_estimates)
    for (ii in 1:replicates) {
      X <- rnorm(n, 1, 1)
      W1 <- X + rnorm(n,0,0.85)
      W2 <- X + rnorm(n,0,0.85)
      Wr <- RC(cbind(W1,W2), NULL, F)
      A1 <- rbinom(n, 1, expit(s$a10+s$a11*W1))
      
      Y <- X + A1*(1 + X) + rnorm(n,0,2)
      
      wts <- abs(A1 - fitted(glm(A1~Wr, family = binomial)))
      coeffs <- lm(Y~A1*Wr, weights = wts)$coefficients
      results[ii, ] <- c(coeffs['A1'], coeffs['A1:Wr'])
    }
    results
  }
  
  save(l1, file="sensitivity1_results.RData")
  
  jj<-1
  par(mfrow=c(5,5))
  for (li in l1) {
    boxplot(li, outline=F,
            main=bquote(~alpha[0] == .(scenarios[[jj]]$a10) ~ " and " ~alpha[1]== .(scenarios[[jj]]$a11)))
    abline(h=1, lty=1)
    jj <- jj+1
  }
}

# Sensitivity Analysis #2
# This is run to test how sensitive the results are to the heteroskedastic errors
# induced by the pseudo outcomes in a 2-stage setting.
# This is accomplished by varying the amount of variance that the error-term would
# have, while fixing the overall variance of the observations.
if(sensitivity2) {
  scenarios <- list(list(sig.u = 0.25, sig.eps = 1.875),
                    list(sig.u = 0.5, sig.eps = 1.75),
                    list(sig.u = 1, sig.eps = 1.5),
                    list(sig.u = 1.25, sig.eps = 1.375),
                    list(sig.u = 1.5, sig.eps = 1.25),
                    list(sig.u = 2, sig.eps = 1),
                    list(sig.u = 2.25, sig.eps = 0.875),
                    list(sig.u = 2.5, sig.eps = 0.75),
                    list(sig.u = 2.75, sig.eps = 0.625),
                    list(sig.u = 3, sig.eps = 0.5),
                    list(sig.u = 3.5, sig.eps = 0.25),
                    list(sig.u = 4, sig.eps = 0))
  replicates <- 10000
  registerDoParallel(cores=4)
  
  l1 <- foreach(s = scenarios, .packages = c("matrixcalc","corpcor")) %dopar% {
    n <- 1000
    results <- matrix(nrow=replicates,ncol=2)
    for (ii in 1:replicates){
      X1 <- rnorm(n, 0, 1)
      W11 <- X1 + rnorm(n,0,1)
      W12 <- X1 + rnorm(n,0,1)
      Wr <- RC(cbind(W11,W12),NULL,FALSE)
      
      A1 <- rbinom(n, 1, 0.5)
      A2 <- rbinom(n, 1, 0.5)
      
      Y.tilde <- X1 + A1*(1 + X1) - A2*rnorm(n,0,sqrt(s$sig.u)) + rnorm(n,0,sqrt(s$sig.eps))
      
      wts1 <- abs(A1 - 0.5)
      
      coeffs <- lm(Y.tilde~Wr*A1, weights = wts1)$coefficients
      
      results[ii, ] <- c(coeffs['A1'],coeffs['Wr:A1'])
    }
    
    results
  }
  
  save(l1, file="sensitivity2_results.RData")
  
  jj<-1
  par(mfrow=c(3,4))
  for (li in l1) {
    colnames(li) <- c(0,1)
    boxplot(li, 
            outline=F, 
            ylim=c(0.6,1.4),
            xlab=expression(psi[1]),
            ylab="Blip Estimate",
            main=bquote(~sigma[U]^2 == .(scenarios[[jj]]$sig.u) ~ " ; " ~sigma[epsilon]^2 == .(scenarios[[jj]]$sig.eps)))
    abline(h=1, lty=1)
    jj <- jj+1
  }
  
}

# Sensitivity Analysis #3
# This is run to test how sensitive the results are to non-unit variance ratios
# as is observed in the STAR*D dataset.
if(sensitivity3) {
  
  ratios <- seq(0.5,2,length.out = 16)
  registerDoParallel(cores=4)
  
  
  l1 <- foreach(r = ratios, .packages = c("matrixcalc","corpcor")) %dopar% {
    n <- 1000
    replicates <- 10000
    results <- matrix(nrow=replicates,ncol=2)
    for (ii in 1:replicates) {
      X <- rnorm(n)
      sd_u1 <- 0.5
      sd_u2 <- sqrt(sd_u1^2*r)
      
      W1 <- X + rnorm(n,0,sd_u1)
      W2 <- X + rnorm(n,0,sd_u2)
      Wr <- RC(cbind(W1,W2),NULL,FALSE)
      
      A <- rbinom(n,1,expit(W1))
      
      Y <- 1 + X + A*(1 + X) + rnorm(n,0,4)
      
      results[ii,] <- lm(Y~Wr*A, 
                         weights = abs(A-fitted(glm(A~Wr, family = binomial))))$coefficients[c('A','Wr:A')]
      
    }
    results
  }
  
  save(l1, file="sensitivity3_results.RData")
  
  jj<-1
  par(mfrow=c(4,4))
  for (li in l1) {
    colnames(li) <- c(0,1)
    boxplot(li, 
            outline=F, 
            xlab=expression(psi[1]),
            ylab="Blip Estimate",
            main=paste("W1-to-W2 variance Ratio", ratios[jj]))
    abline(h=1, lty=1)
    jj <- jj+1
  }
  
}