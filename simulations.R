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
linear_tf <- function(x) { return(x) }
quad_tf <- function(x) { return(x + x^2) }
exp_tf <- function(x) { return(x + exp(x)) }
cube_exp_tf <- function(x) { return(x + x^3 + exp(x)) }
comp_tf <- function(x) { return(as.integer(x>-0.5)*exp(x)) }

## Simulation Controls:
##    Set the following variables to 'true' or 'false' based on which settings
##    you are interested in running.
setup1 <- F # One Stage double robustness check
setup2 <- F # Two Stage double robustness (and previous treatment dependence) check
setup3 <- F # Future Treatment Predictions ()

## Simulation Code
# Taken as similar to Simulation #1 in dWOLS Paper (Wallace, Moodie 2015)
# 4 Analysis with 2 Scenarios
  # 1. Nothing Correctly Specified
  # 2. Treatemnt Correct, Treatment Free Incorrect
  # 3. Treatment Incorrect, Treatment Free Correct
  # 4. Both correctly specified
# Ran w/ Naive [based on W1] and corrected [based on Wr] covariates
if (setup1) {
  n <- 1000
  replicates <- 10000
  param_estimates <- 16
  results <- matrix(nrow=replicates, ncol=param_estimates)
  
  for (ii in 1:replicates){
    a0 <- 1
    a1 <- -0.5
    a2 <- 1.5*exp(-1)
    X <- rnorm(n, 0, 1)
    W1 <- X + rnorm(n, 0, 0.5)
    W2 <- X + rnorm(n, 0, 0.5)
    
    A <- rbinom(n, 1, expit(a0+a1*W1+a2*exp(W1)))
    
    Y <- exp(X) - X^3 + X + A*(1+X) + rnorm(n, 0, 1)
    Wr <- RC(cbind(W1, W2), NULL, FALSE)
    
    analysis1.wts <- abs(A - fitted(glm(A~Wr, family = 'binomial')))
    analysis2.wts <- abs(A - fitted(glm(A~Wr+exp(Wr), family = 'binomial')))
    analysis3.wts <- abs(A - fitted(glm(A~Wr, family = 'binomial')))
    analysis4.wts <- abs(A - fitted(glm(A~Wr+exp(Wr), family = 'binomial')))
    
    analysis5.wts <- abs(A - fitted(glm(A~W1, family = 'binomial')))
    analysis6.wts <- abs(A - fitted(glm(A~W1+exp(W1), family = 'binomial')))
    analysis7.wts <- abs(A - fitted(glm(A~W1, family = 'binomial')))
    analysis8.wts <- abs(A - fitted(glm(A~W1+exp(W1), family = 'binomial')))
    
    analysis1.coeffs <- lm(Y~Wr*A, weights = analysis1.wts)$coefficients
    analysis2.coeffs <- lm(Y~Wr*A, weights = analysis2.wts)$coefficients
    analysis3.coeffs <- lm(Y~exp(Wr)+I(Wr^3)+Wr*A, weights = analysis3.wts)$coefficients
    analysis4.coeffs <- lm(Y~exp(Wr)+I(Wr^3)+Wr*A, weights = analysis4.wts)$coefficients
    
    analysis5.coeffs <- lm(Y~W1*A, weights = analysis5.wts)$coefficients
    rm(list=c("l1"))
    analysis6.coeffs <- lm(Y~W1*A, weights = analysis6.wts)$coefficients
    analysis7.coeffs <- lm(Y~exp(W1)+I(W1^3)+W1*A, weights = analysis7.wts)$coefficients
    analysis8.coeffs <- lm(Y~exp(W1)+I(W1^3)+W1*A, weights = analysis8.wts)$coefficients
    
    results[ii,] <- c(analysis1.coeffs["A"], 
                      analysis2.coeffs["A"], 
                      analysis3.coeffs["A"], 
                      analysis4.coeffs["A"], 
                      analysis1.coeffs["Wr:A"], 
                      analysis2.coeffs[
                        "Wr:A"], 
                      analysis3.coeffs["Wr:A"], 
                      analysis4.coeffs["Wr:A"],
                      
                      analysis5.coeffs["A"], 
                      analysis6.coeffs["A"], 
                      analysis7.coeffs["A"], 
                      analysis8.coeffs["A"], 
                      analysis5.coeffs["W1:A"], 
                      analysis6.coeffs["W1:A"], 
                      analysis7.coeffs["W1:A"],  
                      analysis8.coeffs["W1:A"])
  }
  
  
  colnames(results) <- c("1", "2", "3", "4",
                         "1", "2", "3", "4",
                         "1", "2", "3", "4",
                         "1", "2", "3", "4")
  
  save(results, file="simulation1_results.RData")
  
  par(mfrow=c(2,2))
  boxplot(results[,1:4],outline=F, xlab="Analysis", ylab=expression(paste(psi[0], " estimates")))
  abline(h=1, lty=2)
  boxplot(results[,5:8],outline=F, xlab="Analysis", ylab=expression(paste(psi[1], " estimates")))
  abline(h=1, lty=2)
  mtext(paste("With Regression Calibration (n=",n,")"), outer = F, side = 3, line = 2, at = -0.3)
  
  boxplot(results[,9:12], outline=F, xlab="Analysis", ylab=expression(paste(psi[0], " estimates")))
  abline(h=1, lty=2)
  boxplot(results[,13:16],  outline=F, xlab="Analysis", ylab=expression(paste(psi[1], " estimates")))
  abline(h=1, lty=2)
  mtext(paste("Without Regression Calibration (n=",n,")"), outer = F, side = 3, line = 2, at = -0.3)
}


# Taken as similar to Simulation #2 in dWOLS Paper (Wallace, Moodie 2015)
# Two Stage Analysis, with 3 Different Setup
  # 1. Varying the alpha intercept (treatment probability)
  # 2. Varying the psi slope parameter (optimal treatment threshold)
  # 3. Varying the TF model (double robustness in 2 stage)
# Ran w/ Treatment Correctly Specified; Treatment Free is always Linear
  # 5 setups for each scenario above
if(setup2) {
  n <- 1000
  replicates <- 10000
  param_estimates <- 4
    
  scenarios <- list(list(a10=-2, a20=-2, a11=1, a21=1, p10=1, p11=1, p20=1, p21=1, tf=linear_tf),
                     list(a10=-1, a20=-1, a11=1, a21=1, p10=1, p11=1, p20=1, p21=1, tf=linear_tf),
                     list(a10=0, a20=0, a11=1, a21=1, p10=1, p11=1, p20=1, p21=1, tf=linear_tf),
                     list(a10=1, a20=1, a11=1, a21=1, p10=1, p11=1, p20=1, p21=1, tf=linear_tf),
                     list(a10=2, a20=2, a11=1, a21=1, p10=1, p11=1, p20=1, p21=1, tf=linear_tf),
                    
                     list(a10=0, a20=0, a11=1, a21=1, p10=1, p11=-1, p20=1, p21=-1, tf=linear_tf),
                     list(a10=0, a20=0, a11=1, a21=1, p10=1, p11=-.1, p20=1, p21=-.1, tf=linear_tf),
                     list(a10=0, a20=0, a11=1, a21=1, p10=1, p11=0, p20=1, p21=0, tf=linear_tf),
                     list(a10=0, a20=0, a11=1, a21=1, p10=1, p11=.1, p20=1, p21=.1, tf=linear_tf),
                     list(a10=0, a20=0, a11=1, a21=1, p10=1, p11=1, p20=1, p21=1, tf=linear_tf),
                    
                      list(a10=0, a20=0, a11=1, a21=1, p10=1, p11=1, p20=1, p21=1, tf=linear_tf),
                      list(a10=0, a20=0, a11=1, a21=1, p10=1, p11=1, p20=1, p21=1, tf=quad_tf),
                      list(a10=0, a20=0, a11=1, a21=1, p10=1, p11=1, p20=1, p21=1, tf=exp_tf),
                      list(a10=0, a20=0, a11=1, a21=1, p10=1, p11=1, p20=1, p21=1, tf=cube_exp_tf),
                      list(a10=0, a20=0, a11=1, a21=1, p10=1, p11=1, p20=1, p21=1, tf=comp_tf))
  
  registerDoParallel(cores=4)
  
  l1 <- foreach(s = scenarios, .packages = c("matrixcalc","corpcor")) %dopar% {
    results <- matrix(nrow=replicates, ncol=param_estimates)
    for (ii in 1:replicates) {
      # Generate the Data
      X1 <- rnorm(n, 0,1)
      W11 <- X1 + rnorm(n, 0, 0.5)
      W12 <- X1 + rnorm(n, 0, 0.5)
      Wr1 <- RC(cbind(W11, W12), NULL, FALSE)
      A1 <- rbinom(n, 1, expit(s$a10 + s$a11*W11))
      
      X2 <- rnorm(n,A1,1)
      W21 <- X2 + rnorm(n, 0, 1)
      W22 <- X2 + rnorm(n, 0, 1)
      Wr2 <- RC(cbind(W21, W22), as.matrix(A1), FALSE)
      A2 <- rbinom(n, 1, expit(s$a20 + s$a21*W21))
      
      Y <- s$tf(X1) + A1*(s$p10 + s$p11*X1) + A2*(s$p20 + s$p21*X2) + rnorm(n, 0, 2)
      
      # Fit the model
      wts2 <- abs(A2-fitted(glm(A2~Wr2, family=binomial)))
      
      # Stage 2:
      stage2c <- lm(Y~Wr1*A1 + A2 + Wr2:A2, weights = wts2)$coefficients
      
      # Stage 1:
      Ytilde <- Y - A2*(stage2c['A2']+stage2c['A2:Wr2']*W22)
      wts1 <- abs(A1-fitted(glm(A1~Wr1, family=binomial)))
      stage1c <- lm(Ytilde~Wr1*A1, weights = wts1)$coefficients
      
      results[ii,] <- c(stage1c['A1'], stage1c['Wr1:A1'], 
                        stage2c['A2'], stage2c['A2:Wr2'])
    }
    results
  }
  
  save(l1, file="simulation2_results.RData")
  
  jj <- 1
  par(mfrow=c(3,5))
  
  for (li in l1) {
    boxplot(li, outline=F)
    abline(h=scenarios[[jj]]$p10, lty=2)
    abline(h=scenarios[[jj]]$p11, lty=2)
    jj <- jj+1
  }
  
}

# Simulation #3:
# This is a run to predict the predictive strength of our method in the 
# 2-stage setting, as compared to treating this as a naive prediction problem.
# 5 Different Analysis Scenarios are Compared
  # 1. Naive Fitting, Naive Validation
  # 2. Corrected Fitting, Pseudo-Corrected Validation
  # 3. Corrected Fitting, Corrected Validation
  # 4. Naive Fitting, True Validation
  # 5. Corrected Fitting, True Validation
# We use linear blip models with:
  # X1 ~ N(0,1); W ~ X1 + N(0,1)
  # X2 ~ N(A1, 1); W ~ X2 + N(0,0.25)
  # A1 = expit(1 - W1)
  # A2 = expit(1 - W2)
  # Y = X1 + A1(1 - X1) + A2(3-2X2) + N(0,4)
# n = 1000 for fitting; 5000 for validation

if (setup3) {
  n_fit <- 1000
  n_val <- 5000
  replicates <- 1000
  n_param <- 35
  results <- matrix(nrow=replicates,ncol=n_param)
  
  for (ii in 1:replicates){
    X1 <- rnorm(n_fit, 0, 1)
    W11 <- X1 + rnorm(n_fit,0,1)
    W12 <- X1 + rnorm(n_fit,0,1)
    Wr1.obj <- RC(cbind(W11,W12),NULL,TRUE)
    A1 <- rbinom(n_fit,1,expit(1-W11))
    
    X2 <- rnorm(n_fit,A1,1)
    W21 <- X2 + rnorm(n_fit,0,0.5)
    W22 <- X2 + rnorm(n_fit,0,0.5)
    Wr2.obj <- RC(cbind(W21,W22), as.matrix(A1), TRUE)
    A2 <- rbinom(n_fit,1,expit(1-W21))
    
    Y <- X1 + A1*(1-X1) + A2*(3-2*X2) + rnorm(n_fit,0,2)
    
    Wr1 <- Wr1.obj$X.imp
    Wr2 <- Wr2.obj$X.imp
    
    # Get fitted coefficients
    # Corrected
    wts2.corr <- abs(A2-fitted(glm(A2~Wr2, family=binomial)))
    coeffs2.corr <- lm(Y~Wr1*A1+A2+A2:Wr2, weights=wts2.corr)$coefficients
    
    Ytilde.corr <- Y - A2*(coeffs2.corr['A2'] + coeffs2.corr['A2:Wr2']*W22)
    wts1.corr <- abs(A1-fitted(glm(A1~Wr1, family=binomial)))
    coeffs1.corr <- lm(Ytilde.corr~Wr1*A1, weights=wts1.corr)$coefficients
    
    # Naive
    wts1.naive <- abs(A2-fitted(glm(A2~W21, family=binomial)))
    coeffs2.naive <- lm(Y~W11*A1+A2+A2:W21, weights=wts2.corr)$coefficients
    
    Ytilde.naive <- Y - A2*(coeffs2.naive['A2'] + coeffs2.naive['A2:W21']*W21)
    wts1.naive <- abs(A1-fitted(glm(A1~W11, family=binomial)))
    coeffs1.naive <- lm(Ytilde.naive~W11*A1, weights=wts1.naive)$coefficients
    
    ## Generate Validation Data
    X1.v <- rnorm(n_val, 0, 1)
    W11.v <- X1.v + rnorm(n_val, 0, 1)
    W12.v <- X1.v + rnorm(n_val, 0, 1)
    Wr1.v <- RC(cbind(W11.v, W12.v), NULL, FALSE)
    Wp1.v <- Wr1.obj$Mu.W + Wr1.obj$Sigma.XX*(W11.v - Wr1.obj$Mu.W)/(Wr1.obj$Sigma.XX + Wr1.obj$Sigma.UU)
    
    
    A1.opt <- as.integer(1-X1.v > 0)
    
    A1.nn <- as.integer(coeffs1.naive['A1']+coeffs1.naive['W11:A1']*W11.v > 0)
    A1.nt <- as.integer(coeffs1.naive['A1']+coeffs1.naive['W11:A1']*X1.v > 0)
    A1.cc <- as.integer(coeffs1.corr['A1']+coeffs1.corr['Wr1:A1']*Wr1.v > 0)
    A1.ct <- as.integer(coeffs1.corr['A1']+coeffs1.corr['Wr1:A1']*X1.v > 0)
    A1.cp <- as.integer(coeffs1.corr['A1']+coeffs1.corr['Wr1:A1']*Wp1.v > 0)
    
    ### [First Optimal was Taken]
    X2.opt.v <- rnorm(n_val, A1.opt, 1)
    W21.opt.v <- X2.opt.v + rnorm(n_val, 0, 1)
    W22.opt.v <- X2.opt.v + rnorm(n_val, 0, 1)
    Wr2.opt.v <- RC(cbind(W21.opt.v, W22.opt.v), as.matrix(A1.opt), FALSE)
    Wp2.opt.v <- Wr2.obj$Mu.W + cbind(as.matrix(Wr2.obj$Sigma.XX), Wr2.obj$Sigma.XZ)%*%
                              solve(rbind(cbind(as.matrix(Wr2.obj$Sigma.XX + Wr2.obj$Sigma.UU),Wr2.obj$Sigma.XZ),
                                          cbind(t(Wr2.obj$Sigma.XZ), Wr2.obj$Sigma.ZZ)))%*%
                              t(cbind(as.matrix(W21.opt.v - Wr2.obj$Mu.W), as.matrix(A1.opt - Wr2.obj$Z.bar)))
    
    A2.opt.opt <- as.integer(3-2*X2.opt.v > 0)
    A2.opt.nn <- as.integer(coeffs2.naive['A2']+coeffs2.naive['A2:W21']*W21.opt.v > 0)
    A2.opt.nt <- as.integer(coeffs2.naive['A2']+coeffs2.naive['A2:W21']*X2.opt.v > 0)
    A2.opt.cc <- as.integer(coeffs2.corr['A2']+coeffs2.corr['A2:Wr2']*Wr2.opt.v > 0)
    A2.opt.ct <- as.integer(coeffs2.corr['A2']+coeffs2.corr['A2:Wr2']*X2.opt.v > 0)
    A2.opt.cp <- as.integer(coeffs2.corr['A2']+coeffs2.corr['A2:Wr2']*Wp2.opt.v > 0)
    
    ### [First Naive was Taken]
    ### [Naive, Naive]
    X2.nn.v <- rnorm(n_val, A1.nn, 1)
    W21.nn.v <- X2.nn.v + rnorm(n_val, 0, 1)
    
    A2.nn.opt <- as.integer(3-2*X2.nn.v > 0)
    A2.nn <- as.integer(coeffs2.naive['A2']+coeffs2.naive['A2:W21']*W21.nn.v > 0)
    
    ### [Naive, Truth]
    X2.nt.v <- rnorm(n_val, A1.nt, 1)
    
    A2.nt.opt <- as.integer(3-2*X2.nt.v > 0)
    A2.nt <- as.integer(coeffs2.naive['A2']+coeffs2.naive['A2:W21']*X2.nt.v > 0)
    
    
    ## [Corrected, Corrected]
    X2.cc.v <- rnorm(n_val, A1.cc, 1)
    W21.cc.v <- X2.cc.v + rnorm(n_val, 0, 1)
    W22.cc.v <- X2.cc.v + rnorm(n_val, 0, 1)
    Wr2.cc.v <- RC(cbind(W21.cc.v, W22.cc.v), as.matrix(A1.cc), FALSE)
    
    A2.cc.opt <- as.integer(3-2*X2.cc.v > 0)
    A2.cc <- as.integer(coeffs2.corr['A2']+coeffs2.corr['A2:Wr2']*Wr2.cc.v > 0)
    
    ## [Corrected, Truth]
    X2.ct.v <- rnorm(n_val, A1.ct, 1)
    A2.ct.opt <- as.integer(3-2*X2.ct.v > 0)
    A2.ct <- as.integer(coeffs2.corr['A2']+coeffs2.corr['A2:Wr2']*X2.ct.v > 0)
    
    ## [Corrected, Pseudo-corrected]
    X2.cp.v <- rnorm(n_val, A1.cp, 1)
    W21.cp.v <- X2.cp.v + rnorm(n_val, 0, 1)
    
    W2p.cp.v <- Wr2.obj$Mu.W + cbind(as.matrix(Wr2.obj$Sigma.XX), Wr2.obj$Sigma.XZ)%*%
                  solve(rbind(cbind(as.matrix(Wr2.obj$Sigma.XX + Wr2.obj$Sigma.UU), 
                                    Wr2.obj$Sigma.XZ),
                              cbind(t(Wr2.obj$Sigma.XZ), Wr2.obj$Sigma.ZZ)))%*%
                              t(cbind(as.matrix(W21.cp.v - Wr2.obj$Mu.W), as.matrix(A1.cp - Wr2.obj$Z.bar)))
    
    A2.cp.opt <- as.integer(3-2*X2.cp.v > 0)
    A2.cp <- as.integer(coeffs2.corr['A2']+coeffs2.corr['A2:Wr2']*W2p.cp.v > 0)
    
    results[ii,] <- c(sum(A1.opt == A1.nn)/n_val,
                      sum(A1.opt == A1.nt)/n_val,
                      sum(A1.opt == A1.cc)/n_val,
                      sum(A1.opt == A1.ct)/n_val,
                      sum(A1.opt == A1.cp)/n_val,
                      # Stage 2 (Optimal)
                        sum(A2.opt.opt == A2.opt.nn)/n_val,
                        sum(A2.opt.opt == A2.opt.nt)/n_val,
                        sum(A2.opt.opt == A2.opt.cc)/n_val,
                        sum(A2.opt.opt == A2.opt.ct)/n_val,
                        sum(A2.opt.opt == A2.opt.cp)/n_val,
                      # Stage 2 (Given)
                        sum(A2.nn.opt == A2.nn)/n_val,
                        sum(A2.nt.opt == A2.nt)/n_val,
                        sum(A2.cc.opt == A2.cc)/n_val,
                        sum(A2.ct.opt == A2.ct)/n_val,
                        sum(A2.cp.opt == A2.cp)/n_val,
                      # Both Correct
                        sum(as.integer(A1.opt==A1.nn)*as.integer(A2.nn.opt==A2.nn))/n_val,
                        sum((A1.opt==A1.nt)*(A2.nt.opt==A2.nt))/n_val,
                        sum((A1.opt==A1.cc)*(A2.cc.opt==A2.cc))/n_val,
                        sum((A1.opt==A1.ct)*(A2.ct.opt==A2.ct))/n_val,
                        sum((A1.opt==A1.cp)*(A2.cp.opt==A2.cp))/n_val,
                      # Neither Correct
                        sum((A1.opt!=A1.nn)*(A2.nn.opt!=A2.nn))/n_val,
                        sum((A1.opt!=A1.nt)*(A2.nt.opt!=A2.nt))/n_val,
                        sum((A1.opt!=A1.cc)*(A2.cc.opt!=A2.cc))/n_val,
                        sum((A1.opt!=A1.ct)*(A2.ct.opt!=A2.ct))/n_val,
                        sum((A1.opt!=A1.cp)*(A2.cp.opt!=A2.cp))/n_val,
                      # One not Two
                        sum((A1.opt==A1.nn)*(A2.nn.opt!=A2.nn))/n_val,
                        sum((A1.opt==A1.nt)*(A2.nt.opt!=A2.nt))/n_val,
                        sum((A1.opt==A1.cc)*(A2.cc.opt!=A2.cc))/n_val,
                        sum((A1.opt==A1.ct)*(A2.ct.opt!=A2.ct))/n_val,
                        sum((A1.opt==A1.cp)*(A2.cp.opt!=A2.cp))/n_val,
                      # Two not One
                        sum((A1.opt!=A1.nn)*(A2.nn.opt==A2.nn))/n_val,
                        sum((A1.opt!=A1.nt)*(A2.nt.opt==A2.nt))/n_val,
                        sum((A1.opt!=A1.cc)*(A2.cc.opt==A2.cc))/n_val,
                        sum((A1.opt!=A1.ct)*(A2.ct.opt==A2.ct))/n_val,
                        sum((A1.opt!=A1.cp)*(A2.cp.opt==A2.cp))/n_val)
  }
  
  colnames(results) <- rep(c("NN", "NT", "CC", "CT", "CP"), 7)
  
  save(results, file="simulation3_results.RData")
  
  layout(matrix(c(1,1,1,1,2,2,3,3,4,5,6,7), 3, 4, byrow=T))
  boxplot(results[,1:5], main="Stage 1 Correct", outline=F)
  boxplot(results[,6:10], main="Stage 2 Correct (Optimal 1)", outline=F)
  boxplot(results[,11:15], main="Stage 2 Correct (Treated 1)", outline=F)
  boxplot(results[,16:20], main="Both 1 and 2 Correct", outline=F)
  boxplot(results[,21:25], main="Neither 1 nor 2 Correct", outline=F)
  boxplot(results[,26:30], main="1 Correct, 2 Incorrect", outline=F)
  boxplot(results[,31:35], main="2 Correct, 1 Incorrect", outline=F)
}