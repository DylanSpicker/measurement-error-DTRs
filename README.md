**Author:** Dylan Spicker (dylan.spicker@uwaterloo.ca)

**Article:** Measurement error and precision medicine: error-prone tailoring covariates in dynamic treatment regimes.

The enclosed files demonstrate the implementation of the Dynamic WOLS (dWOLS) method for estimating optimal dynamic treatment regimens, in the presence of non-differential measurement error. In particular, simulation results and sensitivity analyses are included showcasing the methods across a wide variety of single- and multi-stage DTRs.

- **simulations.r** contains a sample which demonstrates the estimation procedure in the single (setup1) and multistage (setup2) scenarios, as outlined in the corresponding paper. Additionally, a multistage prediction setup, focused on the prescription of future treatments (setup3) is included. Each scenario can be toggled to run or not using the T/F controls in the preamble of the file.
- **sensitivity.r** contains a sensitivity analyses investigating the impact of effect size (scenario1), heteroskedastic pseudo-outcome error variances (scenario2), and non-unit replicate variance ratios (scenario3). Each scenario can be toggled to run or not using the T/F controls in the preamble of the file.

**rCalibration.r** contains an implementation of of the best linear approximation for regression calibration, a function (called RC) which takes:

- W (Required): The observed data matrix for the error-prone covariate, where columns are replicate observations and rows are individual observations
- Z (Optional): The observed, error-free covariate matrix, which can contain any number of error-free covariates as columns
- return_var (optional): A boolean parameter (default false), which specifies whether the estimated mean and variance structures should be returned or not

RC returns:

- $X.imp:	The estimated imputed X values (always returned. If return_var=false then returned as a vector, otherwise as a part of the object)
- $Sigma.XX, $Sigma.UU, $Sigma.ZZ, $Sigma.XZ: (Returned only if return_var=True). The relevant variance/covariance estimates based on the RC procedure.
- $Mu.W, $Z.bar: (Returned only if return_var=True). The relevant mean estimates based on the RC procedure.
- The variance/mean structures including Z are only included if Z is non null.
