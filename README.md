# R Package: ITRInference
# Estimation and inference for individualized treatment rule using pooled-and-split de-correlated score

Version: 1.0.0

Author: Muxuan Liang <mliang@fredhutch.org>

Maintainer: Muxuan Liang <mliang@fredhutch.org>

Description: The package implements a doubly robust estimation for individualized treatment rule (ITR). User can choose bespoke methods to estimate outcome model and propensity. The package also provides a variaty of methods to estimate these nuisance parameters. It also returns the p-values and CI for the ITR using a pooled-and-split de-correlated score.

License: GPL-3

Imports: 
         glmnet,
         VariableScreening,
         mgcv

Encoding: UTF-8

# Details
This package (function ITRFitAll) implements the EARL estimation for individualized treatment rule and the inference procedure based on the pooled-and-split de-correlated score (see reference).

# Example

### Generate data

n <- 500

p <- 800

x <- array(rnorm(n*p), c(n,p))

beta_pi <- 0.4*c(1,1,1,0,-1,0,-1,rep(0, times=p-7))

trt <- apply(x, 1, function(t){ rbinom(1, 1, prob = exp(t[1]^2+t[2]^2+t[1]*t[2])/(1+exp(t[1]^2+t[2]^2+t[1]*t[2])))})

beta_main <- 0.8*c(-1,-1,1,-1,0,0,rep(0, times=p-6))

main_effect <- exp(x \%*\% beta_main)

beta_inter <- 0.7*c(1,1,-1,-1,rep(0, times=p-4))

interact_effect <- (pnorm(x \%*\% beta_inter)-0.5)*(abs(apply(x[,c(4:8)],1,sum)) + 4 * xi)

e <- rnorm(n)

y <- main_effect + sign(trt - 0.5) * interact_effect + e

### Set up training data

dataTrain <- list(predictor = x, treatment = (trt > 0), outcome = y)

### Fit our approach

infer_itr <- ITRFitInfe(data = dataTrain, outcomeModel = 'kernel', propensityModel = 'kernel')

### Compare with q-learning and implement standard de-correlated score for q-learning

fit_qlearn <- QLearnFit(data = dataTrain, intercept = TRUE)

infer_qlearn <- QFitInfer(fit_qlearn, parallel=FALSE)

# References
Muxuan Liang, Young-Geun Choi, Yang Ning, Maureen Smith, Yingqi Zhao (2020). Estimation and inference on high-dimensional individualized treatment rule in observational data using split-and-pooled de-correlated score. [arXiv](https://arxiv.org/abs/2007.04445)
