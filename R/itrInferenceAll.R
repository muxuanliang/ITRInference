#' Estimation and inference for individualized treatment rule using pooled-and-split de-correlated score
#'
#' This function implements the EARL estimation for individualized treatment rule and the inference procedure based on the pooled-and-split de-correlated score (see reference).
#'
#' @param data A list - list(predictor = x, treatment = trt, outcome = y)), where x is the covariate matrix, trt is 0 or 1 (1 indicates treatment), y is the outcome.
#' @param propensity estimated propensity score p(trt=1|X). This should be used only when the propensity is estimated by a parametric model.
#' @param outcome stimated outcome model E(y|trt, X). This should be used only when the outcome is estimated by a parametric model. It should be a list (list(control=..., treatment=...)).
#' @param outcomeModel this selects method to estimate the outcome when outcome=NULL. Options include lm, glmnet, kernel, and gam. If lm is used, user also need to input the outcomeFormula like y~x used in lm. By default, kernel regression is selected.
#' @param propensityModel Similar to outcomeModel.
#' @param intercept includes intercept or not
#' @param test whether return p-value and sd estimation.
#' @param indexToTest indicates which coefficients to test. By default, c(1:8)
#' @param parallel whether use parallel computing; by default, FALSE.
#' @param screeningMethod when dimension is high, this selects method to preform variable screening for outcomeModel and propensityModel fitting. Methods in VariableScreening package are available.
#' @param standardize whether standardize the input covariant matrix.
#' @return fit is estimated coefficient for ITR on two split datasets. p-values are the p-value for each coefficients included in indexToTest. (betaAN-1.96*sigmaAN/sqrt(sample size), betaAN+1.96*sigmaAN/sqrt(sample size)) provides the 95\% CI for the coefficients.
#' @author Muxuan Liang <mliang@fredhutch.org>
#' @references To be added
#' @examples
#' # generate data
#' n <- 500
#' p <- 800
#' x <- array(rnorm(n*p), c(n,p))
#' beta_pi <- 0.4*c(1,1,1,0,-1,0,-1,rep(0, times=p-7))
#' trt <- apply(x, 1, function(t){ rbinom(1, 1, prob = exp(t[1]^2+t[2]^2+t[1]*t[2])/(1+exp(t[1]^2+t[2]^2+t[1]*t[2])))})
#' beta_main <- 0.8*c(-1,-1,1,-1,0,0,rep(0, times=p-6))
#' main_effect <- exp(x %*% beta_main)
#' beta_inter <- 0.7*c(1,1,-1,-1,rep(0, times=p-4))
#' interact_effect <- (pnorm(x %*% beta_inter)-0.5)*(abs(apply(x[,c(4:8)],1,sum)) + 4 * xi)
#' e <- rnorm(n)
#' y <- main_effect + sign(trt - 0.5) * interact_effect + e
#' # set up training data
#' dataTrain <- list(predictor = x, treatment = (trt > 0), outcome = y)
#' # fit our approach
#' fit_itrAll <- ITRFitAll(data = dataTrain, outcomeModel = 'kernel', propensityModel = 'kernel')
#' # compare with q-learning and implement standard de-correlated score for q0learning
#' fit_qlearn <- QLearnFit(data = dataTrain, intercept = TRUE)
#' score_qlearn <- scoreTestQLearn(fit_qlearn, parallel=FALSE)
#' @export

ITRFitAll <- function(data, propensity = NULL, outcome = NULL, loss = 'logistic', sampleSplitIndex=NULL,
                      type.measure = 'lossFun',
                      outcomeModel='kernel', outcomeFormula = NULL,
                      propensityModel='kernel', propensityFormula = NULL,
                      intercept=FALSE, test=TRUE, indexToTest=c(1:8), parallel=FALSE,
                      screeningMethod="SIRS", outcomeScreeningFamily = 'Gaussian', standardize = TRUE){
  size <- dim(data$predictor)[1]
  if(is.null(sampleSplitIndex)){
    sampleSplitIndex <- (rnorm(size) > 0)
  }
  fit <- NULL
  fit[[1]] <- ITRFit(data = data, propensity = propensity, outcome = outcome, loss = loss, sampleSplitIndex = sampleSplitIndex,
                     type.measure = type.measure,
                     outcomeModel = outcomeModel, outcomeFormula = outcomeFormula, propensityModel = propensityModel,
                     propensityFormula = propensityFormula, intercept = intercept, screeningMethod = screeningMethod,
                     outcomeScreeningFamily = outcomeScreeningFamily, standardize = standardize)
  fit[[2]] <- ITRFit(data = data, propensity = propensity, outcome = outcome, loss = loss, sampleSplitIndex = (!sampleSplitIndex),
                     type.measure = type.measure,
                     outcomeModel = outcomeModel, outcomeFormula = outcomeFormula, propensityModel = propensityModel,
                     propensityFormula = propensityFormula, intercept = intercept, screeningMethod = screeningMethod,
                     outcomeScreeningFamily = outcomeScreeningFamily, standardize = standardize)
  if (test){
    score_1 <- scoreTest(fit[[1]], indexToTest = indexToTest, parallel=parallel, intercept = TRUE)
    score_2 <- scoreTest(fit[[2]], indexToTest = indexToTest, parallel=parallel, intercept = TRUE)
    score <- (score_1$score+ score_2$score)/2
    betaAN <- (score_1$betaAN + score_2$betaAN)/2
    sigma <- sqrt((score_1$sigma^2 + score_2$sigma^2)/2)
    I <- (score_1$I+score_2$I)/2
    sigmaAN <- sigma/I
    res <- list(fit =fit, pvalue=pnorm(-abs(sqrt(size)*score/sigma))*2, betaAN=betaAN, sigmaAN=sigmaAN)
  } else {
    res <- list(fit = fit)
  }
  res
}






