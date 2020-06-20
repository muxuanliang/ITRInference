#' Estimation and inference for itr using pooled-and-split de-correlated score
#'
#' This function
#'
#' @param data A list - list(predictor = x, treatment = trt, outcome = y)), where x is the covariate matrix, trt is 0 or 1 (1 indicates treatment), y is the outcome.
#' @return A matrix of the infile
#' @export

# ITRFitAll returns the pooled data-split de-correlated score. propensity is defined as p(T=1|X)
ITRFitAll <- function(data, propensity = NULL, outcome = NULL, loss = 'logistic', sampleSplitIndex=NULL,
                      type.measure = 'lossFun',
                      outcomeModel=c('lm', 'glmnet', 'kernel', 'others'), outcomeFormula = NULL,
                      propensityModel=c('lm', 'glmnet', 'kernel'), propensityFormula = NULL,
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
    res <- list(fit =fit, score = score, sigma=sigma, pvalue=pnorm(-abs(sqrt(size)*score/sigma))*2, betaAN=betaAN, sigmaAN=sigmaAN)
  } else {
    res <- list(fit = fit)
  }
  res
}






