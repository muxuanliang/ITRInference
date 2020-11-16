#' This function implements the Q-learning estimation for individualized treatment rule and the inference procedure based on the de-correlated score (see reference).
#'
#' @param data A list - list(predictor = x, treatment = trt, outcome = y)), where x is the covariate matrix, trt is 0 or 1 (1 indicates treatment), y is the outcome.
#' @param intercept includes intercept or not
#' @param standardize whether standardize the input covariant matrix.
#' @return A list
#'  \describe{
#' \item{fit}{Use \code{glmnet} to fit a nigh-dimensional Q-learning.}
#' \item{pseudoPredictor}{Augmented design matrix with covariate and treatment interations for Q-learning.}
#' \item{pseudoTreatment}{Treatment coded in 1 and -1.}
#' \item{pseudoOutcome}{Outcome used in the inferenc step.}
#' }
#' @author Muxuan Liang <mliang@fredhutch.org>
#' @references Muxuan Liang, Young-Geun Choi, Yang Ning, Maureen Smith, Yingqi Zhao (2020). Estimation and inference on high-dimensional individualized treatment rule in observational data using split-and-pooled de-correlated score.
#' @export
QLearnFit <- function(data, intercept=FALSE, standardize = TRUE){
  size <- dim(data$predictor)[1]
  Outcome <- data$outcome
  Treatment <- (data$treatment - 0.5) * 2
  pseudoPredictor <- cbind(apply(data$predictor,2,function(t){t*Treatment}), Treatment, data$predictor)
  if(!intercept){
    pseudoPredictor <- cbind(apply(data$predictor,2,function(t){t*Treatment}), data$predictor)
  }
  fit <- glmnet::cv.glmnet(x=pseudoPredictor, y=Outcome, family='gaussian', intercept = TRUE, standardize = standardize)
  list(fit=fit, pseudoPredictor = pseudoPredictor, pseudoTreatment = Treatment, pseudoOutcome = Outcome)
}
