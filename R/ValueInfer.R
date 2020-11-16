#' Estimation and inference for individualized treatment rule using pooled-and-split de-correlated score
#'
#' This function implements the EARL estimation for individualized treatment rule and the inference procedure based on the pooled-and-split de-correlated score (see reference).
#'
#' @param data A list - list(predictor = x, treatment = trt, outcome = y)), where x is the covariate matrix, trt is 0 or 1 (1 indicates treatment), y is the outcome.
#' @param method method can be either 'ITRFit'  or 'QLearn'.
#' @param trainingfrac determines the traning sample size by \code{floor}(total sample size * \code{trainingfrac}).
#' @param resamplingIter determines how many resamplings. By default, it is 1000.
#' @param ... Other parameters in \code{ITRFitInfer} and \code{QLearnFit}
#' @return A list
#' \describe{
#' \item{value}{Estimation of the value function given the decision rule estimated by PEARL or Q-learning}
#' \item{se}{Estimated standard error of the value function}
#' \item{upper}{\code{value+1.96*se}}
#' \item{lower}{\code{value-1.96*se}}
#' }
#' @author Muxuan Liang <mliang@fredhutch.org>
#' @references Muxuan Liang, Young-Geun Choi, Yang Ning, Maureen Smith, Yingqi Zhao (2020). Estimation and inference on high-dimensional individualized treatment rule in observational data using split-and-pooled de-correlated score.
#' @export

ValueInfer <- function(data, method = 'ITRFit', trainingfrac=0.5*log(NROW(data$predictor)), resamplingIter = 2000, ...){
  totalSampleSize <- NROW(data$predictor)
  numberPredictor <- NCOL(data$predictor)
  index_seq <- c(1:totalSampleSize)
  aipw <- NULL
  for(iter in 1:resamplingIter){
    check_balance <- FALSE
    while(!check_balance){
      training_index <- sample(index_seq, floor(totalSampleSize*1/trainingfrac))
      check_balance <- (sum(data$treatment[training_index]) > 5) & (sum(data$treatment[training_index]) < (length(training_index)-5))
    }
    testing_index <- setdiff(index_seq, training_index)
    data_train <- list(predictor = data$predictor[training_index,], treatment = data$treatment[training_index], outcome=data$outcome[training_index])
    stopifnot(method %in% c('ITRFit', 'QLearn'))
    if (method == 'ITRFit'){
      fit <- ITRFitInfer(data_train, test=FALSE,...)
      predictedTreatment <- sign(predict(fit$fit[[1]]$fit, newx = data$predictor[testing_index,], s = fit$fit[[1]]$fit$lambda.min)+predict(fit$fit[[2]]$fit, newx = data$predictor[testing_index,], s = fit$fit[[2]]$fit$lambda.min))
    } else if (method =='QLearn') {
      fit <- QLearnFit(data_train,...)
      predictedTreatment <- sign(data$predictor[testing_index,] %*% fit$fit$glmnet.fit$beta[1:numberPredictor,fit$fit$lambda==fit$fit$lambda.min])
    }
    tmp_training_index <- (index_seq %in% training_index)
    predictedOutcomeAll <- getOutcomeModel(data, method = 'kernel', sampleSplitIndex = tmp_training_index, predictAll = TRUE, screeningMethod = "SIRS")
    predictedOutcome <- NULL
    predictedOutcome$control <- predictedOutcomeAll$control[!tmp_training_index]
    predictedOutcome$treatment <- predictedOutcomeAll$treatment[!tmp_training_index]
    predictedPropensityAll <- getPropensityModel(data, method = 'kernel', sampleSplitIndex = tmp_training_index, predictAll = TRUE, screeningMethod = "SIRS")
    predictedPropensity <- predictedPropensityAll[!tmp_training_index]
    aipw_tmp <- rep(NA, totalSampleSize)
    aipw_tmp[testing_index] <- (predictedTreatment==sign(data$treatment[testing_index]-0.5))/(predictedPropensity*data$treatment[testing_index]+(1-predictedPropensity)*(1-data$treatment[testing_index]))*(data$outcome[testing_index]-data$treatment[testing_index]*predictedOutcome$treatment-(1-data$treatment[testing_index])*predictedOutcome$control)+
      (predictedTreatment>0)*predictedOutcome$treatment+(predictedTreatment<0)*predictedOutcome$control
    aipw <- cbind(aipw, aipw_tmp)
    print(iter)
  }
  if (resamplingIter > 1){
    valueMean <- mean(apply(aipw, 2, function(t){mean(t, na.rm = TRUE)}))
    valueSe <- sd(apply(aipw, 1, function(t){mean(t, na.rm = TRUE)}))/sqrt(totalSampleSize)
  } else {
    valueMean <- mean(aipw, na.rm = TRUE)
    valueSe <- sd(aipw, na.rm = TRUE)/sqrt(length(testing_index))
  }
  return(list(value = valueMean, se = valueSe, upper = valueMean+1.96*valueSe, lower = valueMean-1.96*valueSe))
}
