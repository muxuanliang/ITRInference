#' This function implements the inference for value of individualized treatment rule.
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

ValueInfer <- function(data, method = 'ITRFit', trainingfrac=0.5*log(NROW(data$predictor)), resamplingIter = 2000, inferMethod=c('direct', 'bootstrap','weightedbootstrap'), propensity = NULL,
                       outcome = NULL, ...){
  totalSampleSize <- NROW(data$predictor)
  numberPredictor <- NCOL(data$predictor)
  index_seq <- c(1:totalSampleSize)
  aipw <- NULL
  
  check_balance <- FALSE
  while(!check_balance){
    training_index <- sample(index_seq, floor(totalSampleSize/trainingfrac))
    check_balance <- (sum(data$treatment[training_index]) > 5) & (sum(data$treatment[training_index]) < (length(training_index)-5))
  }
  testing_index <- setdiff(index_seq, training_index)
  data_train <- list(predictor = data$predictor[training_index,], treatment = data$treatment[training_index], outcome=data$outcome[training_index])
  stopifnot(method %in% c('ITRFit', 'QLearn'))
  if (method == 'ITRFit'){
    propensity.train <- NULL
    outcome.train <- NULL
    if(!is.null(propensity)){
      propensity.train <- propensity[training_index]
    }
    if(!is.null(outcome)){
      outcome.train <- list(treatment=outcome$treatment[training_index], control=outcome$control[training_index])
    }
    fit <- ITRFitInfer(data_train, test=FALSE,propensity = propensity.train,
                       outcome = outcome.train,...)
    predictedTreatment <- sign(predict(fit$fit[[1]]$fit, newx = data$predictor[testing_index,], s = fit$fit[[1]]$fit$lambda.min)+predict(fit$fit[[2]]$fit, newx = data$predictor[testing_index,], s = fit$fit[[2]]$fit$lambda.min))
  } else if (method =='QLearn') {
    fit <- QLearnFit(data_train,...)
    predictedTreatment <- sign(data$predictor[testing_index,] %*% fit$fit$glmnet.fit$beta[1:numberPredictor,fit$fit$lambda==fit$fit$lambda.min]+fit$fit$glmnet.fit$a0[fit$fit$lambda==fit$fit$lambda.min])
  }
  tmp_training_index <- (index_seq %in% training_index)
  
  if(is.null(outcome)){
    predictedOutcomeAll <- getOutcomeModel(data, method = 'kernel', sampleSplitIndex = (!tmp_training_index), predictAll = TRUE, screeningMethod = "SIRS")
    predictedOutcome <- NULL
    predictedOutcome$control <- predictedOutcomeAll$control[!tmp_training_index]
    predictedOutcome$treatment <- predictedOutcomeAll$treatment[!tmp_training_index]
  } else {
    predictedOutcome <- NULL
    predictedOutcome$control <- outcome$control[!tmp_training_index]
    predictedOutcome$treatment <- outcome$treatment[!tmp_training_index]
  }
  if(is.null(propensity)){
    predictedPropensityAll <- getPropensityModel(data, method = 'kernel', sampleSplitIndex = (!tmp_training_index), predictAll = TRUE, screeningMethod = "SIRS")
    predictedPropensity <- predictedPropensityAll[!tmp_training_index]
  } else {
    predictedPropensity <- propensity[!tmp_training_index]
  }
  
  predictedTreatment[predictedTreatment==0] <- 1
  aipw <- (predictedTreatment==sign(data$treatment[testing_index]-0.5))/(predictedPropensity*data$treatment[testing_index]+(1-predictedPropensity)*(1-data$treatment[testing_index]))*(data$outcome[testing_index]-data$treatment[testing_index]*predictedOutcome$treatment-(1-data$treatment[testing_index])*predictedOutcome$control)+
    (predictedTreatment>0)*predictedOutcome$treatment+(predictedTreatment<0)*predictedOutcome$control
  value <- mean(aipw, na.rm = TRUE)
  valueSe <- sd(aipw, na.rm = TRUE)/sqrt(length(testing_index))
  
  if(inferMethod=='direct'){
    upper <- value+1.96*valueSe
    lower <- value-1.96*valueSe
  } else if(inferMethod=='bootstrap'){
    aipw_boot <- rep(NA, resamplingIter)
    for(iter in 1:resamplingIter){
      sampling_index <- sample(1:length(testing_index),length(testing_index), replace=TRUE)
      predictedTreatmentBoot <- predictedTreatment[sampling_index]
      observedTreatmenBoot <- data$treatment[testing_index][sampling_index]
      predictedPropensityBoot <- predictedPropensity[sampling_index]
      observedOutcomeBoot <- data$outcome[testing_index][sampling_index]
      predictedOutcomeBoot <- NULL
      predictedOutcomeBoot$treatment <- predictedOutcome$treatment[sampling_index]
      predictedOutcomeBoot$control <- predictedOutcome$control[sampling_index]
      aipw_boot[iter] <- mean((predictedTreatmentBoot==sign(observedTreatmenBoot-0.5))/(predictedPropensityBoot*observedTreatmenBoot+(1-predictedPropensityBoot)*(1-observedTreatmenBoot))*(observedOutcomeBoot-observedTreatmenBoot*predictedOutcomeBoot$treatment-(1-observedTreatmenBoot)*predictedOutcomeBoot$control)+
        (predictedTreatmentBoot>0)*predictedOutcomeBoot$treatment+(predictedTreatmentBoot<0)*predictedOutcomeBoot$control, na.rm=TRUE)
    }
    upper <- quantile(aipw_boot, probs=0.975)
    lower <- quantile(aipw_boot, probs=0.025)
  } else if(inferMethod=='weightedbootstrap'){
    aipw_boot <- rep(NA, resamplingIter)
    for(iter in 1:resamplingIter){
      weightsBoot <- rgamma(length(testing_index),1,1)
      aipw_boot[iter] <- mean(weightsBoot((predictedTreatment==sign(data$treatment[testing_index]-0.5))/(predictedPropensity*data$treatment[testing_index]+(1-predictedPropensity)*(1-data$treatment[testing_index]))*(data$outcome[testing_index]-data$treatment[testing_index]*predictedOutcome$treatment-(1-data$treatment[testing_index])*predictedOutcome$control)+
        (predictedTreatment>0)*predictedOutcome$treatment+(predictedTreatment<0)*predictedOutcome$control))
    }
    upper <- quantile(aipw_boot, probs=0.975)
    lower <- quantile(aipw_boot, probs=0.025)
  }
  return(list(value = value, se = valueSe, upper = upper, lower = lower))
}
