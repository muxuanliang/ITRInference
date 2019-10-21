# # ITRFit obtained the ITR. propensity is defined as p(T=1|X)
# ITRFitAll <- function(data, propensity, loss = c('logistic'), outcomeModel=c('lm', 'glmnet', 'kernel', 'others'), outcomeFormula = NULL, intercept=FALSE){
#   size <- dim(data$predictor)[1]
#   sampleSplitIndex <- (rnorm(size) > 0)
#   data$predictor <- scale(data$predictor)
#   predictedOutcome <- getOutcomeModel(data, method = outcomeModel, sampleSplitIndex = sampleSplitIndex, outcomeFormula = outcomeFormula)
#   workingDataset <- list(predictor = data$predictor[sampleSplitIndex,], treatment = data$treatment[sampleSplitIndex], outcome = data$outcome[sampleSplitIndex])
#
#   robustOutcome_control <- (workingDataset$treatment == FALSE) * (workingDataset$outcome-predictedOutcome$control) / (1-propensity[sampleSplitIndex]) + predictedOutcome$control
#   robustOutcome_treatment <- (workingDataset$treatment == TRUE) * (workingDataset$outcome-predictedOutcome$treatment) / propensity[sampleSplitIndex] + predictedOutcome$treatment
#   robustOutcome <- c(robustOutcome_control, robustOutcome_treatment)
#   pseudoTreatment <- c(-sign(robustOutcome_control), sign(robustOutcome_treatment))
#   pseudoWeight <- c(abs(robustOutcome_control), abs(robustOutcome_treatment))
#   pseudoPredictor <- rbind(workingDataset$predictor, workingDataset$predictor)
#   # We standardize first and set glmnet to not standardize
#   if (loss == 'logistic'){
#     fit <- cv.glmnet(x=pseudoPredictor, y=as.factor(pseudoTreatment), weights=pseudoWeight, family='binomial', intercept = intercept, standardize = FALSE)
#   }
#   list(fit=fit, pseudoPredictor = pseudoPredictor, pseudoWeight = pseudoWeight, pseudoTreatment = pseudoTreatment, sampleSplitIndex=sampleSplitIndex, scaledData=data, propensity = propensity)
# }
#
# # scoreTestAll get the score test for each covariate using all data
# scoreTestAll <- function(itrFit, loss_type='logistic', parallel = TRUE, indexToTest = c(1:8), intercept=FALSE){
#   data <- itrFit$scaledData
#   propensity <- itrFit$propensity
#   link <- predict(itrFit$fit, newx = itrFit$pseudoPredictor, s=itrFit$fit$lambda.min)
#   p <- dim(itrFit$pseudoPredictor)[2]
#   n <- sum(itrFit$sampleSplitIndex)
#   fit_w <- NULL
#   score <- rep(NA, times=length(indexToTest))
#   sigma <- rep(NA, times=length(indexToTest))
#   allPredictor <- rbind(data$predictor, data$predictor)
#   predictedOutcome <- getOutcomeModel(data, method = 'glmnet', sampleSplitIndex = itrFit$sampleSplitIndex, outcomeFormula = NULL, predictAll=TRUE)
#   allOutcome_control <- (data$treatment == FALSE) * (data$outcome-predictedOutcome$control) / (1-propensity) + predictedOutcome$control
#   allOutcome_treatment <- (data$treatment == TRUE) * (data$outcome-predictedOutcome$treatment) / propensity + predictedOutcome$treatment
#   allOutcome <- c(allOutcome_control, allOutcome_treatment)
#   allTreatment <- c(-sign(allOutcome_control), sign(allOutcome_treatment))
#   allWeight <- c(abs(allOutcome_control), abs(allOutcome_treatment))
#   if (!parallel){
#     for (index in indexToTest){
#       pseudoPredictor <- itrFit$pseudoPredictor[,-index]
#       pseudoOutcome <- itrFit$pseudoPredictor[,index]
#       pseudoWeight <- itrFit$pseudoWeight * hessian(itrFit$pseudoTreatment * link,loss_type)
#       fit_w[[index]] <- cv.glmnet(x=pseudoPredictor, y=pseudoOutcome, weights = pseudoWeight, intercept = intercept, standardize = FALSE)
#       link_w <- predict(fit_w[[index]], newx = allPredictor[,-index], s=fit_w[[index]]$lambda.min)
#       # set beta null
#       betaNULL <- array(itrFit$fit$glmnet.fit$beta[,itrFit$fit$lambda==itrFit$fit$lambda.min], c(p,1))
#       betaNULL[index,1] <- 0
#       # get score under null
#       linkNULL <- allPredictor %*% betaNULL + itrFit$fit$glmnet.fit$a0[itrFit$fit$lambda==itrFit$fit$lambda.min]
#       scoreWeight <- allWeight * derivative(allTreatment * linkNULL, loss_type) * allTreatment
#       tmp <- scoreWeight * (allPredictor[,index]-link_w)
#       score[index] <- mean(tmp) * 2
#       sigma[index] <- sqrt(mean((tmp[1:n]+tmp[(n+1):(2*n)])^2))
#     }
#   } else {
#     library(doParallel)
#     n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
#     cl <- makeCluster(min(10, n_cores))
#     registerDoParallel(cl)
#     res <- foreach(index=indexToTest,.packages = 'glmnet') %dopar%{
#       pseudoPredictor <- itrFit$pseudoPredictor[,-index]
#       pseudoOutcome <- itrFit$pseudoPredictor[,index]
#       pseudoWeight <- itrFit$pseudoWeight * hessian(itrFit$pseudoTreatment * link,loss_type)
#       fit_w <- cv.glmnet(x=pseudoPredictor, y=pseudoOutcome, weights = pseudoWeight, intercept = intercept, standardize = FALSE)
#       link_w <- predict(fit_w, newx = allPredictor, s=fit_w$lambda.min)
#       # set beta null
#       betaNULL <- array(itrFit$fit$glmnet.fit$beta[,itrFit$fit$lambda==itrFit$fit$lambda.min],c(p,1))
#       betaNULL[index,1] <- 0
#       # get score under null
#       linkNULL <- allPredictor %*% betaNULL + itrFit$fit$glmnet.fit$a0[itrFit$fit$lambda==itrFit$fit$lambda.min]
#       scoreWeight <- allWeight * derivative(allTreatment * linkNULL, loss_type) * allTreatment
#       tmp <- scoreWeight * (allPredictor[,index]-link_w)
#       score <- mean(tmp) * 2
#       sigma <- sqrt(mean((tmp[1:n]+tmp[(n+1):(2*n)])^2))
#       list(fit_w = fit_w, score=score, sigma=sigma)
#     }
#     stopCluster(cl)
#     for (index in indexToTest){
#       fit_w[[index]] <- res[[index]]$fit_w
#       score[index] <- res[[index]]$score
#       sigma[index] <- res[[index]]$sigma
#     }
#   }
#   list(wFit = fit_w, score = score, sigma=sigma, pvalue=pnorm(-abs(sqrt(n)*score/sigma))*2)
# }

# ITRFitAll obtained the ITR. propensity is defined as p(T=1|X)
ITRFitAll <- function(data, propensity = NULL, outcome = NULL, loss = c('logistic'), sampleSplitIndex=NULL,
                      outcomeModel=c('lm', 'glmnet', 'kernel', 'others'), outcomeFormula = NULL,
                      propensityModel=c('lm', 'glmnet', 'kernel'), propensityFormula = NULL,
                      intercept=FALSE, test=TRUE, parallel=FALSE){
  size <- dim(data$predictor)[1]
  if(is.null(sampleSplitIndex)){
    sampleSplitIndex <- (rnorm(size) > 0)
  }
  fit <- NULL
  fit[[1]] <- ITRFit(data = data, propensity = propensity, outcome = outcome, loss = loss, sampleSplitIndex = sampleSplitIndex,
                     outcomeModel = outcomeModel, outcomeFormula = outcomeFormula, propensityModel = propensityModel,
                     propensityFormula = propensityFormula, intercept = intercept)
  fit[[2]] <- ITRFit(data = data, propensity = propensity, outcome = outcome, loss = loss, sampleSplitIndex = (!sampleSplitIndex),
                     outcomeModel = outcomeModel, outcomeFormula = outcomeFormula, propensityModel = propensityModel,
                     propensityFormula = propensityFormula, intercept = intercept)
  if (test){
    score_1 <- scoreTest(fit[[1]], parallel=parallel)
    score_2 <- scoreTest(fit[[2]], parallel=parallel)
  }
  score <- (score_1$score+ score_2$score)/2
  sigma <- sqrt((score_1$sigma^2 + score_2$sigma^2)/2)
  list(fit =fit, score = score, sigma=sigma, pvalue=pnorm(-abs(sqrt(size)*score/sigma))*2)
}






