# ITRFit obtained the ITR. propensity is defined as p(T=1|X)
ITRFit <- function(data, propensity = NULL, outcome = NULL, loss = 'logistic', sampleSplitIndex=NULL, type.measure = 'lossFun', outcomeModel=c('lm', 'glmnet', 'kernel', 'others'), outcomeFormula = NULL, propensityModel=c('lm', 'glmnet', 'kernel'), propensityFormula = NULL, intercept=FALSE, screeningMethod = "SIRS", outcomeScreeningFamily = 'Gaussian', standardize = TRUE){
  size <- dim(data$predictor)[1]
  if(is.null(sampleSplitIndex)){
    sampleSplitIndex <- (rnorm(size) > 0)
  }
  if (is.null(outcome)){
    predictedOutcomeAll <- getOutcomeModel(data, method = outcomeModel, sampleSplitIndex = sampleSplitIndex, Formula = outcomeFormula, predictAll = TRUE, screeningMethod = screeningMethod, outcomeScreeningFamily = outcomeScreeningFamily)
  } else {
    predictedOutcomeAll <- NULL
    predictedOutcomeAll$control <- outcome$control
    predictedOutcomeAll$treatment <- outcome$treatment
  }
  predictedOutcome <- NULL
  predictedOutcome$control <- predictedOutcomeAll$control[sampleSplitIndex]
  predictedOutcome$treatment <- predictedOutcomeAll$treatment[sampleSplitIndex]
  if (is.null(propensity)){
    predictedPropensityAll <- getPropensityModel(data, method = propensityModel, sampleSplitIndex = sampleSplitIndex, Formula = propensityFormula, predictAll = TRUE, screeningMethod = screeningMethod)
  } else {
    predictedPropensityAll <- propensity
  }
  predictedPropensity <- predictedPropensityAll[sampleSplitIndex]
  
  workingDataset <- list(predictor = data$predictor[sampleSplitIndex,], treatment = data$treatment[sampleSplitIndex], outcome = data$outcome[sampleSplitIndex])
  
  robustOutcome_control <- (workingDataset$treatment == FALSE) * (workingDataset$outcome-predictedOutcome$control) / (1-predictedPropensity) + predictedOutcome$control
  robustOutcome_treatment <- (workingDataset$treatment == TRUE) * (workingDataset$outcome-predictedOutcome$treatment) / predictedPropensity + predictedOutcome$treatment
  robustOutcome <- c(robustOutcome_control, robustOutcome_treatment)
  pseudoTreatment <- c(-sign(robustOutcome_control), sign(robustOutcome_treatment))
  pseudoWeight <- c(abs(robustOutcome_control), abs(robustOutcome_treatment))
  pseudoPredictor <- rbind(workingDataset$predictor, workingDataset$predictor)
  if(standardize){
    pseudoPredictor <- scale(pseudoPredictor)
  }
  
  # We standardize first and set glmnet to not standardize
  if (loss == 'logistic'){
    fitInit <- glmnet::glmnet(x=pseudoPredictor, y=as.numeric(pseudoTreatment==1), weights=pseudoWeight, family='binomial', intercept = intercept, standardize = FALSE)
    lambdaInit <- max(fitInit$lambda)
    lambdaSeq <- lambdaInit * (10/11)^(0:99)
    fit <- fitInit #glmnet::glmnet(x=pseudoPredictor, y=as.numeric(pseudoTreatment==1), weights=pseudoWeight, family='binomial', intercept = intercept, standardize = FALSE, lambda = lambdaSeq)
    validate <- predict(fit, newx = data$predictor[!sampleSplitIndex,])
    if (type.measure=='lossFun'){
      cvm <- apply(validate, 2, function(t){
        predictedOutcome <- NULL
        predictedOutcome$control <- predictedOutcomeAll$control[!sampleSplitIndex]
        predictedOutcome$treatment <- predictedOutcomeAll$treatment[!sampleSplitIndex]
        predictedPropensity <- predictedPropensityAll[!sampleSplitIndex]
        workingDataset <- list(predictor = data$predictor[!sampleSplitIndex,], treatment = data$treatment[!sampleSplitIndex], outcome = data$outcome[!sampleSplitIndex])
        robustOutcome_control <- (workingDataset$treatment == FALSE) * (workingDataset$outcome-predictedOutcome$control) / (1-predictedPropensity) + predictedOutcome$control
        robustOutcome_treatment <- (workingDataset$treatment == TRUE) * (workingDataset$outcome-predictedOutcome$treatment) / predictedPropensity + predictedOutcome$treatment
        robustOutcome <- c(robustOutcome_control, robustOutcome_treatment)
        pseudoTreatment <- c(-sign(robustOutcome_control), sign(robustOutcome_treatment))
        pseudoWeight <- c(abs(robustOutcome_control), abs(robustOutcome_treatment))
        pseudoLink <- c(t, t)
        score <- mean(pseudoWeight*loss(pseudoTreatment * pseudoLink, loss_type='logistic'))
        score
      })
      fit$lambda.min <- fit$lambda[which.min(cvm)]
    } else {
      cvm <- apply(validate, 2, function(t){
        weight <- predictedPropensityAll[!sampleSplitIndex]*(t>=0) + (1-predictedPropensityAll[!sampleSplitIndex])*(t<0)
        predictedOutcome <- NULL
        predictedOutcome$control <- predictedOutcomeAll$control[!sampleSplitIndex]
        predictedOutcome$treatment <- predictedOutcomeAll$treatment[!sampleSplitIndex]
        augment <- predictedOutcome$control*(t<0) + predictedOutcome$treatment * (t>=0)
        score <- mean(data$outcome[!sampleSplitIndex]/weight * (t * (data$treatment[!sampleSplitIndex]-0.5)>0) + augment)
        score
      })
      fit$lambda.min <- fit$lambda[which.max(cvm)]
    }
  }
  list(fit=fit, pseudoPredictor = pseudoPredictor, pseudoWeight = pseudoWeight, pseudoTreatment = pseudoTreatment, sampleSplitIndex=sampleSplitIndex)
}

# scoreTest get the score test for each covariate
scoreTest <- function(itrFit, loss_type='logistic', parallel = TRUE, indexToTest = c(1:8), intercept=FALSE){
  p <- dim(itrFit$pseudoPredictor)[2]
  n <- sum(itrFit$sampleSplitIndex)
  link <- predict(itrFit$fit, newx = itrFit$pseudoPredictor, s=itrFit$fit$lambda.min)
  fit_w <- NULL
  score <- rep(NA, times=length(indexToTest))
  sigma <- rep(NA, times=length(indexToTest))
  betaAN <- rep(NA, times=length(indexToTest))
  I <- rep(NA, times=length(indexToTest))
  if (!parallel){
    for (index in indexToTest){
      pseudoPredictor <- itrFit$pseudoPredictor[,-index]
      pseudoOutcome <- itrFit$pseudoPredictor[,index]
      pseudoWeight <- itrFit$pseudoWeight * hessian(itrFit$pseudoTreatment * link,loss_type)
      fit_w[[index]] <- glmnet::cv.glmnet(x=pseudoPredictor, y=pseudoOutcome, weights = pseudoWeight, intercept = intercept, standardize = FALSE)
      link_w <- predict(fit_w[[index]], newx = pseudoPredictor, s=fit_w[[index]]$lambda.min)
      # set beta null
      betaEst <- array(itrFit$fit$beta[,itrFit$fit$lambda==itrFit$fit$lambda.min], c(p,1))
      betaNULL <- betaEst
      betaNULL[index,1] <- 0
      # get score under null
      linkNULL <- itrFit$pseudoPredictor %*% betaNULL + itrFit$fit$a0[itrFit$fit$lambda==itrFit$fit$lambda.min]
      scoreWeightNULL <- itrFit$pseudoWeight * derivative(itrFit$pseudoTreatment * linkNULL, loss_type) * itrFit$pseudoTreatment
      tmpNULL <- scoreWeightNULL * (pseudoOutcome-link_w)
      score[index] <- mean(tmpNULL) * 2
      # set betaAN
      link <- itrFit$pseudoPredictor %*% betaEst + itrFit$fit$a0[itrFit$fit$lambda==itrFit$fit$lambda.min]
      scoreWeight <- itrFit$pseudoWeight * derivative(itrFit$pseudoTreatment * link, loss_type) * itrFit$pseudoTreatment
      tmp <- scoreWeight * (pseudoOutcome-link_w)
      Itmp <- itrFit$pseudoWeight * hessian(itrFit$pseudoTreatment * link, loss_type) * pseudoOutcome * (pseudoOutcome - link_w)
      betaAN[index] <- betaEst[index]-mean(tmp) * 2/(mean(Itmp)*2)
      sigma[index] <- sqrt(mean((tmp[1:n]+tmp[(n+1):(2*n)])^2))
      I[index] <- (mean(Itmp)*2)
      
      sigma[index] <- sqrt(mean((tmpNULL[1:n]+tmpNULL[(n+1):(2*n)])^2))
    }
  } else {
    library(doParallel)
    n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
    cl <- makeCluster(min(10, n_cores))
    registerDoParallel(cl)
    res <- foreach(index=indexToTest,.packages = 'glmnet') %dopar%{
      pseudoPredictor <- itrFit$pseudoPredictor[,-index]
      pseudoOutcome <- itrFit$pseudoPredictor[,index]
      pseudoWeight <- itrFit$pseudoWeight * hessian(itrFit$pseudoTreatment * link,loss_type)
      fit_w <- glmnet::cv.glmnet(x=pseudoPredictor, y=pseudoOutcome, weights = pseudoWeight, intercept = intercept, standardize = FALSE)
      link_w <- predict(fit_w, newx = pseudoPredictor, s=fit_w$lambda.min)
      # set beta null
      betaEst <- array(itrFit$fit$beta[,itrFit$fit$lambda==itrFit$fit$lambda.min], c(p,1))
      betaNULL <- betaEst
      betaNULL[index,1] <- 0
      # get score under null
      linkNULL <- itrFit$pseudoPredictor %*% betaNULL + itrFit$fit$a0[itrFit$fit$lambda==itrFit$fit$lambda.min]
      scoreWeightNULL <- itrFit$pseudoWeight * derivative(itrFit$pseudoTreatment * linkNULL, loss_type) * itrFit$pseudoTreatment
      tmpNULL <- scoreWeightNULL * (pseudoOutcome-link_w)
      score <- mean(tmpNULL) * 2
      # get betaAN
      link <- itrFit$pseudoPredictor %*% betaEst + itrFit$fit$a0[itrFit$fit$lambda==itrFit$fit$lambda.min]
      scoreWeight <- itrFit$pseudoWeight * derivative(itrFit$pseudoTreatment * link, loss_type) * itrFit$pseudoTreatment
      tmp <- scoreWeight * (pseudoOutcome-link_w)
      Itmp <- itrFit$pseudoWeight * hessian(itrFit$pseudoTreatment * link, loss_type) * pseudoOutcome * (pseudoOutcome - link_w)
      betaAN <- betaEst[index]-mean(tmp) * 2/(mean(Itmp)*2)
      sigma <- sqrt(mean((tmp[1:n]+tmp[(n+1):(2*n)])^2))
      I <- (mean(Itmp)*2)
      
      sigma <- sqrt(mean((tmpNULL[1:n]+tmpNULL[(n+1):(2*n)])^2))
      list(fit_w = fit_w, score=score, sigma=sigma, betaAN=betaAN, I=I)
    }
    stopCluster(cl)
    for (index in indexToTest){
      score[index] <- res[[index]]$score
      sigma[index] <- res[[index]]$sigma
      betaAN[index] <- res[[index]]$betaAN
      I[index] <- res[[index]]$I
    }
  }
  list(score = score, sigma=sigma, pvalue=pnorm(-abs(sqrt(n)*score/sigma))*2, betaAN=betaAN, I=I)
}