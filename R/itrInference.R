# ITRFit obtained the ITR. propensity is defined as p(T=1|X)
ITRFit <- function(data, propensity = NULL, outcome = NULL, loss = c('logistic'), sampleSplitIndex=NULL, outcomeModel=c('lm', 'glmnet', 'kernel', 'others'), outcomeFormula = NULL, propensityModel=c('lm', 'glmnet', 'kernel'), propensityFormula = NULL, intercept=FALSE){
  size <- dim(data$predictor)[1]
  if(is.null(sampleSplitIndex)){
    sampleSplitIndex <- (rnorm(size) > 0)
  }
  if (is.null(outcome)){
    predictedOutcome <- getOutcomeModel(data, method = outcomeModel, sampleSplitIndex = sampleSplitIndex, Formula = outcomeFormula)
  } else {
    predictedOutcome <- NULL
    predictedOutcome$control <- outcome$control[sampleSplitIndex]
    predictedOutcome$treatment <- outcome$treatment[sampleSplitIndex]
  }
  if (is.null(propensity)){
    predictedPropensity <- getPropensityModel(data, method = propensityModel, sampleSplitIndex = sampleSplitIndex, Formula = propensityFormula)
  } else {
    predictedPropensity <- propensity[sampleSplitIndex]
  }

  workingDataset <- list(predictor = data$predictor[sampleSplitIndex,], treatment = data$treatment[sampleSplitIndex], outcome = data$outcome[sampleSplitIndex])

  robustOutcome_control <- (workingDataset$treatment == FALSE) * (workingDataset$outcome-predictedOutcome$control) / (1-predictedPropensity) + predictedOutcome$control
  robustOutcome_treatment <- (workingDataset$treatment == TRUE) * (workingDataset$outcome-predictedOutcome$treatment) / predictedPropensity + predictedOutcome$treatment
  robustOutcome <- c(robustOutcome_control, robustOutcome_treatment)
  pseudoTreatment <- c(-sign(robustOutcome_control), sign(robustOutcome_treatment))
  pseudoWeight <- c(abs(robustOutcome_control), abs(robustOutcome_treatment))
  pseudoPredictor <- rbind(workingDataset$predictor, workingDataset$predictor)
  pseudoPredictor <- scale(pseudoPredictor)
  # We standardize first and set glmnet to not standardize
  if (loss == 'logistic'){
    fit <- cv.glmnet(x=pseudoPredictor, y=as.factor(pseudoTreatment), weights=pseudoWeight, family='binomial', intercept = intercept, standardize = FALSE)
  }
  list(fit=fit, pseudoPredictor = pseudoPredictor, pseudoWeight = pseudoWeight, pseudoTreatment = pseudoTreatment, sampleSplitIndex=sampleSplitIndex)
}

# scoreTest get the score test for each covariate
scoreTest <- function(itrFit, loss_type='logistic', parallel = TRUE, indexToTest = c(1:8), intercept=FALSE){
  link <- predict(itrFit$fit, newx = itrFit$pseudoPredictor, s=itrFit$fit$lambda.min)
  p <- dim(itrFit$pseudoPredictor)[2]
  n <- sum(itrFit$sampleSplitIndex)
  fit_w <- NULL
  score <- rep(NA, times=length(indexToTest))
  sigma <- rep(NA, times=length(indexToTest))
  if (!parallel){
    for (index in indexToTest){
      pseudoPredictor <- itrFit$pseudoPredictor[,-index]
      pseudoOutcome <- itrFit$pseudoPredictor[,index]
      pseudoWeight <- itrFit$pseudoWeight * hessian(itrFit$pseudoTreatment * link,loss_type)
      fit_w[[index]] <- cv.glmnet(x=pseudoPredictor, y=pseudoOutcome, weights = pseudoWeight, intercept = intercept, standardize = FALSE)
      link_w <- predict(fit_w[[index]], newx = pseudoPredictor, s=fit_w[[index]]$lambda.min)
      # set beta null
      betaNULL <- array(itrFit$fit$glmnet.fit$beta[,itrFit$fit$lambda==itrFit$fit$lambda.min], c(p,1))
      betaNULL[index,1] <- 0
      # get score under null
      linkNULL <- itrFit$pseudoPredictor %*% betaNULL + itrFit$fit$glmnet.fit$a0[itrFit$fit$lambda==itrFit$fit$lambda.min]
      scoreWeight <- itrFit$pseudoWeight * derivative(itrFit$pseudoTreatment * linkNULL, loss_type) * itrFit$pseudoTreatment
      tmp <- scoreWeight * (pseudoOutcome-link_w)
      score[index] <- mean(tmp) * 2
      sigma[index] <- sqrt(mean((tmp[1:n]+tmp[(n+1):(2*n)])^2))
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
      fit_w <- cv.glmnet(x=pseudoPredictor, y=pseudoOutcome, weights = pseudoWeight, intercept = FALSE, standardize = FALSE)
      link_w <- predict(fit_w, newx = pseudoPredictor, s=fit_w$lambda.min)
      # set beta null
      betaNULL <- array(itrFit$fit$glmnet.fit$beta[,itrFit$fit$lambda==itrFit$fit$lambda.min],c(p,1))
      betaNULL[index,1] <- 0
      # get score under null
      linkNULL <- itrFit$pseudoPredictor %*% betaNULL + itrFit$fit$glmnet.fit$a0[itrFit$fit$lambda==itrFit$fit$lambda.min]
      scoreWeight <- itrFit$pseudoWeight * derivative(itrFit$pseudoTreatment * linkNULL,loss_type) * itrFit$pseudoTreatment
      tmp <- scoreWeight * (pseudoOutcome-link_w)
      score <- mean(tmp) * 2
      sigma <- sqrt(mean((tmp[1:n]+tmp[(n+1):(2*n)])^2))
      list(fit_w = fit_w, score=score, sigma=sigma)
    }
    stopCluster(cl)
    for (index in indexToTest){
      fit_w[[index]] <- res[[index]]$fit_w
      score[index] <- res[[index]]$score
      sigma[index] <- res[[index]]$sigma
    }
  }
  list(wFit = fit_w, score = score, sigma=sigma, pvalue=pnorm(-abs(sqrt(n)*score/sigma))*2)
}