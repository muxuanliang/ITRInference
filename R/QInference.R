# QLearnFit obtained the ITR by Q learning method.
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

# scoreTest get the score test for each covariate
scoreTestQLearn <- function(qLearnFit, parallel = TRUE, indexToTest = c(1:8), intercept=FALSE){
  p <- dim(qLearnFit$pseudoPredictor)[2]
  n <- dim(qLearnFit$pseudoPredictor)[1]
  fit_w <- NULL
  score <- rep(NA, times=length(indexToTest))
  sigma <- rep(NA, times=length(indexToTest))
  betaAN <- rep(NA, times=length(indexToTest))
  sigmaAN <- rep(NA, times=length(indexToTest))
  if (!parallel){
    for (index in indexToTest){
      pseudoPredictor <- qLearnFit$pseudoPredictor[,-index]
      pseudoOutcome <- qLearnFit$pseudoPredictor[,index]
      fit_w[[index]] <- glmnet::cv.glmnet(x=pseudoPredictor, y=pseudoOutcome, intercept = intercept, standardize = TRUE)
      link_w <- predict(fit_w[[index]], newx = pseudoPredictor, s=fit_w[[index]]$lambda.min)
      # set beta null
      betaEst <- array(qLearnFit$fit$glmnet.fit$beta[,qLearnFit$fit$lambda==qLearnFit$fit$lambda.min], c(p,1))
      betaNULL <- betaEst
      betaNULL[index,1] <- 0
      # get score under null
      linkNULL <- qLearnFit$pseudoPredictor %*% betaNULL + qLearnFit$fit$glmnet.fit$a0[qLearnFit$fit$lambda==qLearnFit$fit$lambda.min]
      scoreWeight <- qLearnFit$pseudoOutcome - linkNULL
      tmp <- -2*scoreWeight * (pseudoOutcome-link_w)
      score[index] <- mean(tmp)
      # set betaAN
      link <- qLearnFit$pseudoPredictor %*% betaEst + qLearnFit$fit$glmnet.fit$a0[qLearnFit$fit$lambda==qLearnFit$fit$lambda.min]
      scoreWeight <- qLearnFit$pseudoOutcome - link
      tmp <- -2 * scoreWeight * (pseudoOutcome-link_w)
      I <- 2 * pseudoOutcome * (pseudoOutcome - link_w)
      betaAN[index] <- betaEst[index]-mean(tmp)/(mean(I))
      sigma[index] <- sqrt(mean(tmp^2))
      sigmaAN[index] <- sigma[index]/sqrt(mean(I))
    }
  } else {
    library(doParallel)
    n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
    cl <- makeCluster(min(10, n_cores))
    registerDoParallel(cl)
    res <- foreach(index=indexToTest,.packages = 'glmnet') %dopar%{
      pseudoPredictor <- qLearnFit$pseudoPredictor[,-index]
      pseudoOutcome <- qLearnFit$pseudoPredictor[,index]
      fit_w <- glmnet::cv.glmnet(x=pseudoPredictor, y=pseudoOutcome, intercept = intercept, standardize = TRUE)
      link_w <- predict(fit_w, newx = pseudoPredictor, s=fit_w$lambda.min)
      # set beta null
      betaEst <- array(qLearnFit$fit$glmnet.fit$beta[,qLearnFit$fit$lambda==qLearnFit$fit$lambda.min], c(p,1))
      betaNULL <- betaEst
      betaNULL[index,1] <- 0
      # get score under null
      linkNULL <- qLearnFit$pseudoPredictor %*% betaNULL + qLearnFit$fit$glmnet.fit$a0[qLearnFit$fit$lambda==qLearnFit$fit$lambda.min]
      scoreWeight <- qLearnFit$pseudoOutcome - linkNULL
      tmp <- -2 * scoreWeight * (pseudoOutcome-link_w)
      score <- mean(tmp)
      # set betaAN
      link <- qLearnFit$pseudoPredictor %*% betaEst + qLearnFit$fit$glmnet.fit$a0[qLearnFit$fit$lambda==qLearnFit$fit$lambda.min]
      scoreWeight <- qLearnFit$pseudoOutcome - link
      tmp <- -2 * scoreWeight * (pseudoOutcome-link_w)
      I <- 2 * pseudoOutcome * (pseudoOutcome - link_w)
      betaAN <- betaEst[index]-mean(tmp)/(mean(I))
      sigma <- sqrt(mean(tmp^2))
      sigmaAN <- sigma/sqrt(mean(I))
      list(fit_w = fit_w, score=score, sigma=sigma, betaAN=betaAN, sigmaAN=sigmaAN)
    }
    stopCluster(cl)
    for (index in indexToTest){
      fit_w[[index]] <- res[[index]]$fit_w
      score[index] <- res[[index]]$score
      betaAN[index] <- res[[index]]$betaAN
      sigma[index] <- res[[index]]$sigma
      sigmaAN[index] <- res[[index]]$sigmaAN
    }
  }
  list(wFit = fit_w, score = score, sigma=sigma, pvalue=pnorm(-abs(sqrt(n)*score/sigma))*2, betaAN=betaAN, sigmaAN=sigmaAN)
}
