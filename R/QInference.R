#' Estimation and inference for individualized treatment rule using pooled-and-split de-correlated score
#' 
#' This function implements the Q-learning estimation for individualized treatment rule and the inference procedure based on the de-correlated score (see reference).
#'
#' @param qLearnFit returns of \code{QLearnFit}
#' @param intercept includes intercept or not
#' @param indexToTest indicates which coefficients to test. By default, c(1:8)
#' @param parallel whether use parallel computing; by default, FALSE.
#' @return p-values are the p-value for each coefficients included in indexToTest. (betaAN-1.96*sigmaAN/sqrt(sample size), betaAN+1.96*sigmaAN/sqrt(sample size)) provides the 95\% CI for the coefficients.
#' @author Muxuan Liang <mliang@fredhutch.org>
#' @references Muxuan Liang, Young-Geun Choi, Yang Ning, Maureen Smith, Yingqi Zhao (2020). Estimation and inference on high-dimensional individualized treatment rule in observational data using split-and-pooled de-correlated score.
#' @export
QFitInfer <- function(qLearnFit, parallel = TRUE, indexToTest = c(1:8), intercept=TRUE){
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
      tmpNULL <- -2*scoreWeight * (pseudoOutcome-link_w)
      score[index] <- mean(tmpNULL)
      # set betaAN
      link <- qLearnFit$pseudoPredictor %*% betaEst + qLearnFit$fit$glmnet.fit$a0[qLearnFit$fit$lambda==qLearnFit$fit$lambda.min]
      scoreWeight <- qLearnFit$pseudoOutcome - link
      tmp <- -2 * scoreWeight * (pseudoOutcome-link_w)
      I <- 2 * pseudoOutcome * (pseudoOutcome - link_w)
      betaAN[index] <- betaEst[index]-mean(tmp)/(mean(I))
      sigma[index] <- sqrt(mean(tmp^2))
      sigmaAN[index] <- sigma[index]/(mean(I))

      sigma[index] <- sqrt(mean(tmpNULL^2))
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
      tmpNULL <- -2 * scoreWeight * (pseudoOutcome-link_w)
      score <- mean(tmpNULL)
      # set betaAN
      link <- qLearnFit$pseudoPredictor %*% betaEst + qLearnFit$fit$glmnet.fit$a0[qLearnFit$fit$lambda==qLearnFit$fit$lambda.min]
      scoreWeight <- qLearnFit$pseudoOutcome - link
      tmp <- -2 * scoreWeight * (pseudoOutcome-link_w)
      I <- 2 * pseudoOutcome * (pseudoOutcome - link_w)
      betaAN <- betaEst[index]-mean(tmp)/(mean(I))
      sigma <- sqrt(mean(tmp^2))
      sigmaAN <- sigma/(mean(I))
      sigma <- sqrt(mean(tmpNULL^2))
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
