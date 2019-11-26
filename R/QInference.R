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
  if (!parallel){
    for (index in indexToTest){
      pseudoPredictor <- qLearnFit$pseudoPredictor[,-index]
      pseudoOutcome <- qLearnFit$pseudoPredictor[,index]
      fit_w[[index]] <- glmnet::cv.glmnet(x=pseudoPredictor, y=pseudoOutcome, intercept = intercept, standardize = TRUE)
      link_w <- predict(fit_w[[index]], newx = pseudoPredictor, s=fit_w[[index]]$lambda.min)
      # set beta null
      betaNULL <- array(qLearnFit$fit$glmnet.fit$beta[,qLearnFit$fit$lambda==qLearnFit$fit$lambda.min], c(p,1))
      betaNULL[index,1] <- 0
      # get score under null
      linkNULL <- qLearnFit$pseudoPredictor %*% betaNULL + qLearnFit$fit$glmnet.fit$a0[qLearnFit$fit$lambda==qLearnFit$fit$lambda.min]
      scoreWeight <- qLearnFit$pseudoOutcome - linkNULL
      tmp <- scoreWeight * (pseudoOutcome-link_w)
      score[index] <- mean(tmp)
      sigma[index] <- sqrt(mean(tmp^2))
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
      betaNULL <- array(qLearnFit$fit$glmnet.fit$beta[,qLearnFit$fit$lambda==qLearnFit$fit$lambda.min],c(p,1))
      betaNULL[index,1] <- 0
      # get score under null
      linkNULL <- qLearnFit$pseudoPredictor %*% betaNULL + qLearnFit$fit$glmnet.fit$a0[qLearnFit$fit$lambda==qLearnFit$fit$lambda.min]
      scoreWeight <- qLearnFit$pseudoOutcome - linkNULL
      tmp <- scoreWeight * (pseudoOutcome-link_w)
      score <- mean(tmp)
      sigma <- sqrt(mean(tmp^2))
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
