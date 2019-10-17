# QLearnFit obtained the ITR by Q learning method.
QLearnFit <- function(data, intercept=FALSE){
  size <- dim(data$predictor)[1]
  Outcome <- data$outcome
  Treatment <- (data$treatment - 0.5) * 2
  pseudoPredictor <- cbind(data$predictor, apply(data$predictor,2,function(t){t*Treatment}))
  fit <- cv.glmnet(x=pseudoPredictor, y=Outcome, weights=pseudoWeight, family='gaussian', intercept = intercept, standardize = TRUE)
  list(fit=fit, pseudoPredictor = pseudoPredictor, pseudoTreatment = Treatment, pseudoOutcome = Outcome)
}

# scoreTest get the score test for each covariate
scoreTest <- function(itrFit, loss_type='logistic', parallel = TRUE){
  link <- predict(itrFit$fit, newx = itrFit$pseudoPredictor, s=itrFit$fit$lambda.min)
  p <- dim(itrFit$pseudoPredictor)[2]
  n <- sum(itrFit$sampleSplitIndex)
  fit_w <- NULL
  score <- rep(NA, times=p)
  sigma <- rep(NA, times=p)
  if (!parallel){
    for (index in 1:p){
      pseudoPredictor <- itrFit$pseudoPredictor[,-index]
      pseudoOutcome <- itrFit$pseudoPredictor[,index]
      pseudoWeight <- itrFit$pseudoWeight * hessian(itrFit$pseudoTreatment * link,loss_type)
      fit_w[[index]] <- cv.glmnet(x=pseudoPredictor, y=pseudoOutcome, weights = pseudoWeight, intercept = FALSE, standardize = FALSE)
      w_est <- fit_w[[index]]$glmnet.fit$beta[,fit_w[[index]]$lambda==fit_w[[index]]$lambda.min]
      # set beta null
      betaNULL <- array(itrFit$fit$glmnet.fit$beta[,itrFit$fit$lambda==itrFit$fit$lambda.min], c(p,1))
      betaNULL[index,1] <- 0
      linkNULL <- itrFit$pseudoPredictor %*% betaNULL + itrFit$fit$glmnet.fit$a0[itrFit$fit$lambda==itrFit$fit$lambda.min]
      scoreWeight <- itrFit$pseudoWeight * derivative(itrFit$pseudoTreatment * linkNULL, loss_type) * itrFit$pseudoTreatment
      tmp <- scoreWeight * (pseudoOutcome-pseudoPredictor%*%w_est)
      score[index] <- mean(tmp) * 2
      sigma[index] <- sqrt(mean((tmp[1:n]+tmp[(n+1):(2*n)])^2))
    }
  } else {
    library(doParallel)
    n_cores <- detectCores(all.tests = FALSE, logical = TRUE)
    cl <- makeCluster(min(10, n_cores))
    registerDoParallel(cl)
    res <- foreach(index=1:p,.packages = 'glmnet') %dopar%{
      pseudoPredictor <- itrFit$pseudoPredictor[,-index]
      pseudoOutcome <- itrFit$pseudoPredictor[,index]
      pseudoWeight <- itrFit$pseudoWeight * hessian(itrFit$pseudoTreatment * link,loss_type)
      fit_w <- cv.glmnet(x=pseudoPredictor, y=pseudoOutcome, weights = pseudoWeight, intercept = FALSE, standardize = FALSE)
      w_est <- fit_w$glmnet.fit$beta[,fit_w$lambda==fit_w$lambda.min]
      # set beta null
      betaNULL <- array(itrFit$fit$glmnet.fit$beta[,itrFit$fit$lambda==itrFit$fit$lambda.min],c(p,1))
      betaNULL[index,1] <- 0
      linkNULL <- itrFit$pseudoPredictor %*% betaNULL + itrFit$fit$glmnet.fit$a0[itrFit$fit$lambda==itrFit$fit$lambda.min]
      scoreWeight <- itrFit$pseudoWeight * derivative(itrFit$pseudoTreatment * linkNULL,loss_type) * itrFit$pseudoTreatment
      tmp <- scoreWeight * (pseudoOutcome-pseudoPredictor%*%w_est)
      score <- mean(tmp) * 2
      sigma <- sqrt(mean((tmp[1:n]+tmp[(n+1):(2*n)])^2))
      list(fit_w = fit_w, score=score, sigma=sigma)
    }
    stopCluster(cl)
    for (index in 1:p){
      fit_w[[index]] <- res[[index]]$fit_w
      score[index] <- res[[index]]$score
      sigma[index] <- res[[index]]$sigma
    }
  }
  list(wFit = fit_w, score = score, sigma=sigma, pvalue=pnorm(-abs(sqrt(n)*score/sigma))*2)
}
