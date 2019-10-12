
# loss gets the loss of the chosen loss
loss <- function(x, loss_type){
  switch(loss_type,
         logistic = log(1+exp(-x)),
         exponential = exp(-x),
         zeroOne = (x<0)
  )
}

# derivative gets the derivative of the chosen loss
derivative <- function(x, loss_type){
  switch(loss_type,
         logistic = -exp(-x)/(1+exp(-x)),
         exponential = -exp(-x)
  )
}

# hessian gets the hessian of the chosen loss
hessian <- function(x, loss_type){
  switch(loss_type,
         logistic = exp(-x)/(1+exp(-x))^2,
         exponential = exp(-x)
  )
}

# ks gets the kernel estimation
ks <- function(xx, yy, xx.test){
  nobs <- nrow(xx)
  nvars <- ncol(xx)
  hopt <- (4/(nvars+2))^(1/(nvars+4)) * (nobs^(-1/(nvars+4)))
  wm <- function(t){
    if (ncol(xx)==1){
      weight <- exp(-0.5 * (as.numeric(t)-xx)^2/(hopt^2)) * hopt
    } else {
      weight <- apply(xx,1,function(s){exp(-0.5 * sum((t-s)^2)/(hopt^2)) * hopt^(ncol(xx))})
    }
    weighted.mean(yy, weight)
  }
  if (nrow(as.matrix(xx.test))==1) {
    yy.test <- wm(xx.test)
  } else {
    yy.test<- array(0,c(nrow(xx.test),1))
    for (index in 1:nrow(xx.test)){
      if (ncol(xx.test)==1){
        yy.test[index] <- wm(xx.test[index])
      } else {
        yy.test[index] <- wm(xx.test[index,])
      }
    }
  }
  yy.test
}

# model.gam gets the gam regression function given data
model.gam <- function(data){
  p <- dim(data$predictor)[2]
  expr <- "gam(outcome~"
  for (i in 1:(p-1)){
    expr <- paste0(expr, "s(predictor[i])+")
  }
  expr <- paste0(expr, "s(predictor[p]),method ='REML'")
  expr
}

# getOutcomeModel contains the outcome regression model
getOutcomeModel <- function(data, method=c('lm', 'glmnet', 'kernel', 'others'), sampleSplitIndex, outcomeFormula = NULL){
  p <- dim(data$predictor)[2]
  size <- dim(data$predictor)[1]
  fit <- NULL
  supp <- NULL
  prediction <- NULL
  dataControl <- list(predictor=data$predictor[(!sampleSplitIndex) & (data$treatment==FALSE),], outcome=data$outcome[(!sampleSplitIndex) & (data$treatment==FALSE)])
  dataTreatment <- list(predictor=data$predictor[(!sampleSplitIndex) & (data$treatment==TRUE),], outcome=data$outcome[(!sampleSplitIndex) & (data$treatment==TRUE)])
  if (0.5*size >= p){
    supp$control <- supp$treatment <- rep(TRUE, times = p)
  }
  if ((0.5*size < p) || (method == 'glmnet')) {
    fit$control <- cv.glmnet(x = dataControl$predictor, y = dataControl$outcome)
    fit$treatment <- cv.glmnet(x = dataTreatment$predictor, y = dataTreatment$outcome)
    prediction$control <- predict(fit$control, newx = data$predictor[sampleSplitIndex,], s=fit$control$lambda.min)
    prediction$treatment <- predict(fit$treatment, newx = data$predictor[sampleSplitIndex,], s=fit$treatment$lambda.min)
    supp$control <- abs(fit$control$glmnet.fit$beta[,fit$control$glmnet.fit$lambda==fit$control$lambda.min])>0
    supp$treatment <- abs(fit$treatment$glmnet.fit$beta[,fit$treatment$glmnet.fit$lambda==fit$treatment$lambda.min])>0
    dataControl$predictor <- dataControl$predictor[,supp$control]
    dataTreatment$predictor <- dataTreatment$predictor[,supp$treatment]
  }

  if (is.null(outcomeFormula)){
    outcomeFormula <- function(support){
      expr <- (outcome ~ predictor)
      if (sum(support)==0){
        expr <- (outcome ~ 1)
      }
    expr
    }
  }
  fit <- NULL
  if ((method == 'lm')||(method == 'glmnet')){
    fit$control <- lm(outcomeFormula(supp$control), data = dataControl)
    fit$treatment <- lm(outcomeFormula(supp$treatment), data = dataTreatment)
    prediction$control <- predict(fit$control, newdata = list(predictor=data$predictor[sampleSplitIndex,supp$control]))
    prediction$treatment <- predict(fit$treatment, newdata = list(predictor=data$predictor[sampleSplitIndex,supp$treatment]))
  } else if (method == 'kernel') {
    prediction$control[sampleSplitIndex] <- ks(dataControl$predictor, dataControl$outcome, data$predictor[sampleSplitIndex,supp$control])
    prediction$treatment[sampleSplitIndex] <- ks(dataTreatment$predictor, dataTreatment$outcome, data$predictor[sampleSplitIndex,supp$control])
  } else {
    fit$control <- eval(parse(text=model.gam(dataControl)))
    fit$treatment <- eval(parse(text=model.gam(dataTreatment)))
    prediction$control[sampleSplitIndex] <- predict(fit$control, newdata = list(predictor=data$predictor[sampleSplitIndex,supp$control]))
    prediction$treatment[sampleSplitIndex] <- predict(fit$treatment, newdata = list(predictor=data$predictor[sampleSplitIndex,supp$treatment]))
  }
  prediction
}

# ITRFit obtained the ITR. propensity is defined as p(T=1|X)
ITRFit <- function(data, propensity, loss = c('logistic'), outcomeModel=c('lm', 'glmnet', 'kernel', 'others'), outcomeFormula = NULL, intercept=FALSE){
  size <- dim(data$predictor)[1]
  sampleSplitIndex <- (rnorm(size) > 0)
  predictedOutcome <- getOutcomeModel(data, method = outcomeModel, sampleSplitIndex = sampleSplitIndex, outcomeFormula = outcomeFormula)
  workingDataset <- list(predictor = data$predictor[sampleSplitIndex,], treatment = data$treatment[sampleSplitIndex], outcome = data$outcome[sampleSplitIndex])

  robustOutcome_control <- (workingDataset$treatment == FALSE) * (workingDataset$outcome-predictedOutcome$control) / (1-propensity[sampleSplitIndex]) + predictedOutcome$control
  robustOutcome_treatment <- (workingDataset$treatment == TRUE) * (workingDataset$outcome-predictedOutcome$treatment) / propensity[sampleSplitIndex] + predictedOutcome$treatment
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
      betaNULL <- itrFit$fit$glmnet.fit$beta[,itrFit$fit$lambda==itrFit$fit$lambda.min]
      betaNULL[index] <- 0
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
      betaNULL <- itrFit$fit$glmnet.fit$beta[,itrFit$fit$lambda==itrFit$fit$lambda.min]
      betaNULL[index] <- 0
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
