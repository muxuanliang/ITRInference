
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
getOutcomeModel <- function(data, method=c('lm', 'glmnet', 'kernel', 'others'), sampleSplitIndex, outcomeFormula = NULL, predictAll = FALSE){
  p <- dim(data$predictor)[2]
  size <- dim(data$predictor)[1]
  fit <- NULL
  supp <- NULL
  dataPredict <- NULL
  dataPredict$control=data$predictor[sampleSplitIndex,]
  dataPredict$treatment=data$predictor[sampleSplitIndex,]
  if (predictAll){
    dataPredict$control=data$predictor
    dataPredict$treatment=data$predictor
  }
  prediction <- NULL
  dataControl <- list(predictor=data$predictor[(!sampleSplitIndex) & (data$treatment==FALSE),], outcome=data$outcome[(!sampleSplitIndex) & (data$treatment==FALSE)])
  dataTreatment <- list(predictor=data$predictor[(!sampleSplitIndex) & (data$treatment==TRUE),], outcome=data$outcome[(!sampleSplitIndex) & (data$treatment==TRUE)])
  if (0.5*size >= p){
    supp$control <- supp$treatment <- rep(TRUE, times = p)
  }
  if ((0.5*size < p) || (method == 'glmnet')) {
    fit$control <- cv.glmnet(x = dataControl$predictor, y = dataControl$outcome)
    fit$treatment <- cv.glmnet(x = dataTreatment$predictor, y = dataTreatment$outcome)
    supp$control <- abs(fit$control$glmnet.fit$beta[,fit$control$glmnet.fit$lambda==fit$control$lambda.1se])>0
    supp$treatment <- abs(fit$treatment$glmnet.fit$beta[,fit$treatment$glmnet.fit$lambda==fit$treatment$lambda.1se])>0
    dataControl$predictor <- dataControl$predictor[,supp$control]
    dataTreatment$predictor <- dataTreatment$predictor[,supp$treatment]
    prediction$control <- predict(fit$control, newx = dataPredict$control)
    prediction$treatment <- predict(fit$treatment, newx = dataPredict$treatment)
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
  dataPredict <- NULL
  dataPredict$control=data$predictor[sampleSplitIndex,supp$control]
  dataPredict$treatment=data$predictor[sampleSplitIndex,supp$treatment]
  if (predictAll){
    dataPredict$control=data$predictor[,supp$control]
    dataPredict$treatment=data$predictor[,supp$treatment]
  }
  if ((method == 'lm')||(method == 'glmnet')){
    if (sum(supp$control) > 0){
      fit$control <- lm(outcomeFormula(supp$control), data = dataControl)
      prediction$control <- predict(fit$control, newdata = list(predictor=dataPredict$control))
    }
    if (sum(supp$treatment) > 0){
      fit$treatment <- lm(outcomeFormula(supp$treatment), data = dataTreatment)
      prediction$treatment <- predict(fit$treatment, newdata = list(predictor=dataPredict$treatment))
    }
  } else if (method == 'kernel') {
    if (sum(supp$control) > 0){
    prediction$control <- ks(dataControl$predictor, dataControl$outcome, dataPredict$control)
    }
    if (sum(supp$treatment) > 0){
    prediction$treatment[sampleSplitIndex] <- ks(dataTreatment$predictor, dataTreatment$outcome, dataPredict$treatment)
    }
  } else {
    if (sum(supp$control) > 0){
    fit$control <- eval(parse(text=model.gam(dataControl)))
    prediction$control <- predict(fit$control, newdata = list(predictor=dataPredict$control))
    }
    if (sum(supp$treatment) > 0){
    fit$treatment <- eval(parse(text=model.gam(dataTreatment)))
    prediction$treatment <- predict(fit$treatment, newdata = list(predictor=dataPredict$treatment))
    }
  }
  prediction
}
