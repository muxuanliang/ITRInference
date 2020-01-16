
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
  nobs <- nrow(as.matrix(xx))
  nvars <- ncol(as.matrix(xx))
  hopt <- (4/(nvars+2))^(1/(nvars+4)) * (nobs^(-1/(nvars+4)))
  wm <- function(t){
    if (ncol(as.matrix(xx))==1){
      weight <- exp(-0.5 * (as.numeric(t)-xx)^2/(hopt^2)) * hopt
    } else {
      weight <- apply(xx,1,function(s){exp(-0.5 * sum((t-s)^2)/(hopt^2)) * hopt^(ncol(xx))})
    }
    weighted.mean(yy, weight)
  }
  if (nrow(as.matrix(xx.test))==1) {
    yy.test <- wm(xx.test)
  } else {
    if (ncol((as.matrix(xx.test)))==1){
      yy.test <- sapply(as.matrix(xx.test), function(t){
        wm(t)
      })
    } else {
      yy.test <- apply(as.matrix(xx.test),1,function(t){
        wm(t)
      })
    }
  }
  yy.test
}

# model.gam gets the gam regression function given data
model.gam <- function(data){
  p <- dim(data$predictor)[2]
  expr <- "mgcv::gam(outcome~"
  for (i in 1:(p-1)){
    expr <- paste0(expr, "s(predictor[,",i,"])+")
  }
  expr <- paste0(expr, "s(predictor[,",p,"]), data = data,method ='REML')")
  expr
}

# getOutcomeModel contains the outcome regression model
getOutcomeModel <- function(data, method=c('lm', 'glmnet', 'kernel', 'others'), sampleSplitIndex, Formula = NULL, predictAll = FALSE, screeningMethod="SIRS", outcomeScreeningFamily='Gaussian'){
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
  if (0.05*size >= p){
    supp$control <- supp$treatment <- rep(TRUE, times = p)
    if ((method != 'lm')&&(method != 'glmnet')){
      ans1 <- screening(data$predictor[data$treatment==FALSE,], data$outcome[data$treatment==FALSE], method = screeningMethod, family = outcomeScreeningFamily)
      ans2 <- screening(data$predictor[data$treatment==TRUE,], data$outcome[data$treatment==TRUE], method = screeningMethod, family = outcomeScreeningFamily)
      if (screeningMethod == 'glmnet'){
        supp$control <- ans1
        supp$treatment <- ans2
      } else {
        supp$control <- supp$treatment <- (ans1 <= floor(p/2))|(ans2 <= floor(p/2))
      }
      dataControl$predictor <- dataControl$predictor[,supp$control]
      dataTreatment$predictor <- dataTreatment$predictor[,supp$treatment]
    }
  }
  if ((0.05*size < p) || (method == 'glmnet')) {
    fit$control <- glmnet::cv.glmnet(x = dataControl$predictor, y = dataControl$outcome)
    fit$treatment <- glmnet::cv.glmnet(x = dataTreatment$predictor, y = dataTreatment$outcome)
    supp$control <- abs(fit$control$glmnet.fit$beta[,fit$control$glmnet.fit$lambda==fit$control$lambda.min])>0
    supp$treatment <- abs(fit$treatment$glmnet.fit$beta[,fit$treatment$glmnet.fit$lambda==fit$treatment$lambda.min])>0
    if ((method != 'lm')&&(method != 'glmnet')){
      ans1 <- screening(data$predictor[data$treatment==FALSE,], data$outcome[data$treatment==FALSE], method = screeningMethod, family = outcomeScreeningFamily)
      ans2 <- screening(data$predictor[data$treatment==TRUE,], data$outcome[data$treatment==TRUE], method = screeningMethod, family = outcomeScreeningFamily)
      if (screeningMethod == 'glmnet'){
        supp$control <- ans1
        supp$treatment <- ans2
      } else {
        supp$control <- supp$treatment <- (ans1 <= 10)|(ans2 <= 10)
      }
    }
    dataControl$predictor <- dataControl$predictor[,supp$control]
    dataTreatment$predictor <- dataTreatment$predictor[,supp$treatment]
    prediction$control <- predict(fit$control, newx = dataPredict$control, s=fit$control$lambda.min)
    prediction$treatment <- predict(fit$treatment, newx = dataPredict$treatment, s=fit$treatment$lambda.min)
  }

  if (is.null(Formula)){
    Formula <- function(support){
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
      fit$control <- lm(Formula(supp$control), data = dataControl)
      prediction$control <- predict(fit$control, newdata = list(predictor=dataPredict$control))
    }
    if (sum(supp$treatment) > 0){
      fit$treatment <- lm(Formula(supp$treatment), data = dataTreatment)
      prediction$treatment <- predict(fit$treatment, newdata = list(predictor=dataPredict$treatment))
    }
  } else if (method == 'kernel') {
    if (sum(supp$control) > 0){
    prediction$control <- ks(dataControl$predictor, dataControl$outcome, dataPredict$control)
    }
    if (sum(supp$treatment) > 0){
    prediction$treatment <- ks(dataTreatment$predictor, dataTreatment$outcome, dataPredict$treatment)
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

# getPropensityModel contains the outcome regression model
getPropensityModel <- function(data, method=c('lm', 'glmnet', 'kernel'), sampleSplitIndex, Formula = NULL, predictAll = FALSE, screeningMethod="SIRS"){
  p <- dim(data$predictor)[2]
  size <- dim(data$predictor)[1]
  fit <- NULL
  supp <- NULL
  dataPredict <- NULL
  dataPredict=data$predictor[sampleSplitIndex,]
  if (predictAll){
    dataPredict=data$predictor
  }
  prediction <- NULL
  dataTrain <- list(predictor=data$predictor[(!sampleSplitIndex),], treatment=data$treatment[(!sampleSplitIndex)])
  if (0.05*size >= p){
    supp$control <- supp$treatment <- rep(TRUE, times = p)
    if ((method != 'lm')&&(method != 'glmnet')){
      ans <- screening(data$predictor, data$treatment, method = screeningMethod, family = 'binomial')
      if (screeningMethod == 'glmnet'){
        supp <- ans
      } else {
        supp <- (ans <= p/2)
      }
    }
  }
  if ((0.05*size < p) || (method == 'glmnet')) {
    fit <- glmnet::cv.glmnet(x = dataTrain$predictor, y = dataTrain$treatment, family='binomial')
    supp <- abs(fit$glmnet.fit$beta[,fit$glmnet.fit$lambda==fit$lambda.min])>0
    if ((method != 'lm')&&(method != 'glmnet')){
      ans <- screening(data$predictor, data$treatment, method = screeningMethod, family = 'binomial')
      if (screeningMethod == 'glmnet'){
        supp <- ans
      } else {
        supp <- (ans <= 5)
      }
    }
    dataTrain$predictor <- dataTrain$predictor[,supp]
    prediction <- predict(fit, newx = dataPredict, type='response', s=fit$lambda.min)
  }

  if (is.null(Formula)){
    Formula <- function(support){
      expr <- (treatment ~ predictor)
      if (sum(support)==0){
        expr <- (treatment ~ 1)
      }
      expr
    }
  }
  fit <- NULL
  dataPredict <- NULL
  dataPredict=data$predictor[sampleSplitIndex,supp]
  if (predictAll){
    dataPredict=data$predictor[,supp]
  }
  if ((method == 'lm')||(method == 'glmnet')){
    if (sum(supp) > 0){
      fit <- glm(Formula(supp), family=binomial, data = dataTrain)
      prediction <- predict(fit, newdata = list(predictor=dataPredict), type="response")
    }
  } else if (method == 'kernel') {
    if (sum(supp) > 0){
      prediction <- ks(dataTrain$predictor, dataTrain$treatment, dataPredict)
      prediction <- (prediction > 0.99) * 0.99 + (prediction < 0.01) * 0.01 + (prediction < 0.99) * (prediction > 0.01) * prediction
    }
  }
  prediction
}

# screening
screening <- function(x, y, method='glmnet', family='Gaussian'){
  var <- apply(x, 2, sd)
  supp <- order(var, decreasing = TRUE)
  if (method=='glmnet'){
    fit <- glmnet::cv.glmnet(x, y, family = family)
    coef <- fit$glmnet.fit$beta[,fit$lambda==fit$lambda.min]
    supp <- (abs(coef)>0)
  } else {
    fit <- VariableScreening::screenIID(x, y, method=method)
    supp <- fit$rank
  }
  supp
}
