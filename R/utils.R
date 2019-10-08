
# loss gets the loss of the chosen loss
loss <- function(x, loss_type){
  switch(loss_type,
         logistic = log(1+exp(-x)),
         exponential = exp(-x),
         zeroOne = (x<0)
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

}

# getOutcomeModel ontains the outcome regression model
getOutcomeModel <- function(data, method=c('lm', 'glmnet', 'kernel', 'others'), sampleSplitIndex, formula = NULL){
  p <- dim(data$predictor)[2]
  size <- dim(data$predictor)[1]
  fit <- NULL
  supp <- NULL
  prediction <- array(0,c(sum(sampleSplitIndex),1))
  dataControl <- data.frame(predictor=data$predictor[(!sampleSplitIndex) && (data$treatment==FALSE),], outcome=data$outcome[(!sampleSplitIndex) && (data$treatment==FALSE)])
  dataTreatment <- data.frame(predictor=data$predictor[(!sampleSplitIndex) && (data$treatment==TRUE),], outcome=data$outcome[(!sampleSplitIndex) && (data$treatment==TRUE)])
  if (0.5*size >= p){
    supp$control <- supp$treatment <- rep(TRUE, times = p)
  }
  if ((0.5*size < p) || (method == 'glmnet')) {
    fit$control <- cv.glmnet(x = predictor, y = outcome, data = dataControl)
    fit$treatment <- cv.glmnet(x = predictor, y = outcome, data = dataTreatment)
    prediction[sampleSplitIndex && (data$treatment==FALSE)] <- predict(fit$control, newx = data.frame(predictor=data$predictor[sampleSplitIndex && (data$treatment==FALSE),]), s=fit$control$lambda.min)
    prediction[sampleSplitIndex && (data$treatment==TRUE)] <- predict(fit$treatment, newx = data.frame(predictor=data$predictor[sampleSplitIndex && (data$treatment==TRUE),]), s=fit$treatment$lambda.min)
    supp$control <- fit$control$glmnet.fit$beta[,fit$control$glmnet.fit$lambda==fit$control$lambda.1se]
    supp$treatment <- fit$treatment$glmnet.fit$beta[,fit$treatment$glmnet.fit$lambda==fit$treatment$lambda.1se]
    dataControl$predictor <- data$predictor[,supp$control]
    dataTreatment$predictor <- data$predictor[,supp$treatment]
  }

  if ((method == 'lm')||(method == 'glmnet')){
    fit$control <- lm(formula, data = dataControl)
    fit$treatment <- lm(formula, data = dataTreatment)
    prediction[sampleSplitIndex && (data$treatment==FALSE)] <- predict(fit$control, newdata = data.frame(predictor=data$predictor[sampleSplitIndex && (data$treatment==FALSE),supp$control]))
    prediction[sampleSplitIndex && (data$treatment==TRUE)] <- predict(fit$treatment, newdata = data.frame(predictor=data$predictor[sampleSplitIndex && (data$treatment==TRUE),supp$treatment]))
  } else if (method == 'kernel') {
    prediction[sampleSplitIndex && (data$treatment==FALSE)] <- ks(dataControl$predictor, dataControl$outcome, data$predictor[sampleSplitIndex && (data$treatment==FALSE),supp$control])
    prediction[sampleSplitIndex && (data$treatment==TRUE)] <- ks(dataTreatment$predictor, dataTreatment$outcome, data$predictor[sampleSplitIndex && (data$treatment==TRUE),supp$control])
  } else {
    fit$control <- eval(parse(text=model.gam(dataControl)))
    fit$treatment <- eval(parse(text=model.gam(dataTreatment)))
  }
}

ITRFit <- function(predictor, treatment, outcome, loss = c('logistic')){

}
