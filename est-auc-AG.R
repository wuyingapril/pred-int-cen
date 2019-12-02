source('weight.f.all.r')
source("auc.f.r")

train.data <- read.csv('trainsample.csv', header = T)
train <- transform.data.f(data = train.data)
head(train)
test.data <- read.csv('testsample.csv', header = T)
test <- transform.data.f(data = test.data)
head(test)

## administrative censoring time
tau <- 4
###############################################################################
## Weibull for event time
### Anderson-Gill model for gap time between assessments
###############################################################################

## Event Process
## From Training data
  x <- train[, c("x", "z")]
  x <- as.matrix(x)
  L <- train$L
  R <- train$R
  status <- ifelse(R == 9999, 0, 3)
  status <- ifelse(L == 0, 2, status)
  
  L <- ifelse(L == 0, NA, L)
  R <- ifelse(R == 9999, NA, R)
  
  fit <- survreg(Surv(L, R, type = "interval2") ~ x, dist="weibull")
  est.lam <- exp(-fit$coef[1])
  est.beta <- -fit$coef[-1]/fit$scale
  est.kappa <- 1/fit$scale
##############################################################################
## Inspection 
  x <- train.data[, c("x", "v")]
  x <- as.matrix(x)
  
  fit <- coxph(Surv(train.data$assess.start, train.data$assess.stop, train.data$assess.status) ~ x, 
               method = 'breslow')
  est.H0.AG <- basehaz(fit, center = FALSE)
  est.alpha1.AG <- fit$coef

###############################################################################
## estimate the area under ROC curve
  est.auc.f <- function(outdata, t0, lam, beta, kappa, H0.A, alpha1, tau){  
    m <- nrow(outdata)
    Delta <- ifelse(outdata$L > t0, 1, 0)
    Delta <- ifelse(outdata$R < t0, 1, Delta)
    
    Y <- rep(NA, nrow(outdata))
    Y <- ifelse(outdata$L > t0, 1, Y)
    Y <- ifelse(outdata$R < t0, 0, Y)
    
    linear.pred <- rep(NA, m)
    Ft0 <- rep(NA, m)
    FLi <- rep(NA, m)
    FRi <- rep(NA, m)
    
    Yhat <- rep(NA, m)
    weight <- rep(NA, m)
    
    for(ith in 1:m){
      x.T <- outdata[ith, c("x", "z")]
      x.A <- outdata[ith, c("x", "v")]
      linear.pred[ith] <- sum(x.T*beta)
      
      Ft0[ith] <- exp(-(lam*t0)^kappa*exp(sum(x.T*beta)))
      FLi[ith] <- exp(-(lam*outdata$L[ith])^kappa*exp(sum(x.T*beta)))
      FRi[ith] <- exp(-(lam*outdata$R[ith])^kappa*exp(sum(x.T*beta)))
      
      if(Delta[ith] == 1){
        weight[ith] <- weight.wei.cox.f(t0 = t0, Y = Y[ith], x.T = x.T, x.A = x.A, 
                                        lam = lam, beta = beta, 
                                        kappa = kappa, H0.A = H0.A, 
                                        alpha1 = alpha1, tau = tau)
        
      }
    }
    
    weight[is.na(weight)] <- 1
    
    outer.loop.f <- function(jth){
      inner.loop.f <- function(ith){
        if(linear.pred[ith] == linear.pred[jth]){
          return(0.5)
        }
        else{
          return(ifelse(linear.pred[ith] < linear.pred[jth], 1, 0))
        }
      }
      return(sapply(1:m, inner.loop.f))
    }
    terms <- lapply(1:m, outer.loop.f)
    terms <- do.call("cbind", terms)
    # sum(sapply(1:m, function(ith){sum(na.omit(terms[ith, ]))}))
    ## Unweighted
    bottom.matrix <- (Delta*Y)%*%t(Delta*(1-Y))
    bottom <- sum(sapply(1:m, function(ith){sum(na.omit(bottom.matrix[ith, ]))}))
    
    top.matrix <- terms*bottom.matrix
    top <- sum(sapply(1:m, function(ith){sum(na.omit(top.matrix[ith, ]))}))
    auc.unweighted <- top/bottom
    
    ## IPW 
    bottom.matrix <- (Delta/weight*Y)%*%t(Delta/weight*(1-Y))
    bottom <- sum(sapply(1:m, function(ith){sum(na.omit(bottom.matrix[ith, ]))}))
    
    top.matrix <- terms*bottom.matrix
    top <- sum(sapply(1:m, function(ith){sum(na.omit(top.matrix[ith, ]))}))
    auc.ipw <- top/bottom
    
    ## AIPW 
    bottom.matrix.aug1 <- (Delta/weight*Y)%*%t((1-Delta/weight)*(1-Ft0))
    agb1 <- sum(sapply(1:m, function(ith){sum(na.omit(bottom.matrix.aug1[ith, ]))}))
    bottom.matrix.aug2 <- ((1-Delta/weight)*Ft0)%*%t((Delta/weight)*(1-Y))
    agb2 <- sum(sapply(1:m, function(ith){sum(na.omit(bottom.matrix.aug2[ith, ]))}))
    aug.term <- Ft0%*%t(1-Ft0)
    bottom.matrix.aug3 <- ((1-Delta/weight)%*%t(1-Delta/weight))*aug.term
    agb3 <- sum(bottom.matrix.aug3)
    
    top.matrix.aug1 <- terms*bottom.matrix.aug1
    agt1 <- sum(sapply(1:m, function(ith){sum(na.omit(top.matrix.aug1[ith, ]))}))
    top.matrix.aug2 <- terms*bottom.matrix.aug2
    agt2 <- sum(sapply(1:m, function(ith){sum(na.omit(top.matrix.aug2[ith, ]))}))
    top.matrix.aug3 <- terms*bottom.matrix.aug3
    agt3 <- sum(sapply(1:m, function(ith){sum(na.omit(top.matrix.aug3[ith, ]))}))
    
    bottom <- bottom + agb1 + agb2 + agb3
    top <- top + agt1 + agt2 + agt3
    auc.aipw <- top/bottom
    
    ## Imputation-based
    Ynew <- ifelse(is.na(Y),999,Y)
    bottom.matrix <- (Delta*Ynew + (1-Delta)*(Ft0-FRi)/(FLi- FRi))%*%
      t(Delta*(1-Ynew)+ (1-Delta)*(FLi-Ft0)/(FLi- FRi))
    bottom <- sum(sapply(1:m, function(ith){sum(na.omit(bottom.matrix[ith, ]))}))
    
    top.matrix <- terms*bottom.matrix
    top <- sum(sapply(1:m, function(ith){sum(na.omit(top.matrix[ith, ]))}))
    auc.imp <- top/bottom
    
    
    result <- list(Unweighted = auc.unweighted, IPW = auc.ipw, 
                   AIPW = auc.aipw, IMP = auc.imp)
    return(result)
  }
  
  
  
  t0 <- 1
  auc.est <- est.auc.f(test, t0 = t0,
                             lam = est.lam, beta = est.beta, kappa = est.kappa,
                             H0.A = est.H0.AG, alpha1 = est.alpha1.AG, tau = tau)
  
  auc.est
  