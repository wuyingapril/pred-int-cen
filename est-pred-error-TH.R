source('weight.f.all.r')

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
### Exponential for gap time between assessments
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
  gap <- train.data$assess.stop - train.data$assess.start
  status <- train.data$assess.status
  x <- train.data[, c("x", "v")]
  x <- as.matrix(x)
  
  fit <- survreg(Surv(gap, status) ~ x, dist="exponential")
  est.alpha0 <- exp(-fit$coef[1])
  est.alpha1 <- -fit$coef[-1]

###############################################################################
## estimate the prediction error
  est.pred.error.f <- function(outdata, t0, lam, beta, kappa, alpha0, alpha1, tau){  
    m <- nrow(outdata)
    Delta <- ifelse(outdata$L > t0, 1, 0)
    Delta <- ifelse(outdata$R < t0, 1, Delta)
    
    Y <- rep(NA, nrow(outdata))
    Y <- ifelse(outdata$L > t0, 1, Y)
    Y <- ifelse(outdata$R < t0, 0, Y)
    
    Yhat <- rep(NA, m)
    weight <- rep(NA, m)
    pe <- rep(NA, m)
    pe.impute <- rep(NA, m)
    for(ith in 1:m){
      x.T <- outdata[ith, c("x", "z")]
      x.A <- outdata[ith, c("x", "v")]
      
      Ft0 <- exp(-(lam*t0)^kappa*exp(sum(x.T*beta)))
      FLi <- exp(-(lam*outdata$L[ith])^kappa*exp(sum(x.T*beta)))
      FRi <- exp(-(lam*outdata$R[ith])^kappa*exp(sum(x.T*beta)))
      
      Yhat[ith] <- ifelse(Ft0 > 0.5, 1, 0)
      pe[ith] <- (1-Yhat[ith])^2*Ft0 + (0-Yhat[ith])^2*(1-Ft0)
      
      if(Delta[ith] == 1){
        weight[ith] <- weight.f(t0 = t0, Y = Y[ith], x.T = x.T, x.A = x.A, 
                                lam = lam, beta = beta, 
                                kappa = kappa, 
                                alpha0 = alpha0, 
                                alpha1 = alpha1, tau = tau)
        
      }
      if(Delta[ith]==0){
        pe.impute[ith] <- ((1-Yhat[ith])^2*(Ft0-FRi) + (0-Yhat[ith])^2*(FLi-Ft0))/(FLi- FRi)
      }
      
    }
    pe.unweighted <- mean(((Y-Yhat)^2)[Delta == 1])
    pe.ipw <- sum(((Y-Yhat)^2/weight)[Delta == 1])/m
    pe.aipw <- pe.ipw + mean(pe) - sum((pe/weight)[Delta == 1])/m
    pe.imp <- sum(((Y-Yhat)^2)[Delta == 1] , (pe.impute)[Delta == 0])/m
    result <- list(Unweighted = pe.unweighted, IPW = pe.ipw, 
                   AIPW = pe.aipw, IMP = pe.imp)
    return(result)
  }
  
  t0 <- 1
  pe.est <- est.pred.error.f(test, t0 = t0,
                          lam = est.lam, beta = est.beta, kappa = est.kappa,
                          alpha0 = est.alpha0, alpha1 = est.alpha1, tau = tau)
  
  pe.est
  