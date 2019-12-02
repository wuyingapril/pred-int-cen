## AUC define from probability view

auc.prob.f <- function(lam, kappa, beta, t0){
  quad <- gauss.quad(n = num.node, kind="hermite")
  x0 <- quad$nodes
  x.continuous <- expand.grid(x0, x0)
  x.continuous <- as.matrix(x.continuous)
  
  prob0 <- expand.grid(quad$weights, quad$weights)
  prob.continuous <- prob0[,1]*prob0[,2]
  
  linear.pred.f <- function(x){ sum(x*beta) }
  linear.pred0 <- apply(x.continuous, 1, linear.pred.f)
  surv0 <- exp(-(lam*t0)^kappa*exp(linear.pred0))
  
  corrector.f <- function(x){
    dnorm(x[1], mean = 0, sd = 1)*dnorm(x[2], mean = 0, sd = 1)*exp(x[1]^2)*exp(x[2]^2)
  }
  corrector0 <- apply(x.continuous, 1, corrector.f)
  num0 <- length(corrector0)
  
  surv <- sum(surv0*corrector0*prob.continuous)
  bottom <- surv*(1-surv)
  
  outer.loop.f <- function(jth){
    inner.loop.f <- function(ith){
      if(linear.pred0[ith] == linear.pred0[jth]){
        return(0.5)
      }
      else{
        return(ifelse(linear.pred0[ith] < linear.pred0[jth], 1, 0))
      }
    }
    inner <- sapply(1:num0, inner.loop.f)
    inner <- inner*surv0*corrector0
    inner <- sum(prob.continuous*inner)*(1-surv0[jth])*corrector0[jth]
    return(inner)
  }
  outer <- sapply(1:num0, outer.loop.f)
  top <- sum(outer*prob.continuous)
  return(top/bottom)
}


## true value based on model
roc.f <- function(lam, kappa, beta, t0, c){
  quad <- gauss.quad(n = num.node, kind="hermite")
  x0 <- quad$nodes
  x.continuous <- expand.grid(x0, x0)
  x.continuous <- as.matrix(x.continuous)
  
  prob0 <- expand.grid(quad$weights, quad$weights)
  prob.continuous <- prob0[,1]*prob0[,2]
  
  linear.pred.f <- function(x){ sum(x*beta) }
  linear.pred0 <- apply(x.continuous, 1, linear.pred.f)
  surv0 <- exp(-(lam*t0)^kappa*exp(linear.pred0))
  
  corrector.f <- function(x){
    dnorm(x[1], mean = 0, sd = 1)*dnorm(x[2], mean = 0, sd = 1)*exp(x[1]^2)*exp(x[2]^2)
  }
  corrector0 <- apply(x.continuous, 1, corrector.f)
  num0 <- length(corrector0)
  
  surv <- sum(surv0*corrector0*prob.continuous)
  inner.f <- function(ith){
    Ft0 <- surv0[ith]
    Yhat <- ifelse(Ft0 > c, 1, 0)
    tp <- ifelse(Yhat == 1, 1, 0)*Ft0
    fp <- ifelse(Yhat == 1, 1, 0)*(1-Ft0)
    fn <- ifelse(Yhat == 0, 1, 0)*Ft0
    tn <- ifelse(Yhat == 0, 1, 0)*(1-Ft0)
    return(c(tp, fp, fn, tn))
  }

  terms <- lapply(1:num0, inner.f)
  terms <- do.call("rbind", terms)

  result <- (corrector0*prob.continuous)%*%terms
  tpr <- result[1]/(result[1]+result[3])
  fpr <- result[2]/(result[2]+result[4])
  return(c(tpr, fpr))
}

##############################################################################################
auc.numeric.f <- function(tpr, fpr){
  ans <- 0
  for (i in 2:length(fpr)) {
    ans <- ans + 0.5 * abs(fpr[i] - fpr[i - 1]) * 
      (tpr[i] + tpr[i - 1])
  }
  return(ans)
}

