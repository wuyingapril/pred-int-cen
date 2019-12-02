#################################################################
require(statmod)
require(parallel)
require(glmnet)
require(survival)
require(mvtnorm)
require(MASS)

#################################################################
### Weibull for event time
### Exponential for gap time between assessments

weight.f <- function(t0, Y, x.T, x.A, lam, beta, kappa, alpha0, alpha1, tau){
  if(Y == 0){
    to.be.integrate <- function(t){
      exp(-alpha0*exp(sum(alpha1*x.A))*(t0 - t))*kappa*lam*(lam*t)^(kappa-1)*exp(sum(x.T*beta))*exp(-(lam*t)^kappa*exp(sum(x.T*beta)))
    }
    top <- integrate(to.be.integrate, 0, t0)$value  
    bottom <- 1-exp(-(lam*t0)^kappa*exp(sum(x.T*beta)))
    weight <- 1-top/bottom
  }
  
  if(Y == 1){
    to.be.integrate <- function(t){
      (alpha0*exp(sum(alpha1*x.A)))*exp(-alpha0*exp(sum(alpha1*x.A))*(t-t0))*exp(-(lam*t)^kappa*exp(sum(x.T*beta)))
    }   
    top <- integrate(to.be.integrate, t0, tau)$value 
    bottom <- exp(-(lam*t0)^kappa*exp(sum(x.T*beta)))
    weight <- top/bottom
  }
  return(weight)
}


#################################################################
## PWC Exponential survivor function 
surv.bar.f <- function(ti, zi, lam, beta, cutpoints){
  cstart <- cutpoints$start
  cstop  <- cutpoints$stop
  npiece <- length(cstart)
  ncov <- length(beta)
  nsubj  <- length(ti)
  
  
  beta <- matrix(beta, nrow=ncov, ncol=1)
  zi   <- matrix(zi, nrow=1, ncol=ncov)
  zi.times.beta <- as.vector( zi %*% beta )
  
  Ht <- 0
  for (k in 1:npiece) {
    wk <- 0
    wk <- ifelse( (ti <= cstop[k]) & (ti >= cstart[k]), ti - cstart[k], wk)
    wk <- ifelse( (ti > cstart[k]) & (ti > cstop[k]), rep(cstop[k] - cstart[k],nsubj), wk)
    Ht <- Ht + ( (wk*lam[k])*exp( zi.times.beta ) )
  }
  Fbar <- exp( (-1)*Ht )
  return(Fbar)
}

## PWC Exponential Hazard function 
haz.bar.f <- function(ti, zi, lam, beta, cutpoints) {
  cstart <- cutpoints$start
  cstop  <- cutpoints$stop
  npiece <- length(cstart)
  ncov <- length(beta)
  nsubj  <- length(ti)
  
  
  beta <- matrix(beta, nrow=ncov, ncol=1)
  zi   <- matrix(zi, nrow=1, ncol=ncov)
  zi.times.beta <- as.vector( zi %*% beta )
  
  Ik <- ifelse(ti > cstart, 1, 0)
  Ik <- ifelse(ti <= cstop, Ik, 0)
  ht <- 0
  for (k in 1:npiece) {
    wk <- 0
    wk <- ifelse( (ti <= cstop[k]) & (ti > cstart[k]), 1, 0)
    ht <- ht + wk*lam[k]
  }
  haz <- ht*exp( zi.times.beta ) 
  return(haz)
}

### PWC Exponential for event time
### Exponential for gap time between assessments

weight.pwc.f <- function(t0, Y, x.T, x.A, lam, beta, cutpoints, alpha0, alpha1, tau){
  if(Y == 0){
    to.be.integrate <- function(t){
      exp(-alpha0*exp(sum(alpha1*x.A))*(t0 - t))*surv.bar.f(ti = t, zi = x.T, lam = lam, beta = beta, cutpoints = cutpoints)*haz.bar.f(ti = t, zi = x.T, lam = lam, beta = beta, cutpoints = cutpoints)
    }
    top <- num.integrate(to.be.integrate, 0, t0)  
    bottom <- 1-surv.bar.f(ti = t0, zi = x.T, lam = lam, beta = beta, cutpoints = cutpoints)
    weight <- 1-top/bottom
  }
  
  if(Y == 1){
    to.be.integrate <- function(t){
      (alpha0*exp(sum(alpha1*x.A)))*exp(-alpha0*exp(sum(alpha1*x.A))*(t-t0))*surv.bar.f(ti = t, zi = x.T, lam = lam, beta = beta, cutpoints = cutpoints)
    }   
    top <- num.integrate(to.be.integrate, t0, tau) 
    bottom <- surv.bar.f(ti = t0, zi = x.T, lam = lam, beta = beta, cutpoints = cutpoints)
    weight <- top/bottom
  }
  return(weight)
}


#################################################################
## Semi-parametric model

#################################################################

num.integrate <- function(f, lower, upper, subdivisions = 100, discrete.times = NULL){
  if(is.null(discrete.times)){
    x <- seq(lower, upper, length.out = subdivisions)
  }
  else{
    x <- discrete.times[discrete.times >= lower & discrete.times <= upper]
  }
  d <- x[-1] - x[-length(x)]
  y <- sapply(x, f)
  value <- (y[-1] + y[-length(y)])/2
  return(sum(value*d))
}


###############
haz.cox.f <- function(ti, zi, H0, beta){ 
  H0$dhazard <- H0$hazard - c(0, H0$hazard[-length(H0$hazard)])  
  value <- ifelse(any(H0$time == ti), H0$dhazard[H0$time == ti], 0)
  value <- value*exp(sum(beta*zi))
  return(value)
}

cum.haz.cox.f <- function(ti, zi, H0, beta){ 
  Ht <- ifelse(all(H0$time > ti), 0, H0$hazard[tail(which(H0$time <= ti),1)])
  Ht <- Ht*exp(sum(beta*zi))
  return(Ht)
}


surv.cox.f <- function(ti, zi, H0, beta){ 
  Ht <- ifelse(all(H0$time > ti), 0, H0$hazard[tail(which(H0$time <= ti),1)])
  Ht <- Ht*exp(sum(beta*zi))
  Fbar <- exp( (-1)*Ht )
  return(Fbar)
}




weight.wei.cox.f <- function(t0, Y, x.T, x.A, lam, beta, kappa, H0.A, alpha1, tau){
  if(Y == 0){
    to.be.integrate <- function(t){
      1/surv.cox.f(ti = t, zi = x.A, H0 = H0.A, beta = alpha1)*kappa*lam*(lam*t)^(kappa-1)*exp(sum(x.T*beta))*exp(-(lam*t)^kappa*exp(sum(x.T*beta)))
    }
    top <- num.integrate(to.be.integrate, 10^-5, t0)
    top <- top*surv.cox.f(ti = t0, zi = x.A, H0 = H0.A, beta = alpha1) 
    bottom <- 1-exp(-(lam*t0)^kappa*exp(sum(x.T*beta)))
    weight <- 1-top/bottom
  }
  
  if(Y == 1){
    to.be.integrate <- function(t){
      surv.cox.f(ti = t, zi = x.A, H0 = H0.A, beta = alpha1)*kappa*lam*(lam*t)^(kappa-1)*exp(sum(x.T*beta))*exp(-(lam*t)^kappa*exp(sum(x.T*beta)))
    }   
    top <- num.integrate(to.be.integrate, t0, tau)
    top <- top/surv.cox.f(ti = t0, zi = x.A, H0 = H0.A, beta = alpha1)
    bottom <- exp(-(lam*t0)^kappa*exp(sum(x.T*beta)))
    weight <- 1-top/bottom
  }
  return(weight)
}



weight.pwc.cox.f <- function(t0, Y, x.T, x.A, lam, beta, cutpoints, H0.A, alpha1, tau){
  if(Y == 0){
    to.be.integrate <- function(t){
      1/surv.cox.f(ti = t, zi = x.A, H0 = H0.A, beta = alpha1)*surv.bar.f(ti = t, zi = x.T, lam = lam, beta = beta, cutpoints = cutpoints)*haz.bar.f(ti = t, zi = x.T, lam = lam, beta = beta, cutpoints = cutpoints)
    }
    top <- num.integrate(to.be.integrate, 0, t0)
    top <- top*surv.cox.f(ti = t0, zi = x.A, H0 = H0.A, beta = alpha1) 
    bottom <- 1-surv.bar.f(ti = t0, zi = x.T, lam = lam, beta = beta, cutpoints = cutpoints)
    weight <- 1-top/bottom
  }
  
  if(Y == 1){
    to.be.integrate <- function(t){
      surv.cox.f(ti = t, zi = x.A, H0 = H0.A, beta = alpha1)*surv.bar.f(ti = t, zi = x.T, lam = lam, beta = beta, cutpoints = cutpoints)*haz.bar.f(ti = t, zi = x.T, lam = lam, beta = beta, cutpoints = cutpoints)
    }   
    top <- num.integrate(to.be.integrate, t0, tau)
    top <- top/surv.cox.f(ti = t0, zi = x.A, H0 = H0.A, beta = alpha1)
    bottom <- surv.bar.f(ti = t0, zi = x.T, lam = lam, beta = beta, cutpoints = cutpoints)
  
    weight <- 1-top/bottom
  }
  return(weight)
}






transform.data.f <- function(data){
  data <- data[data$assess.num != 0, ]
  ## add status_i and L_1 and R_1 to the dataset ###
  ## change the dataset into new format with assessment after t1 and nir ###
  nsubj <- length(unique(data$ID))
  outdata <- lapply(1:nsubj, function(ith){
    subdata <- subset(data, ID == unique(data$ID)[ith])
    ti <- subdata$t[1]
    assess.times <- c(0, subdata$assess.stop)
    
    Li <- max(assess.times[assess.times < ti])
    Ri <- ifelse(all(assess.times < ti), 9999, min(assess.times[assess.times > ti]))
    outdata0 <- data.frame(ID = ith, t = ti, L = Li, R = Ri)
    return(outdata0)
  })
  outdata <- do.call("rbind", outdata)
  outdata <- data.frame(outdata)
  covdata <- unique(data[,c("ID", "x", "z", "v")])
  outdata <- merge(outdata, covdata, by = "ID")
  
  return(outdata)
}







