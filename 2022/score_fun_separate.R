# Models included:
# LS
# LogF(logistf)
# GLS
# EGLS
# PhyLog(or phyloglm)
# PLogML_LR_a (based on phyloglm)
# PLog(based on phyloglm_aaf())
# phyloglmA(same as PLogF, based on phyloglm_aff())
# PLogML_all(based on phyloglm_ML(), includes PLogML, PLogML_LR and PLogML_LR_p)
library(phylolm)
LS <- function(xx,y){
  XX <- t(xx) %*% xx
  XY <- t(xx) %*% y
  b <- solve(XX, XY)
  h <- y - (xx %*% b)
  MSE <- t(h) %*% h/(p - 2)
  iXX <- solve(XX)
  bSE <- (MSE * iXX[2, 2])^0.5
  scoreLS <- b[2]/bSE
  return(scoreLS)
}

LogF <- function(X,y){
  # same as Logistf
  fit_Log <- logistf(y ~ X)
  scoreLog <- sign(fit_Log$coef[2]) * qchisq(fit_Log$prob[2], df=1, lower.tail=F)
  return(scoreLog)
}
GLS <- function(xx,y,iC){
  XiCX <- t(xx) %*% iC %*% xx
  XiCY <- t(xx) %*% iC %*% y
  b <- solve(XiCX, XiCY)
  h <- y - (xx %*% b)
  MSE <- t(h) %*% iC %*% h/(p - 2)
  iXiCX <- solve(XiCX)
  bSE <- (MSE * iXiCX[2, 2])^0.5
  scoreGLS <- b[2]/bSE
  return(scoreGLS)
}

EGLS <- function(xx,y,C) {
  #estimated generalized least squares
  egls <- function(par = par, xx, y, C) {
    p <- dim(xx)[1]
    CC <- (1 - par) * diag(p) + par * C
    CC <- CC/det(CC)^(1/p) 
    iCC <- solve(CC)
    XiCCX <- t(xx) %*% iCC %*% xx
    XiCCY <- t(xx) %*% iCC %*% y
    b <- solve(XiCCX, XiCCY)
    h <- y - (xx %*% b)
    MSE <- t(h) %*% iCC %*% h/(p - 2)
    return(MSE)
  }
  lambda <- optim(fn=egls, par=par, xx=xx, y=y, C=C, method = "Brent", upper = 1, lower = 0)
  
  CC <- (1 - lambda$par) * diag(dim(xx)[1]) + lambda$par * C
  iCC <- solve(CC)
  XiCCX <- t(xx) %*% iCC %*% xx
  XiCCY <- t(xx) %*% iCC %*% y
  b <- solve(XiCCX, XiCCY)
  h <- y - (xx %*% b)
  MSE <- t(h) %*% iCC %*% h/(p - 2)
  iXiCCX <- solve(XiCCX)
  bSE <- (MSE * iXiCCX[2, 2])^0.5
  scoreEGLS <- b[2]/bSE
  return(scoreEGLS)
}

PhyLog <- function(X,y,phy.extend){
  # in the ms it was called phyloglm but that would be confused with phylolm::phyloglm
  fit <- phylolm::phyloglm(y ~ X, phy=phy.extend,method = "logistic_MPLE")
  # using b1/se{b1} as the scoring statistic
  scorePhyLog <- fit$coef[2]/fit$sd[2] 
  return(scorePhyLog)
}

PLogML_LR_a <- function(X,y,phy.extend){
  # PLogML_LR_a is the same as PLogML_LR_p except that a is using phyloglm and p is using phyloglm_ML. 
  # Also mentioned as PLogML_a in SimulationScoring_EGLS_Tony_15Dec16.R.
  if(any(X == y) & any(X != y)) {
    # The default method phyloglm uses is logistic_MPLE
    fit_ML_a <- try(phyloglm(y ~ X, phy=phy.extend), silent=T)
    fit_ML0_a <- try(phyloglm(y ~ 1, phy=phy.extend), silent=T)
    if(is.null(attr(fit_ML_a,'condition')) & is.null(attr(fit_ML0_a,'condition'))){
      LR_a <- fit_ML_a$penlogLik - fit_ML0_a$penlogLik
      if(LR_a > 0){
        scorePLogML_LR_a <- sign(fit_ML_a$coef[2]) * 2 * LR_a
      }else{
        scorePLogML_LR_a <- 0
      }
      # output$convergencePLogML[i] <- fit_ML$convergence
      # output$convergencePLogML0[i] <- fit_ML0$convergence
    }else{
      scorePLogML_LR_a <- 0
      # output$convergencePLogML[i] <- 100
      # output$convergencePLogML0[i] <- 100
    }
  }else{
    if(any(X == y)){
      scorePLogML_LR_a <- Inf
    }else{
      scorePLogML_LR_a <- -Inf
    }
  }   
  return(scorePLogML_LR_a)
}

PLog_old <- function(X,y,phy.extend){
  fit <- phyloglm_aaf(y ~ X, phy=phy.extend, Firth=F)
  scorePLog <- fit$zB[2]
  return(scorePLog)
}

PLog <- function(X,y,phy.extend){
  # trying to catch some warnings
  if(any(X == y) & any(X != y)) {
    fit_F <- try(phyloglm_aaf(y ~ X, phy=phy.extend, Firth=F))
    if(is.null(attr(fit_F,'condition'))){
      scorePLog <- fit_F$zB[2]
      # output$convergencePLogF[i]<- fit_F$convergence
    }else{
      scorePLog <- 0
      # output$convergencePLogF[i] <- 100
    }
  }else{
    if(any(X == y)){
      scorePLog <- Inf
    }else{
      scorePLog <- -Inf
    }
  }
  return(scorePLog)
}

phyloglmA <- function(X,y,phy.extend){
  # same as PLogF
  # depends on phyloglm_aaf
  if(any(X == y) & any(X != y)) {
    fit_F <- try(phyloglm_aaf(y ~ X, phy=phy.extend, Firth=T, ultrametric=F), silent=T)
    if(is.null(attr(fit_F,'condition'))){
      scorePLogF <- fit_F$zB[2]
      # output$convergencePLogF[i]<- fit_F$convergence
    }else{
      scorePLogF <- 0
      # output$convergencePLogF[i] <- 100
    }
  }else{
    if(any(X == y)){
      scorePLogF <- Inf
    }else{
      scorePLogF <- -Inf
    }
  }
  return(scorePLogF)
}

PLogML_all <- function(X,y,phy.extend){
  # PLogML_all includes PLogML, PLogML_LR and PLogML_LR_p
  # PLogML_all depends on phyloglm_ML
  if(any(X == y) & any(X != y)) {
    fit_ML <- try(phyloglm_ML(y ~ X, phy=phy.extend), silent=T)
    fit_ML0 <- try(phyloglm_ML(y ~ 1, phy=phy.extend), silent=T)
    if(is.null(attr(fit_ML,'condition')) & is.null(attr(fit_ML0,'condition'))){
      scorePLogML <- fit_ML$coef[2]/fit_ML$sd[2]
      LR <- fit_ML$logLik - fit_ML0$logLik
      LR_p <- fit_ML$penlogLik - fit_ML0$penlogLik
      if(LR > 0){
        scorePLogML_LR <- sign(fit_ML$coef[2]) * 2 * LR
      }else{
        scorePLogML_LR <- 0
      }
      if(LR_p > 0){
        scorePLogML_LR_p <- sign(fit_ML$coef[2]) * 2 * LR_p
      }else{
        scorePLogML_LR_p <- 0
      }
    }else{
      scorePLogML <- 0
      scorePLogML_LR <- 0
      scorePLogML_LR_p <- 0
    }
  }else{
    if(any(X == y)){
      scorePLogML <- Inf
      scorePLogML_LR <- Inf
      scorePLogML_LR_p <- Inf
    }else{
      scorePLogML <- Inf
      scorePLogML_LR <- -Inf
      scorePLogML_LR_p <- -Inf
    }
  }
  #R functions don't return multiple objects in the strict sense.
    return(c(scorePLogML,scorePLogML_LR,scorePLogML_LR_p))
}

