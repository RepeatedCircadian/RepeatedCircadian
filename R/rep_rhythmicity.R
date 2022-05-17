library(lme4)
library(MASS)
library(lmtest)
library(pbkrtest)
##' Likelihood-based test for circadian pattern detection with repeated measurement.
##'
##' Test the significance of circadian curve fitting using mixed model with random intercept and likelihood-based test.
##' @title Likelihood-based Test for Detecting Circadian Pattern with Repeated Measurement.
##' @param tt Time vector
##' @param yy Expression vector
##' @param id Subject ID number
##' @param period Period of the since curve. Default is 24.
##' @param method Testing methods can be "LR" or "F". Default is "LR".
##' @return A list of test statistic and pvalue.

##' Formula 1: \eqn{yy = amp * sin(2\pi/period * (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A * sin(2\pi/period * tt) + B * cos(2*\pi/period * tt) + offset}
##' \item{stat}{Test statistic.}
##' \item{pvalue}{P-value from the test.}
##' @author Haocheng Ding, Zhiguang Huo
##' @export
##' @examples
##' Example 1
##'
##' set.seed(32611)
##' m <- 10
##' n <- 12
##' id <- rep(1:n,each=m)
##' rho <- 0.2
##' offset <- runif(1,0,3)
##' sigmaMat <- ifelse(diag(m)==1,1,rho)
##' tt <- rep(runif(m,0,24),n)
##' yy <- as.vector(t(mvrnorm(n,rep(offset,m),sigmaMat)))

##' rpt_rhythmicity(tt, yy, id, period=24, method="LR")
##'
##' Example 2
##' set.seed(32611)
##'m <- 10
##'n <- 12
##'id <- rep(1:n,each=m)
##'rho <- 0.2
##'amp <- 1
##'phase <- 3
##'offset <- runif(1,0,3)
##'sigmaMat <- ifelse(diag(m)==1,1,rho)
##'t <- runif(m,0,24)
##'tt <- rep(t,n)
##'yy <- as.vector(t(mvrnorm(n,amp*sin(2*pi/24*(phase+t))+offset,sigmaMat)))
##'rpt_rhythmicity(tt, yy, id, period=24, method="LR")


rpt_rhythmicity <- function(tt,yy,id,period=24,method="LR"){
  temp_data <- data.frame("y"=yy,"EE"=sin(2*pi/period*tt),"FF"=cos(2*pi/period*tt),"id"=id)

  # lmer with random intercept
  rdm_int_h1 <- lmer(y~EE+FF+(1|id),data=temp_data,REML=F,control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol),check.conv.hess     = .makeCC(action = "ignore", tol = 1e-6),check.conv.grad     = .makeCC("ignore", tol = 2e-3, relTol = NULL)))
  rdm_int_h0 <- lmer(y~(1|id),data=temp_data,REML=F,control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol),check.conv.hess     = .makeCC(action = "ignore", tol = 1e-6),check.conv.grad     = .makeCC("ignore", tol = 2e-3, relTol = NULL)))

  C0 <- coef(summary(rdm_int_h1))[ 1, "Estimate"]
  EE <- coef(summary(rdm_int_h1))[ 2, "Estimate"]
  FF <- coef(summary(rdm_int_h1))[ 3, "Estimate"]
  A <- sqrt(EE^2+FF^2)
  phi <- atan2(FF,EE)/(2*pi/period)
  tempVcov <- as.data.frame(summary(rdm_int_h1)$varcor)
  sigma_alpha <- tempVcov[1,"sdcor"]
  sigma_0 <- tempVcov[2,"sdcor"]


  if(method=="LR"){
    LR_rdm_int <- lrtest(rdm_int_h0,rdm_int_h1)
    stat <- LR_rdm_int[2,4]
    pvalue <- LR_rdm_int[2,5]
  }
  else if(method=="F"){
    stat <- getKR(KRmodcomp(rdm_int_h0,rdm_int_h1),name='Fstat')
    pvalue <- getKR(KRmodcomp(rdm_int_h0,rdm_int_h1),name='p.value')
  }
  else{
    cat("Please check your input! Method only supports 'LR' or 'F'.")}

  return(list(statistic=stat,pvalue=pvalue,A=A,phi=phi,basal=C0,sigma_alpha=sigma_alpha,sigma_0=sigma_0))
}
