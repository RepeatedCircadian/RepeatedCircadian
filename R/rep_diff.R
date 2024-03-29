##' Likelihood-based test for differential circadian pattern detection with repeated measurement.
##'
##' Test differential rhythmicity of circadian curve fitting using random intercept linear mixed model and likelihood-based tests.
##' @title Likelihood-based Tests for Detecting Differential Circadian Pattern with Repeated Measurement
##' @param tt1 Time vector of condition 1.
##' @param yy1 Expression vector of condition 1.
##' @param tt2 Time vector of condition 2.
##' @param yy2 Expression vector of condition 2.
##' @param id1 Subject ID numbers in group 1.
##' @param id2 Subject ID numbers in group 2.
##' @param group Group indicator of all subjects.
##' @param period Period of the since curve. Default is 24.
##' @param method Test used to detect differential circadian pattern. It can be chosen either "LR" or "F". Default is LR.

##' @return A list, see details below.
##'

##' Formula 1: \eqn{yy = amp * sin(2\pi/period * (phase + tt)) + offset.}
##' Formula 2: \eqn{yy = A * sin(2\pi/period * tt) + B * cos(2*\pi/period * tt) + offset.}

##' \item{stat}{Test statistic.}
##' \item{pvalue}{P-value from the test.}
##' \item{A1}{Amplitude for group 1.}
##' \item{A2}{Amplitude for group 2.}
##' \item{phi1}{Phase for group 1.}
##' \item{phi2}{Phase for group 2.}
##' \item{basal1}{Basal level for group 1.}
##' \item{basal2}{Basal level for group 2.}
##' \item{sigma_0}{Standard deviation for the fixed part of intercept.}
##' \item{sigma_alpha}{Standard deviation for the random part of intercept.}
##' @author Haocheng Ding, Zhiguang Huo
##' @import lme4
##' @import MASS
##' @import lmtest
##' @import pbkrtest
##' @export
##' @examples
##' Example 1
##'
##'set.seed(32611)
##'n1 <- 12
##'n2 <- 12
##'m1 <- 10
##'m2 <- 10
##'id1 <- rep(1:n1,each=m1)
##'id2 <- rep((1+n1):(n1+n2),each=m2)
##'rho <- 0.2
##'offset1 <- runif(1,0,3)
##'offset2 <- runif(1,2,5)
##'amp1 <- 1
##'amp2 <- 2
##'phase1 <- 3
##'phase2 <- 6
##'sigmaMat1 <- ifelse(diag(m1)==1,1,rho)
##'sigmaMat2 <- ifelse(diag(m2)==1,1,rho)
##'t1 <- runif(m1,0,24)
##'tt1 <- rep(t1,n1)
##'yy1 <- as.vector(t(mvrnorm(n1,amp1*sin(2*pi/24*(phase1+t1))+offset1,sigmaMat1)))

##'t2 <- runif(m2,0,24)
##'tt2 <- rep(t2,n2)
##'yy2 <- as.vector(t(mvrnorm(n2,amp2*sin(2*pi/24*(phase2+t2))+offset2,sigmaMat2)))

##'group <- c(rep("group1",n1*m1),rep("group2",n2*m2))
##'rpt_diff(tt1,yy1,tt2,yy2,id1,id2,group)
##'
##'
##' Example 2
##'
##'set.seed(32611)
##'n1 <- 12
##'n2 <- 12
##'m1 <- 10
##'m2 <- 10
##'id1 <- rep(1:n1,each=m1)
##'id2 <- rep((1+n1):(n1+n2),each=m2)
##'rho <- 0.2
##'offset1 <- offset2 <- runif(1,0,3)
##'sigmaMat1 <- ifelse(diag(m1)==1,1,rho)
##'sigmaMat2 <- ifelse(diag(m2)==1,1,rho)
##'t1 <- runif(m1,0,24)
##'tt1 <- rep(t1,n1)
##'yy1 <- as.vector(t(mvrnorm(n1,rep(offset1,m1),sigmaMat1)))

##'t2 <- runif(m2,0,24)
##'tt2 <- rep(t2,n2)
##'yy2 <- as.vector(t(mvrnorm(n2,rep(offset2,m2),sigmaMat2)))

##'group <- c(rep("group1",n1*m1),rep("group2",n2*m2))
##'rpt_diff(tt1,yy1,tt2,yy2,id1,id2,group)

rpt_diff <- function(tt1,yy1,tt2,yy2,id1,id2,group,period=24,method="LR"){
  EE <- c(sin(2*pi/period*tt1),sin(2*pi/period*tt2))
  FF <- c(cos(2*pi/period*tt1),cos(2*pi/period*tt2))
  id <- c(id1,id2)
  Y <- c((yy1),(yy2))
  testData <- data.frame(EE,FF,group,Y,id)

  rdm_int_h0 <- lmer(Y~EE+FF+group+(1|id),data=testData,REML=F,control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol),check.conv.hess     = .makeCC(action = "ignore", tol = 1e-6),check.conv.grad     = .makeCC("ignore", tol = 2e-3, relTol = NULL)))

  rdm_int_h1 <- lmer(Y~EE+FF+group+EE*group+FF*group+(1|id),data=testData,REML=F,control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol),check.conv.hess     = .makeCC(action = "ignore", tol = 1e-6),check.conv.grad     = .makeCC("ignore", tol = 2e-3, relTol = NULL)))


  summh1 <- summary(rdm_int_h1)
  C01 <- coef(summh1)[ 1, "Estimate"]
  C02 <- coef(summh1)[ 1, "Estimate"]+coef(summh1)[4, "Estimate"]
  EE1 <- coef(summh1)[ 2, "Estimate"]
  FF1 <- coef(summh1)[ 3, "Estimate"]
  EE2 <- coef(summh1)[ 2, "Estimate"]+coef(summh1)[5, "Estimate"]
  FF2 <- coef(summh1)[ 3, "Estimate"]+coef(summh1)[6, "Estimate"]

  A1 <- sqrt(EE1^2+FF1^2)
  phi1 <- atan2(FF1,EE1)/(2*pi/period)

  A2 <- sqrt(EE2^2+FF2^2)
  phi2 <- atan2(FF2,EE2)/(2*pi/period)
  tempVcov2 <- as.data.frame(summh1$varcor)
  sigma_alpha <- tempVcov2[1,"sdcor"]
  sigma_0 <- tempVcov2[2,"sdcor"]



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
    cat("Please check your input! Method only supports 'LR' or 'F'.")
  }
  return(list(statistic=stat,pvalue=pvalue,A1=A1,A2=A2,phi1=phi1,phi2=phi2,basal1=C01,basal2=C02,sigma_alpha=sigma_alpha,sigma_0=sigma_0))
}
