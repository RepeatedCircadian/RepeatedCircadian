##' Likelihood-based test for circadian pattern detection with repeated measurement.
##'
##' Test the significance of circadian curve fitting using mixed model with random intercept and likelihood-based test.
##' @title Likelihood-based Test for Detecting Circadian Pattern with Repeated Measurement.
##' @param data A data set with columns: y (gene expression value), t (circadian time), id (subjects IDs) and optional additional covariates.
##' @param type Types of testing model. 'lm' for linear model, 'int' for mixed model with random intercept only and 'slope' for mixed model with random slope. Default is 'int'.
##' @param period Period of the since curve. Default is 24.
##' @param method Testing methods can be "LR" or "F". Default is "LR".
##' @param cov Optional vector of additional covariates adding to model fitting (e.g. cov=c('age','gender')). Default is NULL.
##' @return A list of test statistic, pvalue and parameter estimates.

##' Formula 1: \eqn{yy = amp * sin(2\pi/period * (phase + tt)) + offset}
##' Formula 2: \eqn{yy = A * sin(2\pi/period * tt) + B * cos(2*\pi/period * tt) + offset}
##' \item{stat}{Test statistic.}
##' \item{pvalue}{P-value from the test.}
##' \item{A}{Amplitude estimate.}
##' \item{phi}{Phase estimate.}
##' \item{basal}{Basal level estimate.}
##' \item{sigma_null}{Standard deviation for the null model under linear model setting (type='lm').}
##' \item{sigma_alt}{Standard deviation for the alternative model unde linear model setting (type='lm').}
##' \item{sigma_0}{Standard deviation for the fixed part of intercept under random intercept or random slope setting (type='int' or 'slope').}
##' \item{sigma_alpha}{Standard deviation for the random part of intercept under random intercept or random slope setting (type='int' or 'slope').}
##' \item{coefficients}{Coefficient table of fitted model.}
##' @author Haocheng Ding, Zhiguang Huo
##' @import lme4
##' @import MASS
##' @import lmtest
##' @import pbkrtest
##' @import lmerTest
##' @export
##' @examples
##' Example 1
##'
##' library(MASS)
##' set.seed(32611)
##' m <- 10
##' n <- 12
##' id <- rep(1:n,each=m)
##' rho <- 0.2
##' offset <- runif(1,0,3)
##' sigmaMat <- ifelse(diag(m)==1,1,rho)
##' tt <- rep(runif(m,0,24),n)
##' yy <- as.vector(t(mvrnorm(n,rep(offset,m),sigmaMat)))
##' group <- factor(rep(c('A','B'),each=60),levels = c('A','B'))
##' age <- rep(sample(20:80,n,replace = T),each=m)
##' data <- data.frame(y=yy,t=tt,id=id,group=group,age=age)
##' rpt_rhythmicity_new(data, period=24, type='int', method="LR", cov=NULL)
##'
##' Example 2
##' 
##' library(MASS)
##' set.seed(32611)
##' m <- 10
##' n <- 12
##' id <- rep(1:n,each=m)
##' rho <- 0.2
##' amp <- 1
##' phase <- 3
##' offset <- runif(1,0,3)
##' sigmaMat <- ifelse(diag(m)==1,1,rho)
##' t <- runif(m,0,24)
##' tt <- rep(t,n)
##' yy <- as.vector(t(mvrnorm(n,amp*sin(2*pi/24*(phase+t))+offset,sigmaMat)))
##' group <- factor(rep(c('A','B'),each=60),levels = c('A','B'))
##' age <- rep(sample(20:80,n,replace = T),each=m)
##' data <- data.frame(y=yy,t=tt,id=id,group=group,age=age)
##' rpt_rhythmicity_new(data, type='int', period=24, method="LR",cov=c('group','age'))


rpt_rhythmicity_new <- function(data, period=24, type="int", method="LR", cov=NULL){
  temp_data <- data.frame(data,"EE"=sin(2*pi/period*data$t),"FF"=cos(2*pi/period*data$t))
  
  if(type=="lm"){
    ## linear model without additional covariates.
    if(is.null(cov)){
      model_h0 <- y~1
      model_h1 <- y~EE+FF
    }
    ## linear model under h0 with additional covariates.
    else if(is.null(cov)=="FALSE"){
      model_h0 <- as.formula(paste('y', paste(cov, collapse=" + "), sep=" ~ "))
      model_h1 <- as.formula(paste('y', paste(c('EE','FF',cov), collapse=" + "), sep=" ~ "))
    }
    ## Fitting the linear models
    fitlm_h0 <- lm(model_h0,data=temp_data)
    fitlm_h1 <- lm(model_h1,data=temp_data)
    
    C0 <- coef(summary(fitlm_h1))["(Intercept)", "Estimate"]
    EE <- coef(summary(fitlm_h1))["EE", "Estimate"]
    FF <- coef(summary(fitlm_h1))["FF", "Estimate"]
    A <- sqrt(EE^2+FF^2)
    phi <- atan2(FF,EE)/(2*pi/period)
    sigma_null <- summary(fitlm_h0)$sigma
    sigma_alt <- summary(fitlm_h1)$sigma
    coefficients <- coef(summary(fitlm_h1))
    
    if(method=="LR"){
      LR_stat <- lrtest(fitlm_h0,fitlm_h1)
      stat <- LR_stat[2,4]
      pvalue <- LR_stat[2,5]
    }
    else if(method=="F"){
      F_stat <- anova(fitlm_h0,fitlm_h1)
      stat <- F_stat[2,5]
      pvalue <- F_stat[2,6]
    }
    else{
      cat("Please check your input! Method only supports 'LR' or 'F'.")}
    return(list(statistic=stat,pvalue=pvalue,A=A,phi=phi,basal=C0,sigma_null=sigma_null,sigma_alt=sigma_alt,coefficients=coefficients))
  }
  
  
  
  
  else if(type=="int"){
    ## random intercept linear mixed model without covariates.
    if(is.null(cov)){
      model_h0 <- y~1+(1|id)
      model_h1 <- y~EE+FF+(1|id)
    }
    ## random intercept linear mixed model with covariates.
    else if(!is.null(cov)){
      model_h0 <- as.formula(paste('y', paste(c('(1|id)',cov), collapse=" + "), sep=" ~ "))
      model_h1 <- as.formula(paste('y', paste(c('EE','FF','(1|id)',cov), collapse=" + "), sep=" ~ "))
    }
    
    rdm_int_h0 <- lmer(model_h0,data=temp_data,REML=F,control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol),check.conv.hess     = .makeCC(action = "ignore", tol = 1e-6),check.conv.grad     = .makeCC("ignore", tol = 2e-3, relTol = NULL)))
    rdm_int_h1 <- lmer(model_h1,data=temp_data,REML=F,control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol),check.conv.hess     = .makeCC(action = "ignore", tol = 1e-6),check.conv.grad     = .makeCC("ignore", tol = 2e-3, relTol = NULL)))
    
    C0 <- coef(summary(rdm_int_h1))["(Intercept)", "Estimate"]
    EE <- coef(summary(rdm_int_h1))["EE", "Estimate"]
    FF <- coef(summary(rdm_int_h1))["FF", "Estimate"]
    A <- sqrt(EE^2+FF^2)
    phi <- atan2(FF,EE)/(2*pi/period)
    tempVcov <- as.data.frame(summary(rdm_int_h1)$varcor)
    sigma_alpha <- tempVcov[1,"sdcor"]
    sigma_0 <- tempVcov[nrow(tempVcov),"sdcor"]
    coefficients <- coef(summary(rdm_int_h1))
    
    
    
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
    
    return(list(statistic=stat,pvalue=pvalue,A=A,phi=phi,basal=C0,sigma_alpha=sigma_alpha,sigma_0=sigma_0,coefficients=coefficients))
    
  }
  
  
  
  
  
  else if(type=="slope"){
    ## random slope linear mixed model without covariates.
    if(is.null(cov)){
      model_h0 <- y~1+(1|id)
      model_h1 <- y~EE+FF+(1+EE+FF|id)
    }
    ## random slope linear mixed model with covariates.
    else if(!is.null(cov)){
      model_h0 <- as.formula(paste('y', paste(c('(1|id)',cov), collapse=" + "), sep=" ~ "))
      model_h1 <- as.formula(paste('y', paste(c('EE','FF','(1+EE+FF|id)',cov), collapse=" + "), sep=" ~ "))
    }
    
    rdm_slope_h0 <- lmer(model_h0,data=temp_data,REML=F,control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol),check.conv.hess     = .makeCC(action = "ignore", tol = 1e-6),check.conv.grad     = .makeCC("ignore", tol = 2e-3, relTol = NULL)))
    rdm_slope_h1 <- lmer(model_h1,data=temp_data,REML=F,control = lmerControl(check.conv.singular = .makeCC(action = "ignore", tol = formals(isSingular)$tol),check.conv.hess     = .makeCC(action = "ignore", tol = 1e-6),check.conv.grad     = .makeCC("ignore", tol = 2e-3, relTol = NULL)))
    
    C0 <- coef(summary(rdm_slope_h1))["(Intercept)", "Estimate"]
    EE <- coef(summary(rdm_slope_h1))["EE", "Estimate"]
    FF <- coef(summary(rdm_slope_h1))["FF", "Estimate"]
    A <- sqrt(EE^2+FF^2)
    phi <- atan2(FF,EE)/(2*pi/period)
    tempVcov <- as.data.frame(summary(rdm_slope_h1)$varcor)
    sigma_alpha <- tempVcov[1,"sdcor"]
    sigma_0 <- tempVcov[nrow(tempVcov),"sdcor"]
    coefficients <- coef(summary(rdm_slope_h1))
    
    if(method=="LR"){
      LR_rdm_slope <- lrtest(rdm_slope_h0,rdm_slope_h1)
      stat <- LR_rdm_slope[2,4]
      pvalue <- LR_rdm_slope[2,5]
    }
    else if(method=="F"){
      stat <- getKR(KRmodcomp(rdm_slope_h0,rdm_slope_h1),name='Fstat')
      pvalue <- getKR(KRmodcomp(rdm_slope_h0,rdm_slope_h1),name='p.value')
    }
    else{
      cat("Please check your input! Method only supports 'LR' or 'F'.")}
    
    return(list(statistic=stat,pvalue=pvalue,A=A,phi=phi,basal=C0,sigma_alpha=sigma_alpha,sigma_0=sigma_0,coefficients=coefficients))
    
  }
  
  
  
  else{
    cat("Type only supports 'lm','int' and 'slope'.")
  }
  
    
}
  
  
