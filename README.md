# RepeatedCircadian
Repeated measured omics data circadian and differential circadian analysis

## Install This Package from github
* In R console

```{R}
library(devtools)
install_github("RepeatedCircadian/RepeatedCircadian") 
```

## Citation


## Full tutorial
http://htmlpreview.github.io/?https://github.com/RepeatedCircadian/RepeatedCircadian/blob/main/vignettes/RepeatedCircadian_tutorial.html

## Short tutorial for circadian pattern detection with repeated measured data

```{R}
library(RepeatedCircadian)
require(MASS)

#Example 1

set.seed(32611)
m <- 10
n <- 12
id <- rep(1:n,each=m)
rho <- 0.2
offset <- runif(1,0,3)
sigmaMat <- ifelse(diag(m)==1,1,rho)
tt <- rep(runif(m,0,24),n)
yy <- as.vector(t(mvrnorm(n,rep(offset,m),sigmaMat)))
rpt_rhythmicity(tt, yy, id, period=24, method="LR")

#Example 2

set.seed(32611)
m <- 10
n <- 12
id <- rep(1:n,each=m)
rho <- 0.2
amp <- 1
phase <- 3
offset <- runif(1,0,3)
sigmaMat <- ifelse(diag(m)==1,1,rho)
t <- runif(m,0,24)
tt <- rep(t,n)
yy <- as.vector(t(mvrnorm(n,amp*sin(2*pi/24*(phase+t))+offset,sigmaMat)))
rpt_rhythmicity(tt, yy, id, period=24, method="LR")

```

## Short tutorial for differential circadian analysis with repeated measured data

```{R}
library(RepeatedCircadian)

#Example 1

set.seed(32611)
n1 <- 12
n2 <- 12
m1 <- 10
m2 <- 10
id1 <- rep(1:n1,each=m1)
id2 <- rep((1+n1):(n1+n2),each=m2)
rho <- 0.2
offset1 <- runif(1,0,3)
offset2 <- runif(1,2,5)
amp1 <- 1
amp2 <- 2
phase1 <- 3
phase2 <- 6
sigmaMat1 <- ifelse(diag(m1)==1,1,rho)
sigmaMat2 <- ifelse(diag(m2)==1,1,rho)
t1 <- runif(m1,0,24)
tt1 <- rep(t1,n1)
yy1 <- as.vector(t(mvrnorm(n1,amp1*sin(2*pi/24*(phase1+t1))+offset1,sigmaMat1)))
t2 <- runif(m2,0,24)
tt2 <- rep(t2,n2)
yy2 <- as.vector(t(mvrnorm(n2,amp2*sin(2*pi/24*(phase2+t2))+offset2,sigmaMat2)))
group <- c(rep("group1",n1*m1),rep("group2",n2*m2))
rpt_diff(tt1,yy1,tt2,yy2,id1,id2,group)


#Example 2

set.seed(32611)
n1 <- 12
n2 <- 12
m1 <- 10
m2 <- 10
id1 <- rep(1:n1,each=m1)
id2 <- rep((1+n1):(n1+n2),each=m2)
rho <- 0.2
offset1 <- offset2 <- runif(1,0,3)
sigmaMat1 <- ifelse(diag(m1)==1,1,rho)
sigmaMat2 <- ifelse(diag(m2)==1,1,rho)
t1 <- runif(m1,0,24)
tt1 <- rep(t1,n1)
yy1 <- as.vector(t(mvrnorm(n1,rep(offset1,m1),sigmaMat1)))
t2 <- runif(m2,0,24)
tt2 <- rep(t2,n2)
yy2 <- as.vector(t(mvrnorm(n2,rep(offset2,m2),sigmaMat2)))
group <- c(rep("group1",n1*m1),rep("group2",n2*m2))
rpt_diff(tt1,yy1,tt2,yy2,id1,id2,group)
```



