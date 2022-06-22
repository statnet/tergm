#  File tests/dynamic_MLE_blockdiag.bipartite.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2022 Statnet Commons
################################################################################
library(statnet.common)
opttest({
library(tergm)
options(tergm.eval.loglik=FALSE)

tolerance<-0.05
n<-20
m<-13
theta<--.4

z.error <- function(truth, est, variance){
  if(truth==est) 0 # Infinite case
  else abs(truth-est)/sqrt(variance)
}

prop.weights <- c("default", "random")

logit<-function(p) log(p/(1-p))

block.dyadcount<-function(y, na.omit=TRUE){
  a <- y %v% "a"
  M <- outer(a,a,"==")
  M[1:m,1:m]<-0
  M[(m+1):n,(m+1):n]<-0
  M[lower.tri(M, TRUE)]<-0
  if(na.omit) M[as.edgelist(is.na(y))] <- 0
  sum(M)
}

form.mle<-function(y0,y1){
  logit(network.edgecount(y1-y0,na.omit=TRUE)/(block.dyadcount(y1)-network.edgecount(y0-is.na(y1))))
}

diss.mle<-function(y0,y1){
  -logit(network.edgecount(y0-y1,na.omit=TRUE)/(network.edgecount(y0-is.na(y1))))
}

y0 <- network.initialize(n, directed=FALSE, bipartite=m)
a <- rep(1:20,1:20)[1:n]
a <- unlist(split(a, rep(1:2, n/2)))
a <- c(sort(a[1:m]), sort(a[-(1:m)]))
y0 %v% "a" <- a

set.seed(321)
y0<-simulate(y0~edges, constraints=~blockdiag("a"), coef=theta, control=control.simulate(MCMC.burnin=n^2*2), dynamic=FALSE)

cat("Complete data:\n")

set.seed(125)
y1<-simulate(y0~edges, constraints=~blockdiag("a"), coef=theta, control=control.simulate(MCMC.burnin=n^2*2), dynamic=FALSE)

# Force CMPLE
set.seed(543)
fit<-tergm(list(y0,y1) ~ Form(~edges) + Persist(~edges), constraints=~blockdiag("a"), estimate="CMPLE", times=c(1,2))

stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(form.mle(y0,y1), coef(fit)[1], vcov(fit)[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1), coef(fit)[2], vcov(fit)[2,2]) <= tolerance)

# Autodetected CMPLE
set.seed(543)
fit<-tergm(list(y0,y1) ~ Form(~edges) + Persist(~edges), constraints=~blockdiag("a"), estimate="CMLE", times=c(1,2))

stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(form.mle(y0,y1), coef(fit)[1], vcov(fit)[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1), coef(fit)[2], vcov(fit)[2,2]) <= tolerance)

# Force CMLE
for(prop.weight in prop.weights){
cat("====",prop.weight,"====\n")
set.seed(5432)
fit<-tergm(list(y0,y1) ~ Form(~edges) + Persist(~edges), constraints=~blockdiag("a"), estimate="CMLE", control=control.tergm(CMLE.ergm=control.ergm(MCMLE.effectiveSize = NULL, MCMC.samplesize = 2*1024, MCMC.burnin=10000, MCMC.interval = 1024, force.main=TRUE, MCMC.prop.weights=prop.weight)), times=c(1,2))

stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(form.mle(y0,y1), coef(fit)[1], vcov(fit)[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1), coef(fit)[2], vcov(fit)[2,2]) <= tolerance)
}

cat("Missing data:\n")

y1m<-network.copy(y1)
set.seed(765)
e <- as.edgelist(y1)[1,]
y1m[e[1], e[2]] <- NA
y1m[m,n] <- NA

# Force CMPLE
set.seed(765)
fit<-tergm(list(y0,y1m) ~ Form(~edges) + Persist(~edges), constraints=~blockdiag("a"), estimate="CMPLE", times=c(1,2))

stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(form.mle(y0,y1m), coef(fit)[1], vcov(fit)[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1m), coef(fit)[2], vcov(fit)[2,2]) <= tolerance)

# Autodetected CMPLE
set.seed(765)
fit<-tergm(list(y0,y1m) ~ Form(~edges) + Persist(~edges), constraints=~blockdiag("a"), estimate="CMLE", times=c(1,2))

stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(form.mle(y0,y1m), coef(fit)[1], vcov(fit)[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1m), coef(fit)[2], vcov(fit)[2,2]) <= tolerance)

# Force CMLE
for(prop.weight in prop.weights){
cat("====",prop.weight,"====\n")
set.seed(123456)
fit<-tergm(list(y0,y1m) ~ Form(~edges) + Persist(~edges), constraints=~blockdiag("a"), estimate="CMLE", control=control.tergm(CMLE.ergm=control.ergm(MCMLE.effectiveSize = NULL,  MCMC.samplesize = 2*1024, MCMC.burnin=10000, MCMC.interval = 1024, force.main=TRUE, MCMC.prop.weights=prop.weight)), times=c(1,2))

stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(form.mle(y0,y1m), coef(fit)[1], vcov(fit)[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1m), coef(fit)[2], vcov(fit)[2,2]) <= tolerance)
}

}, "dynamic MLE with block-diagonal constraints")
