#  File tests/dynamic_MLE_blockdiag.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################
library(statnet.common)
#opttest({
library(tergm)
options(tergm.eval.loglik=FALSE)

tolerance<-0.05
n<-16
ns <- c(7,9)
theta<--1.5

z.error <- function(truth, est, variance){
  if(truth==est) 0 # Infinite case
  else abs(truth-est)/sqrt(variance)
}

logit<-function(p) log(p/(1-p))

block.dyadcount<-function(y){
  a <- rle(y %v% "a")
  sum(a$lengths*(a$lengths-1)) / (if(is.directed(y)) 1 else 2)
}

form.mle<-function(y0,y1){
  logit(network.edgecount(y1-y0,na.omit=TRUE)/(block.dyadcount(y0)-network.edgecount(y0-is.na(y1))-network.naedgecount(y1)))
}

diss.mle<-function(y0,y1){
  -logit(network.edgecount(y0-y1,na.omit=TRUE)/(network.edgecount(y0-is.na(y1))))
}

do.run <- function(dir, prop.weights="default"){
y0<-network.initialize(n,dir=dir)
y0 %v% "a" <- rep(seq_along(ns), ns)
set.seed(3213)
y0<-simulate(y0~edges, constraints=~blockdiag("a"), coef=theta, control=control.simulate(MCMC.burnin=n^2*2), dynamic=FALSE)

cat("Complete data:\n")

set.seed(123)
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
set.seed(543)
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
y1m[1,2] <- NA

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
set.seed(234)
fit<-tergm(list(y0,y1m) ~ Form(~edges) + Persist(~edges), constraints=~blockdiag("a"), estimate="CMLE", control=control.tergm(CMLE.ergm=control.ergm(MCMLE.effectiveSize = NULL, MCMC.samplesize = 2*1024, MCMC.burnin=10000, MCMC.interval = 1024, force.main=TRUE, MCMC.prop.weights=prop.weight)), times=c(1,2))

stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(form.mle(y0,y1m), coef(fit)[1], vcov(fit)[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1m), coef(fit)[2], vcov(fit)[2,2]) <= tolerance)
}
}

cat("=========== Directed test ===========\n")
do.run(TRUE, prop.weights=c("default","random"))
cat("=========== Undirected test ===========\n")
do.run(FALSE, prop.weights=c("default","random"))

#}, "dynamic MLE with block-diagonal constraints")
