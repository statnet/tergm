#  File tests/dynamic_MLE.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2021 Statnet Commons
################################################################################
library(tergm)
options(tergm.eval.loglik=FALSE)

tolerance<-3
n<-10
m<-6
theta<--1.5

z.error <- function(truth, est, variance){
  if(abs(truth-est)<1e-6) 0 # Infinite case
  else abs(truth-est)/sqrt(variance)
}

logit<-function(p) log(p/(1-p))

form.mle<-function(y0,y1){
  logit(network.edgecount(y1-y0,na.omit=TRUE)/(network.dyadcount(y1)-network.edgecount(y0-is.na(y1))))
}

diss.mle<-function(y0,y1){
  -logit(network.edgecount(y0-y1,na.omit=TRUE)/(network.edgecount(y0-is.na(y1))))
}

cross.mle<-function(y0,y1){
  logit(network.edgecount(y1, na.omit=TRUE)/network.dyadcount(y1, na.omit=TRUE))
}

change.mle<-function(y0,y1){
  logit(network.edgecount((y0-y1)|(y1-y0), na.omit=TRUE)/network.dyadcount(y1, na.omit=TRUE))
}

do.run <- function(dir, bip=FALSE, prop.weights="default"){
if(bip){ # Extreme theta creates networks with too few ties to properly test.
  theta <- theta/2
}
  
y0<-network.initialize(n,dir=dir,bipartite=bip)
set.seed(321)
y0<-simulate(y0~edges, coef=theta, control=control.simulate(MCMC.burnin=n^2*2), dynamic=FALSE)

cat("Complete data:\n")

set.seed(123)
y1<-simulate(y0~edges, coef=theta, control=control.simulate(MCMC.burnin=n^2*2), dynamic=FALSE)

# Force CMPLE
set.seed(543)
fit<-tergm(list(y0,y1) ~ Form(~edges) + Persist(~edges), estimate="CMPLE", times=c(1,2), eval.loglik=FALSE)

stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(form.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

fit<-tergm(list(y0,y1) ~ edges, estimate="CMPLE", times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(cross.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1) ~ Cross(~edges), estimate="CMPLE", times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(cross.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1) ~ Change(~edges), estimate="CMPLE", times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(change.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

# Autodetected CMPLE
set.seed(543)
fit<-tergm(list(y0,y1) ~ Form(~edges) + Persist(~edges), estimate="CMLE", times=c(1,2), eval.loglik=FALSE)

stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(form.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

fit<-tergm(list(y0,y1) ~ edges, estimate="CMLE", times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1) ~ Cross(~edges), estimate="CMLE", times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1) ~ Change(~edges), estimate="CMLE", times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(change.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

# Force CMLE
for(prop.weight in prop.weights){
cat("====",prop.weight,"====\n")

ctrl <- control.tergm(CMLE.ergm=control.ergm(force.main=TRUE, MCMC.prop.weights=prop.weight))

set.seed(543)
fit<-tergm(list(y0,y1) ~ Form(~edges) + Persist(~edges), estimate="CMLE", control=ctrl, times=c(1,2), eval.loglik=FALSE)

stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(form.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

fit<-tergm(list(y0,y1) ~ edges, estimate="CMLE", control=ctrl, times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1) ~ Cross(~edges), estimate="CMLE", control=ctrl, times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1) ~ Change(~edges), estimate="CMLE", control=ctrl, times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(change.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
}

cat("Missing data:\n")

y1m<-network.copy(y1)
set.seed(765)
e <- as.edgelist(y1)[1,]
y1m[e[1], e[2]] <- NA
y1m[1,m+1] <- NA

# Force CMPLE
set.seed(765)
fit<-tergm(list(y0,y1m) ~ Form(~edges) + Persist(~edges), estimate="CMPLE", times=c(1,2), eval.loglik=FALSE)

stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(form.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1m), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

fit<-tergm(list(y0,y1m) ~ edges, estimate="CMPLE", times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(cross.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1m) ~ Cross(~edges), estimate="CMPLE", times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(cross.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1m) ~ Change(~edges), estimate="CMPLE", times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(change.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

# Autodetected CMPLE
set.seed(765)
fit<-tergm(list(y0,y1m) ~ Form(~edges) + Persist(~edges), estimate="CMLE", times=c(1,2), eval.loglik=FALSE)

stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(form.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1m), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

fit<-tergm(list(y0,y1m) ~ edges, estimate="CMLE", times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1m) ~ Cross(~edges), estimate="CMLE", times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1m) ~ Change(~edges), estimate="CMLE", times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(change.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

# Force CMLE
for(prop.weight in prop.weights){
cat("====",prop.weight,"====\n")

ctrl <- control.tergm(CMLE.ergm=control.ergm(force.main=TRUE, MCMC.prop.weights=prop.weight))

set.seed(234)
fit<-tergm(list(y0,y1m) ~ Form(~edges) + Persist(~edges), estimate="CMLE", control=ctrl, times=c(1,2), eval.loglik=FALSE)

stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(form.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1m), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)


fit<-tergm(list(y0,y1m) ~ edges, estimate="CMLE", control=ctrl, times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1m) ~ Cross(~edges), estimate="CMLE", control=ctrl, times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1m) ~ Change(~edges), estimate="CMLE", control=ctrl, times=c(1,2), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(change.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

}
}

cat("=========== Directed test ===========\n")
do.run(TRUE, prop.weights=c("default","random"))
cat("=========== Undirected test ===========\n")
do.run(FALSE, prop.weights=c("default","random"))
cat("=========== Undirected bipartite test ===========\n")
do.run(FALSE, m, prop.weights=c("default","random"))
