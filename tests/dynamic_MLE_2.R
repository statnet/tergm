#  File tests/dynamic_MLE_2.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2021 Statnet Commons
################################################################################
library(statnet.common)
#opttest({
library(tergm)
options(tergm.eval.loglik=FALSE)

tolerance<-3
n<-10
m<-6
theta<--1.5

z.error <- function(truth, est, variance){
  if(abs(truth-est)<1e6) 0 # Infinite case
  else abs(truth-est)/sqrt(variance)
}

logit<-function(p) log(p/(1-p))

form.mle<-function(y0,y1,y2){
  logit((network.edgecount(y1-y0,na.omit=TRUE) +
         network.edgecount(y2-y1,na.omit=TRUE))/
        (network.dyadcount(y1)-network.edgecount(y0-is.na(y1)) +
         network.dyadcount(y2)-network.edgecount(y1-is.na(y2))))
}

diss.mle<-function(y0,y1,y2){
  -logit((network.edgecount(y0-y1,na.omit=TRUE) +
          network.edgecount(y1-y2,na.omit=TRUE))/
         (network.edgecount(y0-is.na(y1)) +
          network.edgecount(y1-is.na(y2))))
}

cross.mle<-function(y0,y1,y2){
  logit((network.edgecount(y1, na.omit=TRUE) +
         network.edgecount(y2, na.omit=TRUE))/
        (network.dyadcount(y1, na.omit=TRUE) +
         network.dyadcount(y2, na.omit=TRUE)))
}

change.mle<-function(y0,y1,y2){
  logit((network.edgecount((y0-y1)|(y1-y0), na.omit=TRUE) +
         network.edgecount((y1-y2)|(y2-y1), na.omit=TRUE))/
        (network.dyadcount(y1, na.omit=TRUE) +
         network.dyadcount(y2, na.omit=TRUE)))
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
y2<-simulate(y1~edges, coef=theta, control=control.simulate(MCMC.burnin=n^2*2), dynamic=FALSE)

# Force CMPLE
set.seed(543)
fit<-tergm(list(y0,y1,y2) ~ Form(~edges) + Persist(~edges), estimate="CMPLE", times=c(1,2,3))

stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(form.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1,y2), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

fit<-tergm(list(y0,y1,y2) ~ edges, estimate="CMPLE", times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(cross.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1,y2) ~ Cross(~edges), estimate="CMPLE", times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(cross.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1,y2) ~ Change(~edges), estimate="CMPLE", times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(change.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

# Autodetected CMPLE
set.seed(543)
fit<-tergm(list(y0,y1,y2) ~ Form(~edges) + Persist(~edges), estimate="CMLE", times=c(1,2,3))

stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(form.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1,y2), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

fit<-tergm(list(y0,y1,y2) ~ edges, estimate="CMLE", times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1,y2) ~ Cross(~edges), estimate="CMLE", times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1,y2) ~ Change(~edges), estimate="CMLE", times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(change.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

# Force CMLE
for(prop.weight in prop.weights){
cat("====",prop.weight,"====\n")

ctrl <- control.tergm(CMLE.ergm=control.ergm(force.main=TRUE, MCMC.prop.weights=prop.weight))

set.seed(543)
fit<-tergm(list(y0,y1,y2) ~ Form(~edges) + Persist(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3))

stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(form.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1,y2), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

fit<-tergm(list(y0,y1,y2) ~ edges, estimate="CMLE", control=ctrl, times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1,y2) ~ Cross(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1,y2) ~ Change(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(change.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
}

cat("Missing data:\n")

y2m<-network.copy(y2)
set.seed(765)
e <- as.edgelist(y2)[1,]
y2m[e[1], e[2]] <- NA
y2m[1,m+1] <- NA

# Force CMPLE
set.seed(765)
fit<-tergm(list(y0,y1,y2m) ~ Form(~edges) + Persist(~edges), estimate="CMPLE", times=c(1,2,3))

stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(form.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1,y2m), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

fit<-tergm(list(y0,y1,y2m) ~ edges, estimate="CMPLE", times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(cross.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1,y2m) ~ Cross(~edges), estimate="CMPLE", times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(cross.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1,y2m) ~ Change(~edges), estimate="CMPLE", times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMPLE")
stopifnot(z.error(change.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

# Autodetected CMPLE
set.seed(765)
fit<-tergm(list(y0,y1,y2m) ~ Form(~edges) + Persist(~edges), estimate="CMLE", times=c(1,2,3))

stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(form.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1,y2m), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

fit<-tergm(list(y0,y1,y2m) ~ edges, estimate="CMLE", times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1,y2m) ~ Cross(~edges), estimate="CMLE", times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1,y2m) ~ Change(~edges), estimate="CMLE", times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(change.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

# Force CMLE
for(prop.weight in prop.weights){
cat("====",prop.weight,"====\n")
set.seed(234)
fit<-tergm(list(y0,y1,y2m) ~ Form(~edges) + Persist(~edges), estimate="CMLE", control=control.tergm(CMLE.ergm=control.ergm(force.main=TRUE, MCMC.prop.weights=prop.weight)), times=c(1,2,3))

stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(form.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
stopifnot(z.error(diss.mle(y0,y1,y2m), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

fit<-tergm(list(y0,y1,y2m) ~ edges, estimate="CMLE", control=ctrl, times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1,y2m) ~ Cross(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(cross.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

fit<-tergm(list(y0,y1,y2m) ~ Change(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3), eval.loglik=FALSE)
stopifnot(fit$estimate=="CMLE")
stopifnot(z.error(change.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
}
}

cat("=========== Directed test ===========\n")
do.run(TRUE, prop.weights=c("default","random"))
cat("=========== Undirected test ===========\n")
do.run(FALSE, prop.weights=c("default","random"))
cat("=========== Undirected bipartite test ===========\n")
do.run(FALSE, m, prop.weights=c("default","random"))



#}, "dynamic MLE with two transitions")
