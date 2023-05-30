#  File tests/degree.mean.age.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2023 Statnet Commons
################################################################################
library(statnet.common)
#opttest({
library(tergm)
set.seed(0)
logit<-function(p)log(p/(1-p))

coef.form.f<-function(coef.diss,density) -log(((1+exp(coef.diss))/(density/(1-density)))-1)

S<-300

n<-40
target.stats<-edges<-40
duration<-12
coef.diss<-logit(1-1/duration)

### Undirected

dyads<-n*(n-1)/2
density<-edges/dyads
coef.form<-coef.form.f(coef.diss,density)

cat("\nUndirected:\n")

g0<-network.initialize(n,dir=FALSE)

g0 %v% "a" <- rep(1:2, c(1,3)/4*n)

print(coef.form)
print(coef.diss)

# Simulate from the fit.
dynsim<-simulate(g0 ~ Form(~edges) + Persist(~edges),coef=c(coef.form,coef.diss),time.burnin=S, time.slices=S,verbose=TRUE,output="stats",
                 monitor=~edges+mean.age
                 +degree.mean.age(1:3)+degrange.mean.age(1:2,3:4)+degrange.mean.age(1:2)
                 +degree.mean.age(1:3,"a")+degrange.mean.age(1:2,3:4,"a")+degrange.mean.age(1:2,by="a"), dynamic=TRUE,
                 constraints=~.
                 )

dynsim.dup <- duplicated(as.data.frame(t(dynsim)))
dynsim <- dynsim[,!dynsim.dup]

targets <- c(edges,rep(12, ncol(dynsim)-1))
test <- approx.hotelling.diff.test(dynsim,mu0=targets)
if(test$p.value < 0.001){
  print(test)
  stop("At least one statistic differs from target.")
}
               

#dynsim<-simulate(g0,formation=~edges,dissolution=~edges+edges.ageinterval(7),coef.form=coef.form,coef.diss=c(logit(0.8),logit(0.8)-logit(0.7)),time.slices=S,verbose=TRUE,statsonly=TRUE,monitor=~edges+edges.ageinterval(1:40,2:41))


#print.sim.stats(dynsim,target.stats,duration)
#dynsim<-simulate(g0,formation=~edges,dissolution=~edges,coef.form=coef.form,coef.diss=coef.diss,time.slices=S,verbose=TRUE,statsonly=TRUE,monitor=~degrange(0:2,2:4,"a"))
#}, "degree mean age terms simulation")
