#  File tests/dynamic_EGMME.simple.R in package tergm, part of the
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
n<-20
do.plot <- FALSE
g0<-network.initialize(n,dir=FALSE)

#                    edges, mean.age
target.stats<-c(     n*1/2,       10)
logit<-function(p)log(p/(1-p))
coef.exact<-function(density,duration)
    list(form=-log(((1+exp(logit(1-1/duration)))/(density/(1-density)))-1),
         diss=logit(1-1/duration))


truth <- coef.exact(target.stats[1]/network.dyadcount(g0),
                    target.stats[2])

# Get a deliberately bad starting network.
set.seed(0)
g1<-san(g0~meandeg,target.stats=target.stats[1],verbose=TRUE)

# Fit the model with very poor starting values.
set.seed(1)
dynfit<-tergm(g1 ~ Form(~edges) + Persist(~edges), targets=~edges+mean.age, estimate="EGMME",target.stats=target.stats[-3], constraints=~., verbose=TRUE,control=control.tergm(SA.plot.progress=do.plot,SA.phase2.levels.min=2, SA.phase2.levels.max=4, SA.phase2.repeats = 10, SA.restart.on.err=FALSE,init=c(-log(.95/.05), 1)))

print(summary(dynfit))
mcmc.diagnostics(dynfit)

stopifnot(all.equal(unlist(truth),coef(dynfit),tol=0.05,check.attributes=FALSE))
#},"simple EGMME")
