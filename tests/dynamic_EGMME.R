#  File tests/dynamic_EGMME.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2019 Statnet Commons
#######################################################################
library(statnet.common)
opttest({
library(tergm)
n<-40
do.plot <- TRUE
g0<-network.initialize(n,dir=FALSE)

#                    edges, degree(1), mean.age
target.stats<-c(      n*1/2,    n*0.6,       20)

# Get a deliberately bad starting network.
g1<-san(g0~meandeg+degree(1),target.stats=target.stats[-3],verbose=TRUE)

coef.form <- c(-6.57, 1.01)
coef.diss <- c(2.944439)

# Fit the model with very poor starting values.
set.seed(3)
dynfit<-tergm(g1 ~ FormE(~edges + degree(1)) + offset(DissE(~edges)), targets=~edges + degree(1), estimate="EGMME", constraints="discordTNT"~., offset.coef=coef.diss,target.stats=target.stats[-3],verbose=TRUE,control=control.tergm(SA.plot.progress=do.plot,SA.restart.on.err=FALSE,init=c(-log(.95/.05),0, coef.diss)))

#print(summary(dynfit))
#mcmc.diagnostics(dynfit)

stopifnot(all.equal(c(coef.form,coef.diss),dynfit$fit$coef,tol=0.01,check.attributes=FALSE))

# All parameters free, edges, degree(1), and edge.ages as target.
set.seed(5)
dynfit2<-tergm(g1 ~ FormE(~edges + degree(1)) + DissE(~edges), targets=~edges+degree(1)+mean.age, estimate="EGMME", constraints="discordTNT"~., target.stats=target.stats,control=control.tergm(SA.plot.progress=do.plot,SA.plot.stats=TRUE))

#print(summary(dynfit2))
#mcmc.diagnostics(dynfit2)

stopifnot(all.equal(c(coef.form,coef.diss),dynfit2$fit$coef,tol=0.01,check.attributes=FALSE))
}, "EGMME")
