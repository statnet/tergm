library(ergm)
library(coda)

n<-50

do.plot <- FALSE

g0<-network.initialize(n,dir=FALSE)

#            meandeg, degree(1)
target.stats<-c(      n*1/2,    n*0.6)

# Get a reasonably close starting network.
set.seed(9)
g1<-san(g0~meandeg+degree(1),target.stats=target.stats,verbose=TRUE)

coef.form <- c(-6.801597, 1.016334)
coef.diss <- c(2.944439)

# Fit the model with very poor starting values.
set.seed(123)
dynfit<-stergm(g1,formation=~edges+degree(1),dissolution=~offset(edges), targets="formation", estimate="EGMME", offset.coef.diss=log(.95/.05),target.stats=target.stats,verbose=TRUE,control=control.stergm(SA.plot.progress=do.plot,init.form=c(-log(.95/.05),0)))

print(summary(dynfit))
mcmc.diagnostics(dynfit)

# Simulate from the fit.
set.seed(321)
dynsim<-simulate(dynfit,time.slices=1000,monitor="all",verbose=TRUE)

# Print out the resulting target.stats and their t-value w.r.t. the
# target target.stats.
print(target.stats)
stats.sim <- attr(dynsim,"stats")
target.stats.sim<-apply(stats.sim,2,mean)
print(target.stats.sim)
print(effectiveSize(stats.sim))
print((target.stats.sim-target.stats)/sqrt(apply(stats.sim,2,var)/effectiveSize(stats.sim)))

# Simulate from an equivalent fit.
set.seed(456)
dynsim<-simulate(g1, time.slices=1000, formation=~edges+degree(1), dissolution=~dyadcov(matrix(1,n,n))+edges, monitor=~edges+edge.ages,
                 coef.form=coef.form, coef.diss=c(1,coef.diss-1), statsonly=TRUE, verbose=TRUE)

# Average age of an extant edge.
print(colMeans(dynsim)[2]/colMeans(dynsim)[1])

# Fit a model with some unusual offset and target configurations.

# Fix formation, fit dissolution.
set.seed(654)
dynfit<-stergm(g1,formation=~offset(edges)+offset(degree(1)),dissolution=~edges, targets="dissolution", estimate="EGMME", offset.coef.form=coef.form,target.stats=target.stats[1],control=control.stergm(SA.plot.progress=do.plot, init.diss=1))

print(summary(dynfit))
mcmc.diagnostics(dynfit)

# Fix formation edges and dissolution, fit degree(1) with edges as targets.
set.seed(987)
dynfit<-stergm(g1,formation=~offset(edges)+degree(1),dissolution=~offset(edges), targets=~edges, estimate="EGMME", offset.coef.diss=coef.diss,target.stats=target.stats[1],control=control.stergm(SA.plot.progress=do.plot, init.form=c(coef.form[1],0)))

print(summary(dynfit))
mcmc.diagnostics(dynfit)

# Fix formation edges and dissolution, fit degree(1) with edges and degree(1) as target (overspecified model).
set.seed(789)
dynfit<-stergm(g1,formation=~offset(edges)+degree(1),dissolution=~offset(edges), targets="formation", estimate="EGMME", offset.coef.diss=coef.diss,target.stats=target.stats,control=control.stergm(SA.plot.progress=do.plot, init.form=c(coef.form[1],0)))

print(summary(dynfit))
mcmc.diagnostics(dynfit)

# Fix formation edges, fit degree(1) formation and edges dissolution with edges and degree(1) as target.
set.seed(12345)
dynfit<-stergm(g1,formation=~offset(edges)+degree(1),dissolution=~edges, targets="formation", estimate="EGMME", target.stats=target.stats,control=control.stergm(SA.plot.progress=do.plot,init.form=c(coef.form[1],0),init.diss=1))
print(summary(dynfit))
mcmc.diagnostics(dynfit)

# All parameters free, edges, degree(1), and edge.ages as target.
# I am not 100% sure the target statistic is correct.
set.seed(56789)
dynfit<-stergm(g1,formation=~edges+degree(1),dissolution=~edges, targets=~edges+degree(1)+edge.ages, estimate="EGMME", target.stats=c(target.stats,target.stats[1]*20),control=control.stergm(SA.plot.progress=do.plot,init.form=c(-1,0),init.diss=1))

print(summary(dynfit))
mcmc.diagnostics(dynfit)