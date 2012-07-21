### R code from vignette source 'STERGM.Snw'

###################################################
### code chunk number 1: STERGM.Snw:22-23
###################################################
options(width=60)


###################################################
### code chunk number 2: STERGM.Snw:26-27
###################################################
options(continue=" ")


###################################################
### code chunk number 3: STERGM.Snw:44-48
###################################################
#install.packages("ergm")
#install.packages("networkDynamic")
library(ergm)
library(networkDynamic)


###################################################
### code chunk number 4: STERGM.Snw:83-86
###################################################
library(tergm)
data("florentine")
ls()


###################################################
### code chunk number 5: flobusplot
###################################################
plot(flobusiness)


###################################################
### code chunk number 6: STERGM.Snw:104-106
###################################################
fit1 <- ergm(flobusiness~edges+gwesp(0,fixed=T))
summary(fit1)


###################################################
### code chunk number 7: STERGM.Snw:111-112
###################################################
sim1 <- simulate(fit1,nsim=1,control=control.simulate.ergm(MCMC.burnin=1000))


###################################################
### code chunk number 8: sim1plot
###################################################
 plot(sim1)


###################################################
### code chunk number 9: STERGM.Snw:328-329
###################################################
theta.diss <- log(9)


###################################################
### code chunk number 10: STERGM.Snw:337-345
###################################################

stergm.fit.1 <- stergm(flobusiness,formation= ~edges+gwesp(0,fixed=T),dissolution = ~offset(edges),targets="formation",offset.coef.diss = theta.diss,estimate = "EGMME",verbose=2,control=control.stergm(SA.plot.progress=TRUE))


###################################################
### code chunk number 11: fit1diag
###################################################
mcmc.diagnostics(stergm.fit.1)


###################################################
### code chunk number 12: STERGM.Snw:364-369
###################################################
stergm.fit.1
names(stergm.fit.1)
stergm.fit.1$formation
stergm.fit.1$formation.fit
summary(stergm.fit.1)


###################################################
### code chunk number 13: STERGM.Snw:380-381
###################################################
stergm.sim.1 <- simulate.stergm(stergm.fit.1, nsim=1, time.slices = 1000)


###################################################
### code chunk number 14: STERGM.Snw:397-398
###################################################
stergm.sim.1


###################################################
### code chunk number 15: STERGM.Snw:405-406
###################################################
network.extract(stergm.sim.1,at=429)


###################################################
### code chunk number 16: simex
###################################################
plot(network.extract(stergm.sim.1,at=882))


###################################################
### code chunk number 17: STERGM.Snw:425-427
###################################################
summary(flobusiness~edges+gwesp(0,fixed=T))
colMeans(attributes(stergm.sim.1)$stats)


###################################################
### code chunk number 18: statsform1
###################################################
plot(attributes(stergm.sim.1)$stats)


###################################################
### code chunk number 19: STERGM.Snw:450-453
###################################################
stergm.sim.1.dm <- as.data.frame(stergm.sim.1)
names(stergm.sim.1.dm)
mean(stergm.sim.1.dm$duration)


###################################################
### code chunk number 20: STERGM.Snw:464-465
###################################################
stergm.sim.1$mel[[25]]$atl$active


###################################################
### code chunk number 21: STERGM.Snw:470-471
###################################################
stergm.sim.1 %e% "active"


###################################################
### code chunk number 22: STERGM.Snw:531-532
###################################################
theta.diss.100 <- log(99)


###################################################
### code chunk number 23: STERGM.Snw:538-541
###################################################
summary(fit1)
theta.form <- fit1$coef 
theta.form


###################################################
### code chunk number 24: STERGM.Snw:547-548
###################################################
theta.form[1] <- theta.form[1] - theta.diss.100


###################################################
### code chunk number 25: STERGM.Snw:553-560
###################################################
stergm.sim.2 <- simulate(flobusiness,
	formation=~edges+gwesp(0,fixed=T),
	dissolution=~edges,
	monitor="all",
	coef.form=theta.form,
	coef.diss=theta.diss.100,
	time.slices=10000)


###################################################
### code chunk number 26: STERGM.Snw:565-569
###################################################
summary(flobusiness~edges+gwesp(0,fixed=T))
colMeans(attributes(stergm.sim.2)$stats)
stergm.sim.dm.2 <- as.data.frame(stergm.sim.2)
mean(stergm.sim.dm.2$duration)


###################################################
### code chunk number 27: simform
###################################################
plot(attributes(stergm.sim.2)$stats)


###################################################
### code chunk number 28: STERGM.Snw:593-598
###################################################
data(samplk)
ls(pattern="samp*")
samp <- list()
samp[[1]] <- samplk1
samp[[2]] <- samplk2


###################################################
### code chunk number 29: STERGM.Snw:615-616
###################################################
plot(samplk1)


###################################################
### code chunk number 30: STERGM.Snw:634-639
###################################################
stergm.fit.3 <- stergm(samp,
	formation= ~edges+mutual+ctriad+ttriad,
	dissolution = ~edges+mutual+ctriad+ttriad,
	estimate = "CMLE"
	)


###################################################
### code chunk number 31: STERGM.Snw:649-650
###################################################
summary(stergm.fit.3)


###################################################
### code chunk number 32: STERGM.Snw:687-690
###################################################
msm.net <- network.initialize(500, directed=F)	
msm.net %v% 'race' <- c(rep(0,250),rep(1,250))
msm.net


###################################################
### code chunk number 33: STERGM.Snw:696-699
###################################################
msm.form.formula <- ~edges+nodematch('race')+degree(0)+concurrent
msm.target.stats <- c(225,187,180,90)



###################################################
### code chunk number 34: STERGM.Snw:709-710
###################################################
msm.diss.formula <- ~offset(edges)+offset(nodematch("race"))


###################################################
### code chunk number 35: STERGM.Snw:744-745
###################################################
msm.theta.diss <- c(2.944, -0.747) 


###################################################
### code chunk number 36: STERGM.Snw:750-759
###################################################
msm.fit <- stergm(msm.net,
	formation= msm.form.formula,
	dissolution= msm.diss.formula,
	targets="formation",
	target.stats= msm.target.stats,
	offset.coef.diss = msm.theta.diss,
	estimate = "EGMME",
	control=control.stergm(SA.plot.progress=TRUE,SA.init.gain=0.005)
)


###################################################
### code chunk number 37: msmdiag
###################################################
mcmc.diagnostics(msm.fit)


###################################################
### code chunk number 38: STERGM.Snw:776-777
###################################################
summary(msm.fit)


###################################################
### code chunk number 39: STERGM.Snw:782-783
###################################################
msm.sim <- simulate(msm.fit,time.slices=1000)


###################################################
### code chunk number 40: STERGM.Snw:788-790
###################################################
colMeans(attributes(msm.sim)$stats)
msm.target.stats


###################################################
### code chunk number 41: STERGM.Snw:795-796
###################################################
msm.sim.dm <- as.data.frame(msm.sim)


###################################################
### code chunk number 42: msmht
###################################################
plot(msm.sim.dm$head,msm.sim.dm$tail)


###################################################
### code chunk number 43: STERGM.Snw:811-817
###################################################
names(msm.sim.dm)
msm.sim.dm$race1 <- msm.sim.dm$head>250
msm.sim.dm$race2 <- msm.sim.dm$tail>250
msm.sim.dm$homoph <- msm.sim.dm$race1 == msm.sim.dm$race2
mean(msm.sim.dm$duration[msm.sim.dm$homoph==T & msm.sim.dm$left.censored==F & msm.sim.dm$right.censored==F ])
mean(msm.sim.dm$duration[msm.sim.dm$homoph==F & msm.sim.dm$left.censored==F & msm.sim.dm$right.censored==F ])



data(florentine)
theta.diss <- log(9)
stergm.fit.1 <- stergm(flomarriage, formation= ~edges+nodecov('wealth'), dissolution = ~offset(edges), targets="formation", offset.coef.diss = theta.diss, estimate = "EGMME",control=control.stergm(SA.plot.progress=TRUE),verbose=2)
