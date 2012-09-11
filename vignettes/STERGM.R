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
### code chunk number 3: STERGM.Snw:114-116
###################################################
data(samplk)
ls(pattern="samp*")


###################################################
### code chunk number 4: STERGM.Snw:120-124
###################################################
samp <- list()
samp[[1]] <- samplk1
samp[[2]] <- samplk2



###################################################
### code chunk number 5: STERGM.Snw:130-131
###################################################
plot(samplk1)


###################################################
### code chunk number 6: STERGM.Snw:145-150
###################################################
samp.stergm.fit <- stergm(samp,
	formation= ~edges+mutual+ctriad+ttriad,
	dissolution = ~edges+mutual+ctriad+ttriad,
	estimate = "CMLE"
	)


###################################################
### code chunk number 7: STERGM.Snw:164-167
###################################################

summary(samp.stergm.fit)



###################################################
### code chunk number 8: STERGM.Snw:189-192
###################################################
data(florentine)
ls()
flobusiness


###################################################
### code chunk number 9: STERGM.Snw:196-197
###################################################
plot(flobusiness)


###################################################
### code chunk number 10: STERGM.Snw:228-229
###################################################
theta.diss <- log(9)


###################################################
### code chunk number 11: STERGM.Snw:238-247
###################################################

flo.stergm.fit <- stergm(flobusiness,
	formation= ~edges+gwesp(0,fixed=T),
	dissolution = ~offset(edges),
	targets="formation",
	offset.coef.diss = theta.diss,
	estimate = "EGMME",
	control=control.stergm(SA.plot.progress=TRUE)
)


###################################################
### code chunk number 12: fit1diag
###################################################
mcmc.diagnostics(flo.stergm.fit)


###################################################
### code chunk number 13: STERGM.Snw:273-286
###################################################
flo.stergm.fit
names(flo.stergm.fit)
flo.stergm.fit$formation
flo.stergm.fit$formation.fit
summary(flo.stergm.fit)
\end{enumerate}

We have now obtained estimates for the coefficients of a formation model that, 
conditional on the stated dissolution model, yields simulated targets that matched those observed.
Something very useful we have also gained in the process is the ability to simulate networks with the 
desired cross-sectional structure and mean relational duration.  This ability serves us well for any application 
areas that requires us to simulate phenomena on dynamic networks, whether they entail the diffusion of information or disease, or some other process.  



###################################################
### code chunk number 14: STERGM.Snw:287-289
###################################################
flo.stergm.sim <- simulate.stergm(flo.stergm.fit, nsim=1, 
    time.slices = 1000)


###################################################
### code chunk number 15: STERGM.Snw:304-305
###################################################
flo.stergm.sim


###################################################
### code chunk number 16: STERGM.Snw:312-314
###################################################
net <- network.extract(flo.stergm.sim,at=429)
net


###################################################
### code chunk number 17: simex
###################################################
plot(network.extract(flo.stergm.sim,at=882))


###################################################
### code chunk number 18: STERGM.Snw:333-335
###################################################
summary(flobusiness~edges+gwesp(0,fixed=T))
colMeans(attributes(flo.stergm.sim)$stats)


###################################################
### code chunk number 19: statsform1
###################################################
plot(attributes(flo.stergm.sim)$stats)


###################################################
### code chunk number 20: STERGM.Snw:355-358
###################################################
flo.stergm.sim.dm <- as.data.frame(flo.stergm.sim)
names(flo.stergm.sim.dm)
mean(flo.stergm.sim.dm$duration)


###################################################
### code chunk number 21: STERGM.Snw:369-370
###################################################
get.edge.value(flo.stergm.sim, "active", unlist=FALSE)[[25]]


###################################################
### code chunk number 22: STERGM.Snw:382-383
###################################################
theta.diss.100 <- log(99)


###################################################
### code chunk number 23: STERGM.Snw:389-393
###################################################
flo.stergm.approx <- ergm(flobusiness~edges+gwesp(0,fixed=T))
summary(flo.stergm.approx)
theta.form <- flo.stergm.approx$coef 
theta.form


###################################################
### code chunk number 24: STERGM.Snw:399-400
###################################################
theta.form[1] <- theta.form[1] - theta.diss.100


###################################################
### code chunk number 25: STERGM.Snw:405-412
###################################################
flo.stergm.sim.2 <- simulate(flobusiness,
	formation=~edges+gwesp(0,fixed=T),
	dissolution=~edges,
	monitor="all",
	coef.form=theta.form,
	coef.diss=theta.diss.100,
	time.slices=10000)


###################################################
### code chunk number 26: simform
###################################################
summary(flobusiness~edges+gwesp(0,fixed=T))
colMeans(attributes(flo.stergm.sim.2)$stats)
stergm.sim.dm.2 <- as.data.frame(flo.stergm.sim.2)
mean(stergm.sim.dm.2$duration)
plot(attributes(flo.stergm.sim.2)$stats)


###################################################
### code chunk number 27: STERGM.Snw:458-461
###################################################
msm.net <- network.initialize(500, directed=F)	
msm.net %v% 'race' <- c(rep(0,250),rep(1,250))
msm.net


###################################################
### code chunk number 28: STERGM.Snw:467-470
###################################################
msm.form.formula <- ~edges+nodematch('race')+degree(0)+
    concurrent
msm.target.stats <- c(225,187,180,90)


###################################################
### code chunk number 29: STERGM.Snw:480-481
###################################################
msm.diss.formula <- ~offset(edges)+offset(nodematch("race"))


###################################################
### code chunk number 30: STERGM.Snw:509-510
###################################################
msm.theta.diss <- c(2.944, -0.747) 


###################################################
### code chunk number 31: STERGM.Snw:515-526
###################################################
set.seed(0)
msm.fit <- stergm(msm.net,
	formation= msm.form.formula,
	dissolution= msm.diss.formula,
	targets="formation",
	target.stats= msm.target.stats,
	offset.coef.diss = msm.theta.diss,
	estimate = "EGMME",
	control=control.stergm(SA.plot.progress=TRUE)
#		SA.phase2.levels.max=1)
)


###################################################
### code chunk number 32: msmdiag
###################################################
mcmc.diagnostics(msm.fit)


###################################################
### code chunk number 33: STERGM.Snw:547-548
###################################################
summary(msm.fit)


###################################################
### code chunk number 34: STERGM.Snw:553-554
###################################################
msm.sim <- simulate(msm.fit,time.slices=1000)


###################################################
### code chunk number 35: STERGM.Snw:559-561
###################################################
colMeans(attributes(msm.sim)$stats)
msm.target.stats


###################################################
### code chunk number 36: msmht
###################################################
msm.sim.dm <- as.data.frame(msm.sim)
plot(msm.sim.dm$head,msm.sim.dm$tail)


###################################################
### code chunk number 37: STERGM.Snw:579-587
###################################################
names(msm.sim.dm)
msm.sim.dm$race1 <- msm.sim.dm$head>250
msm.sim.dm$race2 <- msm.sim.dm$tail>250
msm.sim.dm$homoph <- msm.sim.dm$race1 == msm.sim.dm$race2
mean(msm.sim.dm$duration[msm.sim.dm$homoph==T & 
  msm.sim.dm$left.censored==F & msm.sim.dm$right.censored==F ])
mean(msm.sim.dm$duration[msm.sim.dm$homoph==F & 
  msm.sim.dm$left.censored==F & msm.sim.dm$right.censored==F ])


