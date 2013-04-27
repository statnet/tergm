library(statnet.common)
opttest({
rm(list=ls())
library(tergm)
{
data(samplk)

# Fit a transition from Time 1 to Time 2
samplk12 <- stergm(list(samplk1, samplk2),
                   formation=~edges+mutual+transitiveties+cyclicalties,
                   dissolution=~edges+mutual+transitiveties+cyclicalties,
                   estimate="CMLE")

samplk12.gof <- gof(samplk12)

samplk12.gof

summary(samplk12.gof)

plot(samplk12.gof)

plot(samplk12.gof, plotlogodds=TRUE)
}
},"gof.stergm.Rd")
opttest({
rm(list=ls())
library(tergm)
{
# EGMME Example
par(ask=FALSE)
n<-30
g0<-network.initialize(n,dir=FALSE)

#                     edges, degree(1), mean.age
target.stats<-c(      n*1/2,    n*0.6,        20)

dynfit<-stergm(g0,formation = ~edges+degree(1), dissolution = ~edges,
               targets = ~edges+degree(1)+mean.age,
               target.stats=target.stats, estimate="EGMME",
               control=control.stergm(SA.plot.progress=TRUE))

par(ask=TRUE)
mcmc.diagnostics(dynfit)
summary(dynfit)

# CMLE Example
data(samplk)

# Fit a transition from Time 1 to Time 2
samplk12 <- stergm(list(samplk1, samplk2),
                   formation=~edges+mutual+transitiveties+cyclicalties,
                   dissolution=~edges+mutual+transitiveties+cyclicalties,
                   estimate="CMLE")

mcmc.diagnostics(samplk12)
summary(samplk12)

# Fit a transition from Time 1 to Time 2 and from Time 2 to Time 3 jointly
samplk123 <- stergm(list(samplk1, samplk2, samplk3),
                    formation=~edges+mutual+transitiveties+cyclicalties,
                    dissolution=~edges+mutual+transitiveties+cyclicalties,
                    estimate="CMLE")

mcmc.diagnostics(samplk123)
summary(samplk123)
}
},"stergm.Rd")
