# test simulate.networkDynamic using example code
require(networkDynamic)
require(tergm)

logit<-function(p)log(p/(1-p))
coef.form.f<-function(coef.diss,density) -log(((1+exp(coef.diss))/(density/(1-density)))-1)

# Construct a network with 20 nodes and 20 edges
n<-20
target.stats<-edges<-20
g0<-network.initialize(n,dir=TRUE)
g1<-san(g0~edges,target.stats=target.stats,verbose=TRUE)

S<-10

# To get an average duration of 10...
duration<-10
coef.diss<-logit(1-1/duration)

# To get an average of 20 edges...
dyads<-network.dyadcount(g1)
density<-edges/dyads
coef.form<-coef.form.f(coef.diss,density)



# Simulate a networkDynamic for five steps
dynsim<-simulate(g1,formation=~edges,dissolution=~edges,coef.form=coef.form,coef.diss=coef.diss,time.slices=S)

# ----- check net.obs.period encoding ----

# should have a net.obs.period object from 0 to 5
if (!all(unlist((dynsim%n%'net.obs.period')$observations)==c(0,1,1,11))){
  stop("simulate.network did not encode net.obs.period$observation correctly")
}

# "Resume" the simulation for five steps
dynsim2<-simulate(dynsim,time.slices=S)

if (!all(unlist((dynsim2%n%'net.obs.period')$observations)==c(0,1,1,11,11,21))){
  stop("simulate.networkDynamic did not encode net.obs.period$observation correctly")
}

