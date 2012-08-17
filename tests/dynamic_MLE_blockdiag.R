library(tergm)

tolerance<-0.05
n<-10
ns <- c(4,6)
theta<--1.5

logit<-function(p) log(p/(1-p))

block.dyadcount<-function(y){
  a <- rle(y %v% "a")
  sum(a$lengths*(a$lengths-1)) / (if(is.directed(y)) 1 else 2)
}

form.mle<-function(y0,y1){
  logit(network.edgecount(y1-y0,na.omit=TRUE)/(block.dyadcount(y0)-network.edgecount(y0-is.na(y1))-network.naedgecount(y1)))
}

diss.mle<-function(y0,y1){
  -logit(network.edgecount(y0-y1,na.omit=TRUE)/(network.edgecount(y0-is.na(y1))))
}

do.run <- function(dir){
y0<-network.initialize(n,dir=dir)
y0 %v% "a" <- rep(seq_along(ns), ns)
set.seed(321)
y0<-standardize.network(simulate(y0~edges, constraints=~blockdiag("a"), coef=theta, control=control.simulate(MCMC.burnin=n^2*2)))

cat("Complete data:\n")

set.seed(432)
y1<-standardize.network(simulate(y0~edges, constraints=~blockdiag("a"), coef=theta, control=control.simulate(MCMC.burnin=n^2*2)))

# Force CMPLE
set.seed(543)
fit<-stergm(list(y0,y1), formation=~edges, dissolution=~edges, constraints=~blockdiag("a"), estimate="CMPLE")

stopifnot(fit$estimate=="CMPLE", fit$formation.fit$estimate=="MPLE", fit$dissolution.fit$estimate=="MPLE")
stopifnot(all.equal(form.mle(y0,y1), coef(fit$formation.fit), tolerance=tolerance, check.attributes=FALSE))
stopifnot(all.equal(diss.mle(y0,y1), coef(fit$dissolution.fit), tolerance=tolerance, check.attributes=FALSE))

# Autodetected CMPLE
set.seed(543)
fit<-stergm(list(y0,y1), formation=~edges, dissolution=~edges, constraints=~blockdiag("a"), estimate="CMLE")

stopifnot(fit$estimate=="CMLE", is.null(fit$formation.fit$sample), is.null(fit$dissolution.fit$sample))
stopifnot(all.equal(form.mle(y0,y1), coef(fit$formation.fit), tolerance=tolerance, check.attributes=FALSE))
stopifnot(all.equal(diss.mle(y0,y1), coef(fit$dissolution.fit), tolerance=tolerance, check.attributes=FALSE))

# Force CMLE
set.seed(543)
fit<-stergm(list(y0,y1), formation=~edges, dissolution=~edges, constraints=~blockdiag("a"), estimate="CMLE", control=control.stergm(force.main=TRUE))

stopifnot(fit$estimate=="CMLE", fit$formation.fit$estimate=="MLE", fit$dissolution.fit$estimate=="MLE")
stopifnot(all.equal(form.mle(y0,y1), coef(fit$formation.fit), tolerance=tolerance, check.attributes=FALSE))
stopifnot(all.equal(diss.mle(y0,y1), coef(fit$dissolution.fit), tolerance=tolerance, check.attributes=FALSE))

cat("Missing data:\n")

y1m<-network.copy(y1)
set.seed(654)
y1m[as.edgelist(simulate(y0~edges, constraints=~blockdiag("a"), coef=theta, control=control.simulate(MCMC.burnin=n^2*2)))]<-NA

# Force CMPLE
set.seed(765)
fit<-stergm(list(y0,y1m), formation=~edges, dissolution=~edges, constraints=~blockdiag("a"), estimate="CMPLE")

stopifnot(fit$estimate=="CMPLE", fit$formation.fit$estimate=="MPLE", fit$dissolution.fit$estimate=="MPLE")
stopifnot(all.equal(form.mle(y0,y1m), coef(fit$formation.fit), tolerance=tolerance, check.attributes=FALSE))
stopifnot(all.equal(diss.mle(y0,y1m), coef(fit$dissolution.fit), tolerance=tolerance, check.attributes=FALSE))

# Autodetected CMPLE
set.seed(765)
fit<-stergm(list(y0,y1m), formation=~edges, dissolution=~edges, constraints=~blockdiag("a"), estimate="CMLE")

stopifnot(fit$estimate=="CMLE", is.null(fit$formation.fit$sample), is.null(fit$dissolution.fit$sample))
stopifnot(all.equal(form.mle(y0,y1m), coef(fit$formation.fit), tolerance=tolerance, check.attributes=FALSE))
stopifnot(all.equal(diss.mle(y0,y1m), coef(fit$dissolution.fit), tolerance=tolerance, check.attributes=FALSE))

# Force CMLE
set.seed(765)
fit<-stergm(list(y0,y1m), formation=~edges, dissolution=~edges, constraints=~blockdiag("a"), estimate="CMLE", control=control.stergm(force.main=TRUE))

stopifnot(fit$estimate=="CMLE", fit$formation.fit$estimate=="MLE", fit$dissolution.fit$estimate=="MLE")
stopifnot(all.equal(form.mle(y0,y1m), coef(fit$formation.fit), tolerance=tolerance, check.attributes=FALSE))
stopifnot(all.equal(diss.mle(y0,y1m), coef(fit$dissolution.fit), tolerance=tolerance, check.attributes=FALSE))
}

# Directed test
do.run(TRUE)
# Undirected test
do.run(FALSE)
