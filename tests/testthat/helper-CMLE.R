o <- options(tergm.eval.loglik=FALSE)

CMLE.tools <- new.env()

CMLE.tools$tolerance <- 3
CMLE.tools$n <- 10
CMLE.tools$m <- 6
CMLE.tools$theta <- -1.5

CMLE.tools$z.error <- function(truth, est, variance){
  if(abs(truth-est)<1e-6) 0 # Infinite case
  else abs(truth-est)/sqrt(variance)
}

CMLE.tools$logit<-function(p) log(p/(1-p))

CMLE.tools$form.mle<-function(y0,y1,y2){
  if(missing(y2)) logit(network.edgecount(y1-y0,na.omit=TRUE)/(network.dyadcount(y1)-network.edgecount(y0-is.na(y1))))
  else logit((network.edgecount(y1-y0,na.omit=TRUE) +
              network.edgecount(y2-y1,na.omit=TRUE))/
             (network.dyadcount(y1)-network.edgecount(y0-is.na(y1)) +
              network.dyadcount(y2)-network.edgecount(y1-is.na(y2))))
}

CMLE.tools$diss.mle<-function(y0,y1,y2){
  if(missing(y2)) -logit(network.edgecount(y0-y1,na.omit=TRUE)/(network.edgecount(y0-is.na(y1))))
  else -logit((network.edgecount(y0-y1,na.omit=TRUE) +
               network.edgecount(y1-y2,na.omit=TRUE))/
              (network.edgecount(y0-is.na(y1)) +
               network.edgecount(y1-is.na(y2))))
}

CMLE.tools$cross.mle<-function(y0,y1,y2){
  if(missing(y2)) logit(network.edgecount(y1, na.omit=TRUE)/network.dyadcount(y1, na.omit=TRUE))
  else logit((network.edgecount(y1, na.omit=TRUE) +
              network.edgecount(y2, na.omit=TRUE))/
             (network.dyadcount(y1, na.omit=TRUE) +
              network.dyadcount(y2, na.omit=TRUE)))
}

CMLE.tools$change.mle<-function(y0,y1,y2){
  if(missing(y2)) logit(network.edgecount((y0-y1)|(y1-y0), na.omit=TRUE)/network.dyadcount(y1, na.omit=TRUE))
  else logit((network.edgecount((y0-y1)|(y1-y0), na.omit=TRUE) +
              network.edgecount((y1-y2)|(y2-y1), na.omit=TRUE))/
             (network.dyadcount(y1, na.omit=TRUE) +
              network.dyadcount(y2, na.omit=TRUE)))
}

CMLE.tools$do.run_2 <- function(dir, bip=FALSE, prop.weights="default"){
  netdesc <- if(dir) "directed network" else {if(bip) "bipartite undirected network" else "undirected network"}

  if(bip){ # Extreme theta creates networks with too few ties to properly test.
    theta <- theta/2
  }
  
  y0<-network.initialize(n,dir=dir,bipartite=bip)
  set.seed(321)
  y0<-simulate(y0~edges, coef=theta, control=control.simulate(MCMC.burnin=n^2*2), dynamic=FALSE)

  completeness <- "completely observed"

  set.seed(123)
  y1<-simulate(y0~edges, coef=theta, control=control.simulate(MCMC.burnin=n^2*2), dynamic=FALSE)
  y2<-simulate(y1~edges, coef=theta, control=control.simulate(MCMC.burnin=n^2*2), dynamic=FALSE)

  test_that(paste("Force CMPLE on", completeness, netdesc), {

    set.seed(543)
    fit<-tergm(list(y0,y1,y2) ~ Form(~edges) + Persist(~edges), estimate="CMPLE", times=c(1,2,3))

    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(form.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
    expect_true(z.error(diss.mle(y0,y1,y2), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

    fit<-tergm(list(y0,y1,y2) ~ edges, estimate="CMPLE", times=c(1,2,3))
    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(cross.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1,y2) ~ Cross(~edges), estimate="CMPLE", times=c(1,2,3))
    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(cross.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1,y2) ~ Change(~edges), estimate="CMPLE", times=c(1,2,3))
    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(change.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

  })

  test_that(paste("Autodetect CMPLE on", completeness, netdesc), {

    set.seed(543)
    fit<-tergm(list(y0,y1,y2) ~ Form(~edges) + Persist(~edges), estimate="CMLE", times=c(1,2,3))

    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(form.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
    expect_true(z.error(diss.mle(y0,y1,y2), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

    fit<-tergm(list(y0,y1,y2) ~ edges, estimate="CMLE", times=c(1,2,3))
    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(cross.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1,y2) ~ Cross(~edges), estimate="CMLE", times=c(1,2,3))
    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(cross.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1,y2) ~ Change(~edges), estimate="CMLE", times=c(1,2,3))
    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(change.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

  })

  for(prop.weight in prop.weights){

    test_that(paste("Force CMPLE on", completeness, netdesc, "proposal", prop.weight), {

      ctrl <- control.tergm(CMLE.ergm=control.ergm(force.main=TRUE, MCMC.prop.weights=prop.weight))

      set.seed(543)
      fit<-tergm(list(y0,y1,y2) ~ Form(~edges) + Persist(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3))

      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(form.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
      expect_true(z.error(diss.mle(y0,y1,y2), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

      fit<-tergm(list(y0,y1,y2) ~ edges, estimate="CMLE", control=ctrl, times=c(1,2,3))
      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(cross.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

      fit<-tergm(list(y0,y1,y2) ~ Cross(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3))
      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(cross.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

      fit<-tergm(list(y0,y1,y2) ~ Change(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3))
      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(change.mle(y0,y1,y2), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
    })

  }

  completeness <- "partially observed"

  y2m<-network.copy(y2)
  set.seed(765)
  e <- as.edgelist(y2)[1,]
  y2m[e[1], e[2]] <- NA
  y2m[1,m+1] <- NA

  test_that(paste("Force CMPLE on", completeness, netdesc), {

    set.seed(765)
    fit<-tergm(list(y0,y1,y2m) ~ Form(~edges) + Persist(~edges), estimate="CMPLE", times=c(1,2,3))

    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(form.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
    expect_true(z.error(diss.mle(y0,y1,y2m), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

    fit<-tergm(list(y0,y1,y2m) ~ edges, estimate="CMPLE", times=c(1,2,3))
    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(cross.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1,y2m) ~ Cross(~edges), estimate="CMPLE", times=c(1,2,3))
    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(cross.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1,y2m) ~ Change(~edges), estimate="CMPLE", times=c(1,2,3))
    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(change.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

  })

  test_that(paste("Autodetect CMPLE on", completeness, netdesc), {

    set.seed(765)
    fit<-tergm(list(y0,y1,y2m) ~ Form(~edges) + Persist(~edges), estimate="CMLE", times=c(1,2,3))

    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(form.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
    expect_true(z.error(diss.mle(y0,y1,y2m), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

    fit<-tergm(list(y0,y1,y2m) ~ edges, estimate="CMLE", times=c(1,2,3))
    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(cross.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1,y2m) ~ Cross(~edges), estimate="CMLE", times=c(1,2,3))
    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(cross.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1,y2m) ~ Change(~edges), estimate="CMLE", times=c(1,2,3))
    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(change.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

  })
  
  for(prop.weight in prop.weights){

    test_that(paste("Force CMPLE on", completeness, netdesc, "proposal", prop.weight), {

      ctrl <- control.tergm(CMLE.ergm=control.ergm(force.main=TRUE, MCMC.prop.weights=prop.weight))

      set.seed(234)
      fit<-tergm(list(y0,y1,y2m) ~ Form(~edges) + Persist(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3))

      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(form.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
      expect_true(z.error(diss.mle(y0,y1,y2m), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

      fit<-tergm(list(y0,y1,y2m) ~ edges, estimate="CMLE", control=ctrl, times=c(1,2,3))
      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(cross.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

      fit<-tergm(list(y0,y1,y2m) ~ Cross(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3))
      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(cross.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

      fit<-tergm(list(y0,y1,y2m) ~ Change(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3))
      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(change.mle(y0,y1,y2m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
    })

  }
}

CMLE.tools$do.run_1 <- function(dir, bip=FALSE, prop.weights="default"){
  netdesc <- if(dir) "directed network" else {if(bip) "bipartite undirected network" else "undirected network"}

  if(bip){ # Extreme theta creates networks with too few ties to properly test.
    theta <- theta/2
  }
  
  y0<-network.initialize(n,dir=dir,bipartite=bip)
  set.seed(321)
  y0<-simulate(y0~edges, coef=theta, control=control.simulate(MCMC.burnin=n^2*2), dynamic=FALSE)

  completeness <- "completely observed"

  set.seed(123)
  y1<-simulate(y0~edges, coef=theta, control=control.simulate(MCMC.burnin=n^2*2), dynamic=FALSE)

  test_that(paste("Force CMPLE on", completeness, netdesc), {

    set.seed(543)
    fit<-tergm(list(y0,y1) ~ Form(~edges) + Persist(~edges), estimate="CMPLE", times=c(1,2))

    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(form.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
    expect_true(z.error(diss.mle(y0,y1), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

    fit<-tergm(list(y0,y1) ~ edges, estimate="CMPLE", times=c(1,2))
    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(cross.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1) ~ Cross(~edges), estimate="CMPLE", times=c(1,2))
    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(cross.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1) ~ Change(~edges), estimate="CMPLE", times=c(1,2))
    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(change.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

  })

  test_that(paste("Autodetect CMPLE on", completeness, netdesc), {

    set.seed(543)
    fit<-tergm(list(y0,y1) ~ Form(~edges) + Persist(~edges), estimate="CMLE", times=c(1,2))

    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(form.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
    expect_true(z.error(diss.mle(y0,y1), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

    fit<-tergm(list(y0,y1) ~ edges, estimate="CMLE", times=c(1,2))
    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(cross.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1) ~ Cross(~edges), estimate="CMLE", times=c(1,2))
    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(cross.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1) ~ Change(~edges), estimate="CMLE", times=c(1,2))
    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(change.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

  })

  for(prop.weight in prop.weights){

    test_that(paste("Force CMPLE on", completeness, netdesc, "proposal", prop.weight), {

      ctrl <- control.tergm(CMLE.ergm=control.ergm(force.main=TRUE, MCMC.prop.weights=prop.weight))

      set.seed(543)
      fit<-tergm(list(y0,y1) ~ Form(~edges) + Persist(~edges), estimate="CMLE", control=ctrl, times=c(1,2))

      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(form.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
      expect_true(z.error(diss.mle(y0,y1), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

      fit<-tergm(list(y0,y1) ~ edges, estimate="CMLE", control=ctrl, times=c(1,2))
      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(cross.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

      fit<-tergm(list(y0,y1) ~ Cross(~edges), estimate="CMLE", control=ctrl, times=c(1,2))
      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(cross.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

      fit<-tergm(list(y0,y1) ~ Change(~edges), estimate="CMLE", control=ctrl, times=c(1,2))
      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(change.mle(y0,y1), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
    })

  }

  completeness <- "partially observed"

  y1m<-network.copy(y1)
  set.seed(765)
  e <- as.edgelist(y1)[1,]
  y1m[e[1], e[2]] <- NA
  y1m[1,m+1] <- NA

  test_that(paste("Force CMPLE on", completeness, netdesc), {

    set.seed(765)
    fit<-tergm(list(y0,y1m) ~ Form(~edges) + Persist(~edges), estimate="CMPLE", times=c(1,2))

    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(form.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
    expect_true(z.error(diss.mle(y0,y1m), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

    fit<-tergm(list(y0,y1m) ~ edges, estimate="CMPLE", times=c(1,2))
    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(cross.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1m) ~ Cross(~edges), estimate="CMPLE", times=c(1,2))
    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(cross.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1m) ~ Change(~edges), estimate="CMPLE", times=c(1,2))
    expect_true(fit$estimate=="CMPLE")
    expect_true(z.error(change.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

  })

  test_that(paste("Autodetect CMPLE on", completeness, netdesc), {

    set.seed(765)
    fit<-tergm(list(y0,y1m) ~ Form(~edges) + Persist(~edges), estimate="CMLE", times=c(1,2))

    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(form.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
    expect_true(z.error(diss.mle(y0,y1m), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)

    fit<-tergm(list(y0,y1m) ~ edges, estimate="CMLE", times=c(1,2))
    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(cross.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1m) ~ Cross(~edges), estimate="CMLE", times=c(1,2))
    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(cross.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    fit<-tergm(list(y0,y1m) ~ Change(~edges), estimate="CMLE", times=c(1,2))
    expect_true(fit$estimate=="CMLE")
    expect_true(z.error(change.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

  })
  
  for(prop.weight in prop.weights){

    test_that(paste("Force CMPLE on", completeness, netdesc, "proposal", prop.weight), {

      ctrl <- control.tergm(CMLE.ergm=control.ergm(force.main=TRUE, MCMC.prop.weights=prop.weight))

      set.seed(234)
      fit<-tergm(list(y0,y1m) ~ Form(~edges) + Persist(~edges), estimate="CMLE", control=ctrl, times=c(1,2))

      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(form.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)
      expect_true(z.error(diss.mle(y0,y1m), coef(fit)[2], vcov(fit, sources="estimation")[2,2]) <= tolerance)


      fit<-tergm(list(y0,y1m) ~ edges, estimate="CMLE", control=ctrl, times=c(1,2))
      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(cross.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

      fit<-tergm(list(y0,y1m) ~ Cross(~edges), estimate="CMLE", control=ctrl, times=c(1,2))
      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(cross.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

      fit<-tergm(list(y0,y1m) ~ Change(~edges), estimate="CMLE", control=ctrl, times=c(1,2))
      expect_true(fit$estimate=="CMLE")
      expect_true(z.error(change.mle(y0,y1m), coef(fit)[1], vcov(fit, sources="estimation")[1,1]) <= tolerance)

    })

  }
}
