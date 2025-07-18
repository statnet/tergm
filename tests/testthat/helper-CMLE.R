#  File tests/testthat/helper-CMLE.R in package tergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################
logit<-function(p) log(p/(1-p))
ilogit<-function(x) 1/(1+exp(-x))

options(tergm.eval.loglik = FALSE, useFancyQuotes = FALSE)

CMLE.tools <- new.env()

CMLE.tools$tolerance <- 4
CMLE.tools$n <- 10
CMLE.tools$m <- 6
CMLE.tools$theta <- -1.5

expect_logit_ztest <- function(object, expected, variance = 0, alpha = 0.001, p_tolerance = 1e-6){
  eval.parent(call("expect_named", substitute(object), names(expected)))

  act <- quasi_label(rlang::enquo(object), arg = "object")
  exp <- quasi_label(rlang::enquo(expected), arg = "expected")
  var <- quasi_label(rlang::enquo(variance), arg = "variance")

  if(abs(ilogit(exp$val)-ilogit(act$val)) < p_tolerance){
    succeed() # Infinite case
    return(invisible(act$val))
  }

  act$z <- (act$val - exp$val) / sqrt(var$val)
  act$p <- 2 * pnorm(-abs(act$z))

  expect(act$p >= alpha,
         sprintf("z = ( %s - %s ) / sqrt(%s) = ( %f - %f ) / %f = %f; P(|Z|>|z|) = %0.3f < %0.2f",
                 act$lab, exp$lab, var$lab,
                 act$val, exp$val, sqrt(var$val),
                 act$z, act$p, alpha))

  invisible(act$val)
}

CMLE.tools$z.error <- function(truth, est, variance){
  if(abs(truth-est)<1e-6 || abs(ilogit(truth)-ilogit(est))<1e-6) 0 # Infinite case
  else abs(truth-est)/sqrt(variance)
}

CMLE.tools$form.mle<-function(y0,y1,y2){
  setNames(
    if(missing(y2)) logit(network.edgecount(y1-y0,na.omit=TRUE)/(network.dyadcount(y1)-network.edgecount(y0-is.na(y1))))
    else logit((network.edgecount(y1-y0,na.omit=TRUE) +
                network.edgecount(y2-y1,na.omit=TRUE))/
               (network.dyadcount(y1)-network.edgecount(y0-is.na(y1)) +
                network.dyadcount(y2)-network.edgecount(y1-is.na(y2)))),
    "Form(1)~edges")
}

CMLE.tools$diss.mle<-function(y0,y1,y2){
  setNames(
    if(missing(y2)) -logit(network.edgecount(y0-y1,na.omit=TRUE)/(network.edgecount(y0-is.na(y1))))
    else -logit((network.edgecount(y0-y1,na.omit=TRUE) +
                 network.edgecount(y1-y2,na.omit=TRUE))/
                (network.edgecount(y0-is.na(y1)) +
                 network.edgecount(y1-is.na(y2)))),
    "Persist(1)~edges")
}

CMLE.tools$cross.mle<-function(y0,y1,y2){
  setNames(
    if(missing(y2)) logit(network.edgecount(y1, na.omit=TRUE)/network.dyadcount(y1, na.omit=TRUE))
    else logit((network.edgecount(y1, na.omit=TRUE) +
                network.edgecount(y2, na.omit=TRUE))/
               (network.dyadcount(y1, na.omit=TRUE) +
                network.dyadcount(y2, na.omit=TRUE))),
    "Cross(1)~edges"
  )
}

CMLE.tools$plain.mle<-function(y0,y1,y2){
  setNames(CMLE.tools$cross.mle(y0, y1, y2), "edges")
}

CMLE.tools$change.mle<-function(y0,y1,y2){
  setNames(
    if(missing(y2)) logit(network.edgecount((y0-y1)|(y1-y0), na.omit=TRUE)/network.dyadcount(y1, na.omit=TRUE))
    else logit((network.edgecount((y0-y1)|(y1-y0), na.omit=TRUE) +
                network.edgecount((y1-y2)|(y2-y1), na.omit=TRUE))/
               (network.dyadcount(y1, na.omit=TRUE) +
                network.dyadcount(y2, na.omit=TRUE))),
    "Change(1)~edges")
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

    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], form.mle(y0,y1,y2))
    expect_logit_ztest(coef(fit)[2], diss.mle(y0,y1,y2))

    fit<-tergm(list(y0,y1,y2) ~ edges, estimate="CMPLE", times=c(1,2,3))
    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], plain.mle(y0,y1,y2))

    fit<-tergm(list(y0,y1,y2) ~ Cross(~edges), estimate="CMPLE", times=c(1,2,3))
    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], cross.mle(y0,y1,y2))

    fit<-tergm(list(y0,y1,y2) ~ Change(~edges), estimate="CMPLE", times=c(1,2,3))
    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], change.mle(y0,y1,y2))

  })

  test_that(paste("Autodetect CMPLE on", completeness, netdesc), {

    set.seed(543)
    fit<-tergm(list(y0,y1,y2) ~ Form(~edges) + Persist(~edges), estimate="CMLE", times=c(1,2,3))

    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], form.mle(y0,y1,y2))
    expect_logit_ztest(coef(fit)[2], diss.mle(y0,y1,y2))

    fit<-tergm(list(y0,y1,y2) ~ edges, estimate="CMLE", times=c(1,2,3))
    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], plain.mle(y0,y1,y2))

    fit<-tergm(list(y0,y1,y2) ~ Cross(~edges), estimate="CMLE", times=c(1,2,3))
    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], cross.mle(y0,y1,y2))

    fit<-tergm(list(y0,y1,y2) ~ Change(~edges), estimate="CMLE", times=c(1,2,3))
    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], change.mle(y0,y1,y2))

  })

  for(prop.weight in prop.weights){

    test_that(paste("Force CMLE on", completeness, netdesc, "proposal", prop.weight), {

      ctrl <- control.tergm(CMLE.ergm=control.ergm(force.main=TRUE, MCMC.prop.weights=prop.weight))

      set.seed(543)
      fit<-tergm(list(y0,y1,y2) ~ Form(~edges) + Persist(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3))

      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], form.mle(y0,y1,y2), vcov(fit, sources="estimation")[1,1])
      expect_logit_ztest(coef(fit)[2], diss.mle(y0,y1,y2), vcov(fit, sources="estimation")[2,2])

      fit<-tergm(list(y0,y1,y2) ~ edges, estimate="CMLE", control=ctrl, times=c(1,2,3))
      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], plain.mle(y0,y1,y2), vcov(fit, sources="estimation")[1,1])

      fit<-tergm(list(y0,y1,y2) ~ Cross(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3))
      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], cross.mle(y0,y1,y2), vcov(fit, sources="estimation")[1,1])

      fit<-tergm(list(y0,y1,y2) ~ Change(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3))
      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], change.mle(y0,y1,y2), vcov(fit, sources="estimation")[1,1])
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

    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], form.mle(y0,y1,y2m))
    expect_logit_ztest(coef(fit)[2], diss.mle(y0,y1,y2m))

    fit<-tergm(list(y0,y1,y2m) ~ edges, estimate="CMPLE", times=c(1,2,3))
    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], plain.mle(y0,y1,y2m))

    fit<-tergm(list(y0,y1,y2m) ~ Cross(~edges), estimate="CMPLE", times=c(1,2,3))
    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], cross.mle(y0,y1,y2m))

    fit<-tergm(list(y0,y1,y2m) ~ Change(~edges), estimate="CMPLE", times=c(1,2,3))
    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], change.mle(y0,y1,y2m))

  })

  test_that(paste("Autodetect CMPLE on", completeness, netdesc), {

    set.seed(765)
    fit<-tergm(list(y0,y1,y2m) ~ Form(~edges) + Persist(~edges), estimate="CMLE", times=c(1,2,3))

    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], form.mle(y0,y1,y2m))
    expect_logit_ztest(coef(fit)[2], diss.mle(y0,y1,y2m))

    fit<-tergm(list(y0,y1,y2m) ~ edges, estimate="CMLE", times=c(1,2,3))
    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], plain.mle(y0,y1,y2m))

    fit<-tergm(list(y0,y1,y2m) ~ Cross(~edges), estimate="CMLE", times=c(1,2,3))
    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], cross.mle(y0,y1,y2m))

    fit<-tergm(list(y0,y1,y2m) ~ Change(~edges), estimate="CMLE", times=c(1,2,3))
    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], change.mle(y0,y1,y2m))

  })
  
  for(prop.weight in prop.weights){

    test_that(paste("Force CMLE on", completeness, netdesc, "proposal", prop.weight), {

      ctrl <- control.tergm(CMLE.ergm=control.ergm(force.main=TRUE, MCMC.prop.weights=prop.weight))

      set.seed(234)
      fit<-tergm(list(y0,y1,y2m) ~ Form(~edges) + Persist(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3))

      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], form.mle(y0,y1,y2m), vcov(fit, sources="estimation")[1,1])
      expect_logit_ztest(coef(fit)[2], diss.mle(y0,y1,y2m), vcov(fit, sources="estimation")[2,2])

      fit<-tergm(list(y0,y1,y2m) ~ edges, estimate="CMLE", control=ctrl, times=c(1,2,3))
      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], plain.mle(y0,y1,y2m), vcov(fit, sources="estimation")[1,1])

      fit<-tergm(list(y0,y1,y2m) ~ Cross(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3))
      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], cross.mle(y0,y1,y2m), vcov(fit, sources="estimation")[1,1])

      fit<-tergm(list(y0,y1,y2m) ~ Change(~edges), estimate="CMLE", control=ctrl, times=c(1,2,3))
      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], change.mle(y0,y1,y2m), vcov(fit, sources="estimation")[1,1])
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

    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], form.mle(y0,y1))
    expect_logit_ztest(coef(fit)[2], diss.mle(y0,y1))

    fit<-tergm(list(y0,y1) ~ edges, estimate="CMPLE", times=c(1,2))
    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], plain.mle(y0,y1))

    fit<-tergm(list(y0,y1) ~ Cross(~edges), estimate="CMPLE", times=c(1,2))
    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], cross.mle(y0,y1))

    fit<-tergm(list(y0,y1) ~ Change(~edges), estimate="CMPLE", times=c(1,2))
    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], change.mle(y0,y1))

  })

  test_that(paste("Autodetect CMPLE on", completeness, netdesc), {

    set.seed(543)
    fit<-tergm(list(y0,y1) ~ Form(~edges) + Persist(~edges), estimate="CMLE", times=c(1,2))

    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], form.mle(y0,y1))
    expect_logit_ztest(coef(fit)[2], diss.mle(y0,y1))

    fit<-tergm(list(y0,y1) ~ edges, estimate="CMLE", times=c(1,2))
    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], plain.mle(y0,y1))

    fit<-tergm(list(y0,y1) ~ Cross(~edges), estimate="CMLE", times=c(1,2))
    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], cross.mle(y0,y1))

    fit<-tergm(list(y0,y1) ~ Change(~edges), estimate="CMLE", times=c(1,2))
    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], change.mle(y0,y1))

  })

  for(prop.weight in prop.weights){

    test_that(paste("Force CMLE on", completeness, netdesc, "proposal", prop.weight), {

      ctrl <- control.tergm(CMLE.ergm=control.ergm(force.main=TRUE, MCMC.prop.weights=prop.weight))

      set.seed(543)
      fit<-tergm(list(y0,y1) ~ Form(~edges) + Persist(~edges), estimate="CMLE", control=ctrl, times=c(1,2))

      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], form.mle(y0,y1), vcov(fit, sources="estimation")[1,1])
      expect_logit_ztest(coef(fit)[2], diss.mle(y0,y1), vcov(fit, sources="estimation")[2,2])

      fit<-tergm(list(y0,y1) ~ edges, estimate="CMLE", control=ctrl, times=c(1,2))
      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], plain.mle(y0,y1), vcov(fit, sources="estimation")[1,1])

      fit<-tergm(list(y0,y1) ~ Cross(~edges), estimate="CMLE", control=ctrl, times=c(1,2))
      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], cross.mle(y0,y1), vcov(fit, sources="estimation")[1,1])

      fit<-tergm(list(y0,y1) ~ Change(~edges), estimate="CMLE", control=ctrl, times=c(1,2))
      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], change.mle(y0,y1), vcov(fit, sources="estimation")[1,1])
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

    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], form.mle(y0,y1m))
    expect_logit_ztest(coef(fit)[2], diss.mle(y0,y1m))

    fit<-tergm(list(y0,y1m) ~ edges, estimate="CMPLE", times=c(1,2))
    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], plain.mle(y0,y1m))

    fit<-tergm(list(y0,y1m) ~ Cross(~edges), estimate="CMPLE", times=c(1,2))
    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], cross.mle(y0,y1m))

    fit<-tergm(list(y0,y1m) ~ Change(~edges), estimate="CMPLE", times=c(1,2))
    expect_equal(fit$estimate, "CMPLE")
    expect_logit_ztest(coef(fit)[1], change.mle(y0,y1m))

  })

  test_that(paste("Autodetect CMPLE on", completeness, netdesc), {

    set.seed(765)
    fit<-tergm(list(y0,y1m) ~ Form(~edges) + Persist(~edges), estimate="CMLE", times=c(1,2))

    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], form.mle(y0,y1m))
    expect_logit_ztest(coef(fit)[2], diss.mle(y0,y1m))

    fit<-tergm(list(y0,y1m) ~ edges, estimate="CMLE", times=c(1,2))
    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], plain.mle(y0,y1m))

    fit<-tergm(list(y0,y1m) ~ Cross(~edges), estimate="CMLE", times=c(1,2))
    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], cross.mle(y0,y1m))

    fit<-tergm(list(y0,y1m) ~ Change(~edges), estimate="CMLE", times=c(1,2))
    expect_equal(fit$estimate, "CMLE")
    expect_logit_ztest(coef(fit)[1], change.mle(y0,y1m))

  })
  
  for(prop.weight in prop.weights){

    test_that(paste("Force CMLE on", completeness, netdesc, "proposal", prop.weight), {

      ctrl <- control.tergm(CMLE.ergm=control.ergm(force.main=TRUE, MCMC.prop.weights=prop.weight))

      set.seed(1234)
      fit<-tergm(list(y0,y1m) ~ Form(~edges) + Persist(~edges), estimate="CMLE", control=ctrl, times=c(1,2))

      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], form.mle(y0,y1m), vcov(fit, sources="estimation")[1,1])
      expect_logit_ztest(coef(fit)[2], diss.mle(y0,y1m), vcov(fit, sources="estimation")[2,2])


      fit<-tergm(list(y0,y1m) ~ edges, estimate="CMLE", control=ctrl, times=c(1,2))
      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], plain.mle(y0,y1m), vcov(fit, sources="estimation")[1,1])

      fit<-tergm(list(y0,y1m) ~ Cross(~edges), estimate="CMLE", control=ctrl, times=c(1,2))
      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], cross.mle(y0,y1m), vcov(fit, sources="estimation")[1,1])

      fit<-tergm(list(y0,y1m) ~ Change(~edges), estimate="CMLE", control=ctrl, times=c(1,2))
      expect_equal(fit$estimate, "CMLE")
      expect_logit_ztest(coef(fit)[1], change.mle(y0,y1m), vcov(fit, sources="estimation")[1,1])

    })

  }
}
