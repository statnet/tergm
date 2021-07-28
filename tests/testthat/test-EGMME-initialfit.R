#  File tests/testthat/test-EGMME-initialfit.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2021 Statnet Commons
################################################################################

o <- options(ergm.eval.loglik=FALSE)

test_that("EGMME initialfit behaves reasonably", {

  test.EGMME.initialfit <- function(nw, ff, mff, target.stats, init, control = control.tergm()) {
    model <- ergm_model(ff, nw = nw, dynamic = TRUE)
    if(missing(init)) init <- rep(NA, nparam(model, canonical = FALSE))
    model.mon <- ergm_model(mff, nw = nw, dynamic = TRUE)
    model.mon$target.stats <- target.stats
    names(model.mon$target.stats) <- param_names(model.mon, canonical = FALSE)
  
    tergm.EGMME.initialfit(init, nw, model, ff, model.mon, mff, control)
  }
  
  nw <- network.initialize(100, dir = FALSE)
  nw %v% "attr" <- rep(1:5, length.out = 100)
  nw %v% "cov" <- runif(100)
  
  
  
  ff <- ~Form(~edges) + Diss(~edges)
  mff <- ~edges + mean.age
  target.stats <- c(100, 10)
  
  set.seed(0)
  x <- test.EGMME.initialfit(nw, ff, mff, target.stats)
  set.seed(0)
  e <- ergm(nw ~ edges, target.stats = target.stats[1])
  expected <- unname(c(coef(e)[1] - log(11), -log(11)))
  expect_equal(expected, coef(x))
  
  
  
  ff <- ~Form(~edges) + Persist(~edges)
  mff <- ~edges + mean.age
  target.stats <- c(100, 10)
  
  set.seed(0)
  x <- test.EGMME.initialfit(nw, ff, mff, target.stats)
  set.seed(0)
  e <- ergm(nw ~ edges, target.stats = target.stats[1])
  expected <- unname(c(coef(e)[1] - log(11), log(11)))
  expect_equal(expected, coef(x))
  
  
  
  ff <- ~Form(~edges + nodecov("cov") + offset(gwesp(0, fixed = TRUE))) + Persist(~edges)
  mff <- ~edges + nodecov("cov") + mean.age
  target.stats <- c(100, 100, 10)
  init <- c(NA,NA,0.1,NA)
  
  set.seed(0)
  x <- test.EGMME.initialfit(nw, ff, mff, target.stats, init)
  set.seed(0)
  e <- ergm(nw ~ edges + nodecov("cov") + offset(gwesp(0, fixed = TRUE)), target.stats = target.stats[1:2], offset.coef = 0.1)
  expected <- unname(c(coef(e)[1] - log(11), coef(e)[2], 0.1, log(11)))
  expect_equal(expected, coef(x))
  
  
  
  ff <- ~Diss(~edges) + Form(~edges + nodecov("cov") + offset(gwesp(0, fixed = TRUE)))
  mff <- ~edges + nodecov("cov") + mean.age
  target.stats <- c(100, 100, 10)
  init <- c(NA,NA,NA,0.1)
  
  set.seed(0)
  x <- test.EGMME.initialfit(nw, ff, mff, target.stats, init)
  set.seed(0)
  e <- ergm(nw ~ edges + nodecov("cov") + offset(gwesp(0, fixed = TRUE)), target.stats = target.stats[1:2], offset.coef = 0.1)
  expected <- unname(c(-log(11), coef(e)[1] - log(11), coef(e)[2], 0.1))
  expect_equal(expected, coef(x))
  
  
  
  ff <- ~Diss(~offset(edges)) + Form(~edges + nodecov("cov")) + offset(Persist(~nodecov("cov"))) + Form(~offset(gwesp(0, fixed = TRUE)))
  mff <- ~edges + nodecov("cov")
  target.stats <- c(100, 100)
  init <- c(0.2,NA,NA,0.3,0.1)
  
  set.seed(0)
  x <- test.EGMME.initialfit(nw, ff, mff, target.stats, init)
  set.seed(0)
  e <- ergm(nw ~ edges + nodecov("cov") + offset(gwesp(0, fixed = TRUE)), target.stats = target.stats[1:2], offset.coef = 0.1)
  expected <- unname(c(0.2, coef(e)[1] + 0.2, coef(e)[2] - 0.3, 0.3, 0.1))
  expect_equal(expected, coef(x))
  
  
  
  ff <- ~Diss(~concurrent) + Form(~edges + nodecov("cov")) + offset(Persist(~triangle)) + Form(~gwesp(0, fixed = TRUE))
  mff <- ~edges + nodecov("cov") + concurrent + triangle
  target.stats <- c(100, 100, 100, 100)
  init <- c(0.2,0.8,7.3,0.3,0.1)
  
  set.seed(0)
  x <- test.EGMME.initialfit(nw, ff, mff, target.stats, init)
  expect_equal(init, coef(x))
  
  
  
  ff <- ~Diss(~concurrent) + Form(~edges + nodecov("cov")) + offset(Persist(~triangle)) + Form(~gwesp(0, fixed = TRUE)) + Form(~degree(1))
  mff <- ~edges + nodecov("cov") + concurrent + triangle + degree(1)
  target.stats <- c(100, 100, 100, 100, 100)
  init <- c(0.2,0.8,NA,0.3,NA,0.5)
  
  control <- control.tergm(init.method = "zeros")
  
  set.seed(0)
  x <- test.EGMME.initialfit(nw, ff, mff, target.stats, init, control)
  expect_equal(c(0.2,0.8,0,0.3,0,0.5), coef(x))
  
  
  
  ## try some errors
  ff <- ~Diss(~offset(concurrent)) + Form(~edges + nodecov("cov")) + offset(Persist(~nodecov("cov"))) + Form(~offset(gwesp(0, fixed = TRUE)))
  mff <- ~edges + nodecov("cov")
  target.stats <- c(100, 100)
  init <- c(0.2,NA,NA,0.3,0.1)
  
  expect_error(test.EGMME.initialfit(nw, ff, mff, target.stats, init), "No initial parameter method")
  
  
  
  ff <- ~Diss(~edges) + Form(~edges + nodecov("cov")) + offset(Persist(~nodecov("cov"))) + Form(~offset(gwesp(0, fixed = TRUE)))
  mff <- ~edges + nodecov("cov") + concurrent
  target.stats <- c(100, 100, 80)
  init <- c(NA,NA,NA,0.3,0.1)
  
  expect_error(test.EGMME.initialfit(nw, ff, mff, target.stats, init), "No initial parameter method")
  
  
  
  ff <- ~Form(~edges) + Diss(~edges) + Persist(~edges)
  mff <- ~edges + concurrent + mean.age
  target.stats <- c(100, 50, 10)
  
  expect_error(test.EGMME.initialfit(nw, ff, mff, target.stats), "No initial parameter method")
})

options(o)
