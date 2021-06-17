#  File tests/testthat/test-stergm-init.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################

test_that("stergm correctly handles initial coefficient specifications", {
  set.seed(2)
  nw <- network.initialize(100, directed = FALSE)
  nw %v% "attr" <- rep(1:5, length.out = 100)
  nw %v% "cov" <- runif(100)
  nw <- simulate(nw ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "final", time.slices = 50)
  nw1 <- simulate(nw ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "final", time.slices = 1)
  nw2 <- simulate(nw1 ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "final", time.slices = 1)
  nw3 <- simulate(nw2 ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "final", time.slices = 1)
  nwlist <- list(nw1, nw2, nw3)
  
  ## valid specification tests
  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.form = c(1,1),
               offset.coef.diss = c(-1),
               estimate = "CMLE",
               control = control.stergm())

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,-1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.form = c(1,1),
               estimate = "CMLE",
               control = control.stergm(init.diss = c(NA, -1, NA, NA, NA, NA, NA)))

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,-1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.diss = c(-1),
               estimate = "CMLE",
               control = control.stergm(init.form = c(NA, NA, 1, NA, 1)))

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,-1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               estimate = "CMLE",
               control = control.stergm(init.form = c(NA, NA, 1, NA, 1),
                                        init.diss = c(NA, -1, NA, NA, NA, NA, NA)))

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,-1), check.attributes = FALSE)

  ## override tests
  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.form = c(1,1),
               offset.coef.diss = c(2),
               estimate = "CMLE",
               control = control.stergm(init.diss = c(NA, -1, NA, NA, NA, NA, NA)))

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,2), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.form = c(5,-3),
               offset.coef.diss = c(-1),
               estimate = "CMLE",
               control = control.stergm(init.form = c(NA, NA, 1, NA, 1)))

  expect_equal(coef(rv)[c(3,5,7)], c(5,-3,-1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               estimate = "CMLE",
               offset.coef.form = c(5,-3),
               offset.coef.diss = c(1),               
               control = control.stergm(init.form = c(NA, NA, 1, NA, 1),
                                        init.diss = c(NA, -1, NA, NA, NA, NA, NA)))

  expect_equal(coef(rv)[c(3,5,7)], c(5,-3,1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + nodecov("cov") + nodematch("attr", diff = TRUE),
               offset.coef.form = c(4,-7),
               estimate = "CMLE",
               control = control.stergm(init.form = c(NA, NA, 1, NA, 1)))

  expect_equal(coef(rv)[c(3,5)], c(4,-7), check.attributes = FALSE)  
})
