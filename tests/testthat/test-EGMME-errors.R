#  File tests/testthat/test-EGMME-errors.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################

test_that("EGMME errors when it should", {
  nw <- network.initialize(100, directed = FALSE)
  nw %v% "attr" <- rep(1:4, length.out = 100)
  nw <- simulate(nw ~ Form(~edges) + Diss(~edges), coef = c(-5, 1), dynamic = TRUE, output = "final", time.slices = 10)

  expect_error(tergm(nw ~ Form(~edges) + Diss(~edges), 
                     targets = ~edges + mean.age, 
                     target.stats = c(120.9842, 3.718282), 
                     SAN.offsets = c(1),
                     estimate = "EGMME"), "Incorrect number")

  expect_error(tergm(nw ~ Form(~offset(edges)) + Diss(~edges), 
                     targets = ~offset(edges) + mean.age, 
                     target.stats = c(3.718282), 
                     estimate = "EGMME"), "Incorrect number")

  expect_error(tergm(nw ~ Form(~edges + nodefactor("attr", levels = c(2,4))) + Diss(~edges), 
                     targets = ~edges + offset(nodefactor("attr", levels = 2:4), 2) + mean.age, 
                     target.stats = c(120.9842, 50, 70, 3.718282), 
                     SAN.offsets = c(1),
                     estimate = "EGMME"), "Failed to remove")

  expect_error(tergm(nw ~ Form(~offset(edges)) + Diss(~edges), 
                     targets = ~Form(~offset(edges)) + edges, 
                     target.stats = c(120.9842), 
                     estimate = "EGMME"), "Failed to remove")

  expect_error(tergm(nw ~ Form(~edges) + Diss(~edges) + triangle, 
                     targets = ~edges + triangle + mean.age, 
                     target.stats = c(120.9842, 500, 100), 
                     estimate = "EGMME"), "No initial parameter method")

  expect_error(tergm(nw ~ Form(~edges) + Diss(~edges) + offset(triangle), 
                     targets = ~edges + mean.age, 
                     target.stats = c(120.9842, 100), 
                     offset.coef = c(1),
                     estimate = "EGMME"), "No initial parameter method")
})
