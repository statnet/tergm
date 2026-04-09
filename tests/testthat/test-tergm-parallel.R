#  File tests/testthat/test-tergm-parallel.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################

library(statnet.common)
opttest({
  test_that("tergm EGMME runs in parallel (PSOCK)", {
    data(florentine)
    net <- flobusiness
    set.seed(1)

    mod1 <- tergm(
      flobusiness ~ Form(~edges + degree(3)) + Persist(~offset(edges)),
      offset.coef = log(9),
      targets = "formation",
      estimate = "EGMME",
      control = control.tergm(parallel = 2, parallel.type = "PSOCK")
    )
    expect_true(inherits(mod1, "tergm"))
  })
}, testname = "tergm_parallel")

opttest({
  test_that("tergm EGMME runs in parallel (MPI)", {
    data(florentine)
    net <- flobusiness
    set.seed(1)

    mod1 <- tergm(
      flobusiness ~ Form(~edges + degree(3)) + Persist(~offset(edges)),
      offset.coef = log(9),
      targets = "formation",
      estimate = "EGMME",
      control = control.tergm(parallel = 2, parallel.type = "MPI")
    )
    expect_true(inherits(mod1, "tergm"))
  })
}, testname = "tergm_parallel_MPI", testvar = "ENABLE_MPI_TESTS")
