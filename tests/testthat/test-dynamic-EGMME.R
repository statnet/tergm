#  File tests/testthat/test-dynamic-EGMME.R in package tergm, part of the
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
  test_that("dynamic EGMME fits correctly", {
    n <- 40
    do.plot <- FALSE
    g0 <- network.initialize(n, dir = FALSE)

    #                    edges, degree(1), mean.age
    target.stats <- c(n * 1 / 2, n * 0.6, 20)

    # Get a deliberately bad starting network.
    set.seed(1)
    g1 <- san(
      g0 ~ meandeg + degree(1),
      target.stats = target.stats[-3],
      verbose = TRUE
    )

    coef.form <- c(-6.57, 1.01)
    coef.diss <- c(2.944439)

    # Fit the model with very poor starting values.
    set.seed(3)
    dynfit <- tergm(
      g1 ~ Form(~edges + degree(1)) + Diss(~offset(edges)),
      targets = "formation",
      estimate = "EGMME",
      constraints = ~.,
      offset.coef = -coef.diss,
      target.stats = target.stats[-3],
      verbose = TRUE,
      control = control.tergm(
        SA.plot.progress = do.plot,
        SA.restart.on.err = FALSE,
        init = c(-log(.95 / .05), 0, -coef.diss)
      )
    )

    expect_warning(expect_error(print(summary(dynfit)), NA), NA)
    expect_warning(expect_error(mcmc.diagnostics(dynfit), NA), NA)

    expect_equal(
      c(coef.form, -coef.diss),
      coef(dynfit),
      tolerance = 0.01,
      ignore_attr = TRUE
    )

    # All parameters free, edges, degree(1), and edge.ages as target.
    set.seed(5)
    dynfit2 <- tergm(
      g1 ~ Form(~edges + degree(1)) + Diss(~edges),
      targets = ~edges + degree(1) + mean.age,
      estimate = "EGMME",
      constraints = ~.,
      target.stats = target.stats,
      control = control.tergm(
        SA.plot.progress = do.plot,
        SA.plot.stats = TRUE
      )
    )

    expect_warning(expect_error(print(summary(dynfit2)), NA), NA)
    expect_warning(expect_error(mcmc.diagnostics(dynfit2), NA), NA)

    expect_equal(
      c(coef.form, -coef.diss),
      coef(dynfit2),
      tolerance = 0.01,
      ignore_attr = TRUE
    )
  })
}, "EGMME")
