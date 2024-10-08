#  File tests/testthat/test-term-.SubsetStatistics.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2024 Statnet Commons
################################################################################

test_that("the .SubsetStatistics term behaves appropriately", {
  nw0 <- network.initialize(100, dir = FALSE)
  nw <- simulate(nw0 ~ edges, coef = c(-4), dynamic = TRUE, output = "final")

  ## non-durational, non-curved
  ff0 <- ~edges + triangle + gwesp(0, fixed = TRUE) + isolates + concurrent
  stats0 <- c(TRUE, FALSE, TRUE, TRUE, FALSE)
  s1 <- summary(ff0, basis = nw)[stats0]
  s2 <- summary(~.SubsetStatistics(ff0, stats0), basis = nw)
  expect_equal(s1, s2)

  ## non-durational, curved
  ff1 <- ~edges + triangle + gwesp(0, fixed = TRUE) + gwesp(fixed = FALSE, cutoff = 100) + isolates
  stats1 <- sample(which(c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE)))
  s1 <- summary(ff1, basis = nw)[stats1]
  s2 <- summary(~.SubsetStatistics(ff1, stats1), basis = nw)
  expect_equal(s1, s2)

  set.seed(0)
  nw1 <- simulate(nw ~ Form(~edges + triangle) + Persist(~edges + gwesp(0, fixed = TRUE)), coef = c(-6, 0, 4, 0), time.slices = 10, dynamic = TRUE, output = "final")
  set.seed(0)
  nw2 <- simulate(nw ~ .SubsetStatistics(~Form(~edges + triangle) + Persist(~edges + gwesp(0, fixed = TRUE)), c(3,1)), coef = c(4, -6), time.slices = 10, dynamic = TRUE, output = "final")
  attr(nw1, "coef") <- NULL
  attr(nw2, "coef") <- NULL
  expect_equal(nw1, nw2)
  nw <- nw1

  ## durational, non-curved
  ff2 <- ~edges + triangle + gwesp(0, fixed = TRUE) + isolates + concurrent + edges.ageinterval(1) + edges.ageinterval(4) + edges.ageinterval(15) + mean.age + Form(~edges + triangle + mean.age) + Persist(~edges + concurrent + isolates + edges.ageinterval(2) + edges.ageinterval(6))
  stats2 <- sample(which(c(TRUE, FALSE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)))
  s1 <- summary(ff2, basis = nw)[stats2]
  s2 <- summary(~.SubsetStatistics(ff2, I(c(stats2, 27L))), basis = nw)
  expect_equal(s1, s2)

  ## durational, curved
  ff3 <- ~edges + triangle + gwesp(0, fixed = TRUE) + isolates + gwesp(fixed = FALSE, cutoff = 100) + mean.age + Form(~edges + triangle + mean.age) + Diss(~gwesp(fixed = FALSE, cutoff = 100))
  stats3 <- c(TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, TRUE, TRUE)
  s1 <- summary(ff3, basis = nw)[stats3]
  s2 <- summary(~.SubsetStatistics(ff3, stats3), basis = nw)
  expect_equal(s1, s2)

  set.seed(0)
  sim01 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ff0,
                    output = "stats")[,stats0]

  set.seed(0)
  sim02 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~.SubsetStatistics(ff0,stats0),
                    output = "stats")

  expect_equal(sim01, sim02)

  set.seed(1)
  sim11 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ff1,
                    output = "stats")[,stats1]

  set.seed(1)
  sim12 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~.SubsetStatistics(ff1,stats1),
                    output = "stats")

  expect_equal(sim11, sim12)

  set.seed(2)
  sim21 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ff2,
                    output = "stats")[,stats2]

  set.seed(2)
  sim22 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~.SubsetStatistics(ff2,I(c(stats2, 27L))),
                    output = "stats")

  expect_equal(sim21, sim22)

  set.seed(3)
  sim31 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ff3,
                    output = "stats")[,stats3]

  set.seed(3)
  sim32 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~.SubsetStatistics(ff3,stats3),
                    output = "stats")

  expect_equal(sim31, sim32)

  m01 <- ergm_model(ff0, nw = nw)
  m02 <- ergm_model(~.SubsetStatistics(ff0, stats0), nw = nw)
  expect_true(!is.curved(m01))
  expect_true(!is.curved(m02))
  expect_true(!is.durational(m01))
  expect_true(!is.durational(m02))
  expect_identical(param_names(m01, canonical = FALSE)[stats0], param_names(m02, canonical = FALSE))
  expect_identical(param_names(m01, canonical = TRUE)[stats0], param_names(m02, canonical = TRUE))

  m11 <- ergm_model(ff1, nw = nw)
  m12 <- ergm_model(~.SubsetStatistics(ff1, stats1), nw = nw)
  expect_true(is.curved(m11))
  expect_true(is.curved(m12))
  expect_true(!is.durational(m11))
  expect_true(!is.durational(m12))
  expect_identical(param_names(m11, canonical = FALSE), param_names(m12, canonical = FALSE))
  expect_identical(param_names(m11, canonical = TRUE)[stats1], param_names(m12, canonical = TRUE))

  m21 <- ergm_model(ff2, nw = nw)
  m22 <- ergm_model(~.SubsetStatistics(ff2, I(c(stats2, 27L))), nw = nw)
  expect_true(!is.curved(m21))
  expect_true(!is.curved(m22))
  expect_true(is.durational(m21))
  expect_true(is.durational(m22))
  expect_identical(param_names(m21, canonical = FALSE)[stats2], param_names(m22, canonical = FALSE))
  expect_identical(param_names(m21, canonical = TRUE)[stats2], param_names(m22, canonical = TRUE))

  m31 <- ergm_model(ff3, nw = nw)
  m32 <- ergm_model(~.SubsetStatistics(ff3, stats3), nw = nw)
  expect_true(is.curved(m31))
  expect_true(is.curved(m32))
  expect_true(is.durational(m31))
  expect_true(is.durational(m32))
  expect_identical(param_names(m31, canonical = FALSE), param_names(m32, canonical = FALSE))
  expect_identical(param_names(m31, canonical = TRUE)[stats3], param_names(m32, canonical = TRUE))
})
