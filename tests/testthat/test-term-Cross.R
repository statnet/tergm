#  File tests/testthat/test-term-Cross.R in package tergm, part of the Statnet
#  suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################

test_that("term Cross behaves reasonably in dynamic contexts", {
  nw0 <- network.initialize(100, dir = FALSE)
  nw <- simulate(nw0 ~ edges, coef = c(-4), dynamic = TRUE, output = "final")

  ## non-durational, non-curved
  ff0 <- ~edges + triangle + gwesp(0, fixed = TRUE) + isolates + concurrent
  s1 <- summary(ff0, basis = nw)
  s2 <- summary(~Cross(ff0), basis = nw)
  names(s1) <- paste0("Cross~", names(s1))
  expect_equal(s1, s2)

  ## non-durational, curved
  ff1 <- ~edges + triangle + gwesp(0, fixed = TRUE) + gwesp(fixed = FALSE, cutoff = 100) + isolates
  s1 <- summary(ff1, basis = nw)
  s2 <- summary(~Cross(ff1), basis = nw)
  names(s1) <- paste0("Cross~", names(s1))
  expect_equal(s1, s2)

  set.seed(0)
  nw1 <- simulate(nw ~ Form(~edges + triangle) + Persist(~edges + gwesp(0, fixed = TRUE)), coef = c(-6, 0, 4, 0), time.slices = 10, dynamic = TRUE, output = "final")
  set.seed(0)
  nw2 <- simulate(nw ~ Cross(~Form(~edges + triangle) + Persist(~edges + gwesp(0, fixed = TRUE))), coef = c(-6, 0, 4, 0), time.slices = 10, dynamic = TRUE, output = "final")
  expect_equal(nw1, nw2)
  nw <- nw1

  ## durational, non-curved
  ff2 <- ~edges + triangle + gwesp(0, fixed = TRUE) + isolates + concurrent + edges.ageinterval(1) + edges.ageinterval(4) + edges.ageinterval(15) + mean.age + Form(~edges + triangle + mean.age) + Persist(~edges + concurrent + isolates + edges.ageinterval(2) + edges.ageinterval(6))
  s1 <- summary(ff2, basis = nw)
  s2 <- summary(~Cross(ff2), basis = nw)
  names(s1) <- paste0("Cross~", names(s1))
  expect_equal(s1, s2)

  ## durational, curved
  ff3 <- ~edges + triangle + gwesp(0, fixed = TRUE) + isolates + gwesp(fixed = FALSE, cutoff = 100) + mean.age + Form(~edges + triangle + mean.age) + Diss(~gwesp(fixed = FALSE, cutoff = 100))
  s1 <- summary(ff3, basis = nw)
  s2 <- summary(~Cross(ff3), basis = nw)
  names(s1) <- paste0("Cross~", names(s1))
  expect_equal(s1, s2)

  set.seed(0)
  sim01 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ff0,
                    output = "stats")

  set.seed(0)
  sim02 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~Cross(ff0),
                    output = "stats")
                    
  colnames(sim01) <- paste0("Cross~", colnames(sim01))
  expect_equal(sim01, sim02)

  set.seed(1)
  sim11 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ff1,
                    output = "stats")

  set.seed(1)
  sim12 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~Cross(ff1),
                    output = "stats")

  colnames(sim11) <- paste0("Cross~", colnames(sim11))
  expect_equal(sim11, sim12)

  set.seed(2)
  sim21 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ff2,
                    output = "stats")

  set.seed(2)
  sim22 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~Cross(ff2),
                    output = "stats")

  colnames(sim21) <- paste0("Cross~", colnames(sim21))
  expect_equal(sim21, sim22)

  set.seed(3)
  sim31 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ff3,
                    output = "stats")

  set.seed(3)
  sim32 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-6, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~Cross(ff3),
                    output = "stats")

  colnames(sim31) <- paste0("Cross~", colnames(sim31))
  expect_equal(sim31, sim32)

  m01 <- ergm_model(ff0, nw = nw)
  m02 <- ergm_model(~Cross(ff0), nw = nw)
  expect_true(!is.curved(m01))
  expect_true(!is.curved(m02))
  expect_true(!is.durational(m01))
  expect_true(!is.durational(m02))
  expect_identical(paste0("Cross~", param_names(m01, canonical = FALSE)), param_names(m02, canonical = FALSE))
  expect_identical(paste0("Cross~", param_names(m01, canonical = TRUE)), param_names(m02, canonical = TRUE))

  m11 <- ergm_model(ff1, nw = nw)
  m12 <- ergm_model(~Cross(ff1), nw = nw)
  expect_true(is.curved(m11))
  expect_true(is.curved(m12))
  expect_true(!is.durational(m11))
  expect_true(!is.durational(m12))
  expect_identical(paste0("Cross~", param_names(m11, canonical = FALSE)), param_names(m12, canonical = FALSE))
  expect_identical(paste0("Cross~", param_names(m11, canonical = TRUE)), param_names(m12, canonical = TRUE))

  m21 <- ergm_model(ff2, nw = nw)
  m22 <- ergm_model(~Cross(ff2), nw = nw)
  expect_true(!is.curved(m21))
  expect_true(!is.curved(m22))
  expect_true(is.durational(m21))
  expect_true(is.durational(m22))
  expect_identical(paste0("Cross~", param_names(m21, canonical = FALSE)), param_names(m22, canonical = FALSE))
  expect_identical(paste0("Cross~", param_names(m21, canonical = TRUE)), param_names(m22, canonical = TRUE))

  m31 <- ergm_model(ff3, nw = nw)
  m32 <- ergm_model(~Cross(ff3), nw = nw)
  expect_true(is.curved(m31))
  expect_true(is.curved(m32))
  expect_true(is.durational(m31))
  expect_true(is.durational(m32))
  expect_identical(paste0("Cross~", param_names(m31, canonical = FALSE)), param_names(m32, canonical = FALSE))
  expect_identical(paste0("Cross~", param_names(m31, canonical = TRUE)), param_names(m32, canonical = TRUE))
})
