#  File tests/testthat/test-term-Change.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2022 Statnet Commons
################################################################################

test_that("term Change behaves reasonably in dynamic contexts", {
  nw0 <- network.initialize(100, dir = FALSE)
  set.seed(0)
  nw1 <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-4, 2), dynamic = TRUE, output = "final", time.slices = 1)
  set.seed(0)
  nw <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-4, 2), dynamic = TRUE, output = "final", time.slices = 10)
  set.seed(0)
  nwd <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-4, 2), dynamic = TRUE, output = "networkDynamic", time.slices = 10)

  nw10 <- network.collapse(nwd, at = 10)
  nw9 <- network.collapse(nwd, at = 9)
  
  nwdiff <- (nw10 - nw9) + (nw9 - nw10)
  nwdiff %n% "time" <- 10
  nwdiff %n% "lasttoggle" <- cbind(as.edgelist(nwdiff), 10)

  nwdiff01 <- (nw1 - nw0) + (nw0 - nw1)
  nwdiff01 %n% "time" <- 1
  nwdiff01 %n% "lasttoggle" <- cbind(as.edgelist(nwdiff01), 1)

  
  ## non-durational, non-curved
  ff0 <- ~edges + triangle + gwesp(0, fixed = TRUE) + isolates + concurrent
  s1 <- summary(ff0, basis = nwdiff)
  s2 <- summary(~Change(ff0), basis = nw)
  names(s1) <- paste0("Change~", names(s1))
  expect_equal(s1, s2)
  
  s1 <- summary(ff0, basis = nwdiff01)
  s2 <- summary(~Change(ff0), basis = nw1)
  names(s1) <- paste0("Change~", names(s1))
  expect_equal(s1, s2)

  set.seed(0)
  s3 <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-4, 2), dynamic = TRUE, output = "stats", monitor = ~Change(ff0), time.slices = 1)
  expect_equal(s1, s3[1,])

  ## non-durational, curved
  ff1 <- ~edges + triangle + gwesp(0, fixed = TRUE) + gwesp(fixed = FALSE, cutoff = 5) + isolates
  s1 <- summary(ff1, basis = nwdiff)
  s2 <- summary(~Change(ff1), basis = nw)
  names(s1) <- paste0("Change~", names(s1))
  expect_equal(s1, s2)

  s1 <- summary(ff1, basis = nwdiff01)
  s2 <- summary(~Change(ff1), basis = nw1)
  names(s1) <- paste0("Change~", names(s1))
  expect_equal(s1, s2)
  
  set.seed(0)
  s3 <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-4, 2), dynamic = TRUE, output = "stats", monitor = ~Change(ff1), time.slices = 1)
  expect_equal(s1, s3[1,])
  
  ## durational, non-curved
  ff2 <- ~edges + triangle + gwesp(0, fixed = TRUE) + isolates + concurrent + edges.ageinterval(1) + edges.ageinterval(4) + edges.ageinterval(15) + mean.age + Form(~edges + triangle + mean.age) + Persist(~edges + concurrent + isolates + edges.ageinterval(2) + edges.ageinterval(6))
  s1 <- summary(ff2, basis = nwdiff)
  s2 <- summary(~Change(ff2), basis = nw)
  names(s1) <- paste0("Change~", names(s1))
  expect_equal(s1, s2)

  s1 <- summary(ff2, basis = nwdiff01)
  s2 <- summary(~Change(ff2), basis = nw1)
  names(s1) <- paste0("Change~", names(s1))
  expect_equal(s1, s2)

  set.seed(0)
  s3 <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-4, 2), dynamic = TRUE, output = "stats", monitor = ~Change(ff2), time.slices = 1)
  expect_equal(s1, s3[1,])

  ## durational, curved
  ff3 <- ~edges + triangle + gwesp(0, fixed = TRUE) + isolates + gwesp(fixed = FALSE, cutoff = 5) + mean.age + Form(~edges + triangle + mean.age) + Diss(~gwesp(fixed = FALSE, cutoff = 5))
  s1 <- summary(ff3, basis = nwdiff)
  s2 <- summary(~Change(ff3), basis = nw)
  names(s1) <- paste0("Change~", names(s1))
  expect_equal(s1, s2)

  s1 <- summary(ff3, basis = nwdiff01)
  s2 <- summary(~Change(ff3), basis = nw1)
  names(s1) <- paste0("Change~", names(s1))
  expect_equal(s1, s2)

  set.seed(0)
  s3 <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-4, 2), dynamic = TRUE, output = "stats", monitor = ~Change(ff3), time.slices = 1)
  expect_equal(s1, s3[1,])

  set.seed(0)
  sim01 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-4, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ff0,
                    output = "networkDynamic")

  set.seed(0)
  sim02 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-4, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~Change(ff0),
                    output = "stats")
                    
  for(time_step in 10:19) {
    nw_initial <- network.collapse(sim01, at = time_step)
    nw_final <- network.collapse(sim01, at = time_step + 1)
    nw_diff <- (nw_initial - nw_final) + (nw_final - nw_initial)
    
    summ_stats <- summary(ff0, basis = nw_diff)
    names(summ_stats) <- paste0("Change~", names(summ_stats))
    expect_equal(summ_stats, sim02[time_step - 10 + 1,])
  }
                    
  set.seed(1)
  sim11 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-4, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ff1,
                    output = "networkDynamic")

  set.seed(1)
  sim12 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-4, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~Change(ff1),
                    output = "stats")

  for(time_step in 10:19) {
    nw_initial <- network.collapse(sim11, at = time_step)
    nw_final <- network.collapse(sim11, at = time_step + 1)
    nw_diff <- (nw_initial - nw_final) + (nw_final - nw_initial)
    
    summ_stats <- summary(ff1, basis = nw_diff)
    names(summ_stats) <- paste0("Change~", names(summ_stats))
    expect_equal(summ_stats, sim12[time_step - 10 + 1,])
  }
  
  m01 <- ergm_model(ff0, nw = nw)
  m02 <- ergm_model(~Change(ff0), nw = nw)
  expect_true(!is.curved(m01))
  expect_true(!is.curved(m02))
  expect_true(!is.durational(m01))
  expect_true(is.durational(m02))
  expect_identical(paste0("Change~", param_names(m01, canonical = FALSE)), param_names(m02, canonical = FALSE))
  expect_identical(paste0("Change~", param_names(m01, canonical = TRUE)), param_names(m02, canonical = TRUE))

  m11 <- ergm_model(ff1, nw = nw)
  m12 <- ergm_model(~Change(ff1), nw = nw)
  expect_true(is.curved(m11))
  expect_true(is.curved(m12))
  expect_true(!is.durational(m11))
  expect_true(is.durational(m12))
  expect_identical(paste0("Change~", param_names(m11, canonical = FALSE)), param_names(m12, canonical = FALSE))
  expect_identical(paste0("Change~", param_names(m11, canonical = TRUE)), param_names(m12, canonical = TRUE))

  m21 <- ergm_model(ff2, nw = nw)
  m22 <- ergm_model(~Change(ff2), nw = nw)
  expect_true(!is.curved(m21))
  expect_true(!is.curved(m22))
  expect_true(is.durational(m21))
  expect_true(is.durational(m22))
  expect_identical(paste0("Change~", param_names(m21, canonical = FALSE)), param_names(m22, canonical = FALSE))
  expect_identical(paste0("Change~", param_names(m21, canonical = TRUE)), param_names(m22, canonical = TRUE))

  m31 <- ergm_model(ff3, nw = nw)
  m32 <- ergm_model(~Change(ff3), nw = nw)
  expect_true(is.curved(m31))
  expect_true(is.curved(m32))
  expect_true(is.durational(m31))
  expect_true(is.durational(m32))
  expect_identical(paste0("Change~", param_names(m31, canonical = FALSE)), param_names(m32, canonical = FALSE))
  expect_identical(paste0("Change~", param_names(m31, canonical = TRUE)), param_names(m32, canonical = TRUE))
})
