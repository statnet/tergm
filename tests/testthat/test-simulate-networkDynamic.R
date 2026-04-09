#  File tests/testthat/test-simulate-networkDynamic.R in package tergm, part of
#  the Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################

test_that("simulate.networkDynamic encodes net.obs.period correctly", {
  set.seed(0)

  logit <- function(p) log(p / (1 - p))
  coef.form.f <- function(coef.diss, density) {
    -log(((1 + exp(coef.diss)) / (density / (1 - density))) - 1)
  }

  # Construct a network with 20 nodes and 20 edges
  n <- 20
  target.stats <- edges <- 20
  g0 <- network.initialize(n, dir = TRUE)
  g1 <- san(g0 ~ edges, target.stats = target.stats, verbose = TRUE)

  S <- 10

  # To get an average duration of 10...
  duration <- 10
  coef.diss <- logit(1 - 1 / duration)

  # To get an average of 20 edges...
  dyads <- network.dyadcount(g1)
  density <- edges / dyads
  coef.form <- coef.form.f(coef.diss, density)

  # Simulate a networkDynamic for five steps
  dynsim <- simulate(
    g1 ~ Form(~edges) + Persist(~edges),
    coef = c(coef.form, coef.diss),
    time.slices = S,
    dynamic = TRUE
  )

  # should have a net.obs.period object from 0 to 5
  expect_true(
    all(unlist((dynsim %n% "net.obs.period")$observations) == c(0, 1, 1, 11)),
    label = "simulate.network net.obs.period$observation encoding"
  )

  # "Resume" the simulation for five steps
  dynsim2 <- simulate(
    dynsim ~ Form(~edges) + Persist(~edges),
    time.slices = S,
    dynamic = TRUE
  )

  expect_true(
    all(
      unlist((dynsim2 %n% "net.obs.period")$observations) ==
        c(0, 1, 1, 11, 11, 21)
    ),
    label = "simulate.networkDynamic net.obs.period$observation encoding"
  )
})

test_that("simulate.networkDynamic encodes net.obs.period correctly when stopping and restarting", {
  logit <- function(p) log(p / (1 - p))
  coef.form.f <- function(coef.diss, density) {
    -log(((1 + exp(coef.diss)) / (density / (1 - density))) - 1)
  }

  # check a realistic example of starting a sim from scratch
  mean.rel.dur <- 10
  msm.sim <- network.initialize(100, directed = FALSE)
  activate.vertices(msm.sim, -Inf, Inf)
  set.network.attribute(
    msm.sim,
    "net.obs.period",
    list(
      observations = list(c(-1, 0)),
      mode = "discrete",
      time.increment = 1,
      time.unit = "step"
    )
  )
  formation <- ~edges
  dissolution <- ~offset(edges)
  target.stats <- 40
  coef.diss <- log(mean.rel.dur - 1)
  formation.with.stnet <- update.formula(formation, msm.startnet ~ .)
  msm.startnet <- network.collapse(msm.sim, at = 0)
  msm.est <- ergm(formation.with.stnet, target.stats = target.stats)
  coef.form <- coef(msm.est)
  coef.form[1] <- coef.form[1] - coef.diss
  msm.edgelist <- as.edgelist(simulate(msm.est, dynamic = FALSE))
  add.edges(msm.sim, msm.edgelist[, 1], msm.edgelist[, 2])
  activate.edges(msm.sim, -Inf, Inf)

  # first timestep
  msm.sim <- simulate(
    msm.sim ~ Form(~edges) + Persist(~edges),
    coef = c(coef.form, coef.diss),
    time.slices = 1,
    monitor = "all",
    verbose = TRUE,
    dynamic = TRUE
  )
  # second timestep
  msm.sim <- simulate(
    msm.sim ~ Form(~edges) + Persist(~edges),
    coef = c(coef.form, coef.diss),
    time.slices = 1,
    monitor = "all",
    verbose = TRUE,
    dynamic = TRUE
  )

  expect_true(
    all(
      unlist((msm.sim %n% "net.obs.period")$observations) ==
        c(-1, 0, 0, 1, 1, 2)
    ),
    label = "net.obs.period when stopping and restarting simulate.networkDynamic"
  )
})

test_that("simulate.networkDynamic handles vertex dynamics correctly", {
  logit <- function(p) log(p / (1 - p))
  coef.form.f <- function(coef.diss, density) {
    -log(((1 + exp(coef.diss)) / (density / (1 - density))) - 1)
  }

  # test a sim with vertex dynamics applied after the sim stage
  mean.rel.dur <- 10
  msm.sim <- network.initialize(100, directed = FALSE)
  activate.vertices(msm.sim, -Inf, Inf)
  set.network.attribute(
    msm.sim,
    "net.obs.period",
    list(
      observations = list(c(-1, 0)),
      mode = "discrete",
      time.increment = 1,
      time.unit = "step"
    )
  )
  formation <- ~edges
  dissolution <- ~offset(edges)
  target.stats <- 40
  coef.diss <- log(mean.rel.dur - 1)
  formation.with.stnet <- update.formula(formation, msm.startnet ~ .)
  # simulate a set of edges to use as the starting point for the network
  msm.startnet <- network.collapse(msm.sim, at = 0)
  msm.est <- ergm(formation.with.stnet, target.stats = target.stats)
  coef.form <- coef(msm.est)
  coef.form[1] <- coef.form[1] - coef.diss
  msm.edgelist <- as.edgelist(simulate(msm.est, dynamic = FALSE))
  add.edges(msm.sim, msm.edgelist[, 1], msm.edgelist[, 2])
  activate.edges(msm.sim, -Inf, Inf)

  # simulate first timestep (0,1)
  msm.sim <- simulate(
    msm.sim ~ Form(~edges) + Persist(~edges),
    coef = c(coef.form, coef.diss),
    time.slices = 1,
    monitor = "formation",
    verbose = TRUE,
    dynamic = TRUE
  )

  # toggle off vertices for the step just simulated
  msm.sim <- deactivate.vertices(
    msm.sim,
    v = sample(
      which(is.active(msm.sim, v = 1:network.size(msm.sim), at = 0)),
      size = 10
    ),
    onset = 0,
    terminus = Inf,
    deactivate.edges = TRUE
  )
  add.vertices.active(msm.sim, nv = 5, onset = 0, terminus = Inf)


  # simulate second timestep (1,2)
  msm.sim <- simulate(
    msm.sim ~ Form(~edges) + Persist(~edges),
    coef = c(coef.form, coef.diss),
    time.slices = 1,
    monitor = "all",
    verbose = TRUE,
    dynamic = TRUE
  )

  msm.sim <- deactivate.vertices(
    msm.sim,
    v = sample(
      which(is.active(msm.sim, v = 1:network.size(msm.sim), at = 1)),
      size = 10
    ),
    onset = 1,
    terminus = Inf,
    deactivate.edges = TRUE
  )
  add.vertices.active(msm.sim, nv = 1, onset = 1, terminus = Inf)

  # check for correct activity
  expect_true(
    all(
      sapply(-1:2, function(t) network.size.active(msm.sim, at = t)) ==
        c(100, 95, 86, 86)
    ),
    label = "vertex dynamics activity in stergm sim"
  )

  expect_true(
    all(network.dynamic.check(msm.sim, complete = FALSE)$dyad.checks),
    label = "vertex dynamics in stergm sim dyad consistency"
  )
})

test_that("simulate.networkDynamic respects time.start and time.offset", {
  logit <- function(p) log(p / (1 - p))
  coef.form.f <- function(coef.diss, density) {
    -log(((1 + exp(coef.diss)) / (density / (1 - density))) - 1)
  }

  mean.rel.dur <- 10
  msm.sim <- network.initialize(100, directed = FALSE)
  activate.vertices(msm.sim, -Inf, Inf)
  set.network.attribute(
    msm.sim,
    "net.obs.period",
    list(
      observations = list(c(-1, 0)),
      mode = "discrete",
      time.increment = 1,
      time.unit = "step"
    )
  )
  formation <- ~edges
  dissolution <- ~offset(edges)
  target.stats <- 40
  coef.diss <- log(mean.rel.dur - 1)
  formation.with.stnet <- update.formula(formation, msm.startnet ~ .)
  msm.startnet <- network.collapse(msm.sim, at = 0)
  msm.est <- ergm(formation.with.stnet, target.stats = target.stats)
  coef.form <- coef(msm.est)
  coef.form[1] <- coef.form[1] - coef.diss
  msm.edgelist <- as.edgelist(simulate(msm.est, dynamic = FALSE))
  add.edges(msm.sim, msm.edgelist[, 1], msm.edgelist[, 2])
  activate.edges(msm.sim, -Inf, Inf)

  # simulate first timestep (0,1)
  msm.sim <- simulate(
    msm.sim ~ Form(~edges) + Persist(~edges),
    coef = c(coef.form, coef.diss),
    time.slices = 1,
    monitor = "formation",
    verbose = TRUE,
    dynamic = TRUE
  )

  # toggle off vertices for the step just simulated
  msm.sim <- deactivate.vertices(
    msm.sim,
    v = sample(
      which(is.active(msm.sim, v = 1:network.size(msm.sim), at = 0)),
      size = 10
    ),
    onset = 0,
    terminus = Inf,
    deactivate.edges = TRUE
  )
  add.vertices.active(msm.sim, nv = 5, onset = 0, terminus = Inf)

  # simulate second timestep (1,2)
  msm.sim <- simulate(
    msm.sim ~ Form(~edges) + Persist(~edges),
    coef = c(coef.form, coef.diss),
    time.slices = 1,
    monitor = "all",
    verbose = TRUE,
    dynamic = TRUE
  )

  msm.sim <- deactivate.vertices(
    msm.sim,
    v = sample(
      which(is.active(msm.sim, v = 1:network.size(msm.sim), at = 1)),
      size = 10
    ),
    onset = 1,
    terminus = Inf,
    deactivate.edges = TRUE
  )
  add.vertices.active(msm.sim, nv = 1, onset = 1, terminus = Inf)

  # at this point msm.sim has observations as far as 2.  Try sampling at 2
  # and updating from (2,3)
  msm.sim <- simulate(
    msm.sim ~ Form(~edges) + Persist(~edges),
    coef = c(coef.form, coef.diss),
    time.slices = 1,
    time.start = 2,
    time.offset = 0,
    monitor = "formation",
    verbose = TRUE,
    dynamic = TRUE
  )

  expect_lte(
    max(unlist((msm.sim %n% "net.obs.period")$observations)),
    3,
    label = "net.obs.period not updated inappropriately when time.start and time.offset was explicitly set"
  )
  expect_true(
    all(get.change.times(msm.sim) == 0:2),
    label = "edge spells not updated inappropriately when time.start and time.offset was explicitly set"
  )

  # try it again to make sure nothing got mangled
  msm.sim <- simulate(
    msm.sim ~ Form(~edges) + Persist(~edges),
    coef = c(coef.form, coef.diss),
    time.slices = 1,
    time.start = 3,
    time.offset = 0,
    monitor = "all",
    verbose = TRUE,
    dynamic = TRUE
  )
  expect_lte(
    max(unlist((msm.sim %n% "net.obs.period")$observations)),
    4,
    label = "net.obs.period not updated inappropriately (second time) when time.start and time.offset was explicitly set"
  )
  expect_true(
    all(get.change.times(msm.sim) == 0:3),
    label = "edge spells not updated inappropriately (second time) when time.start and time.offset was explicitly set"
  )
})

test_that("simulate.networkDynamic does not return changes involving inactive vertices", {
  # check that vertex ids in changes are correctly translated
  dyn <- as.networkDynamic(network.initialize(4))
  deactivate.vertices(dyn, v = 1)
  # define stergm that should toggle on all ties
  changes <- simulate(
    dyn ~ Form(~edges) + Persist(~edges),
    coef = c(1, 0),
    output = "changes",
    dynamic = TRUE
  )
  # check if any changes involve vertex 1 (shouldn't because it is inactive)
  expect_true(
    !any(changes[, 2:3] == 1),
    label = "simulate.networkDynamic changes do not involve inactive vertices"
  )
})

test_that("simulate.networkDynamic handles time.offset=0 correctly", {
  logit <- function(p) log(p / (1 - p))
  coef.form.f <- function(coef.diss, density) {
    -log(((1 + exp(coef.diss)) / (density / (1 - density))) - 1)
  }

  mean.rel.dur <- 10
  msm.sim <- as.networkDynamic(network.initialize(100, directed = FALSE))
  formation <- ~edges
  dissolution <- ~offset(edges)
  target.stats <- 40
  coef.diss <- log(mean.rel.dur - 1)
  formation.with.stnet <- update.formula(formation, msm.startnet ~ .)
  # simulate a set of edges to use as the starting point for the network
  msm.startnet <- network.collapse(msm.sim, at = 0)
  msm.est <- ergm(formation.with.stnet, target.stats = target.stats)
  coef.form <- coef(msm.est)
  coef.form[1] <- coef.form[1] - coef.diss
  msm.edgelist <- as.edgelist(simulate(msm.est, dynamic = FALSE))
  add.edges(msm.sim, msm.edgelist[, 1], msm.edgelist[, 2])

  msm.sim <- simulate(
    msm.sim ~ Form(~edges) + Persist(~edges),
    coef = c(coef.form, coef.diss),
    time.slices = 1,
    time.start = 0,
    time.offset = 0,
    monitor = "all",
    verbose = TRUE,
    dynamic = TRUE
  )

  expect_lte(
    length((msm.sim %n% "net.obs.period")$observations),
    1,
    label = "simulate.networkDynamic net.obs.period not mangled with time.offset=0 argument"
  )
})

test_that("simulate.networkDynamic handles adding and deleting vertices correctly", {
  set.seed(0)
  logit <- function(p) log(p / (1 - p))
  coef.form.f <- function(coef.diss, density) {
    -log(((1 + exp(coef.diss)) / (density / (1 - density))) - 1)
  }

  g0 <- network.initialize(20, dir = TRUE)
  g1 <- san(g0 ~ edges, target.stats = 20, verbose = TRUE)
  # Simulate a networkDynamic from static
  dynsim <- simulate(
    g1 ~ Form(~edges) + Persist(~edges),
    coef = c(-9.326322, 2.197225),
    time.slices = 1,
    verbose = TRUE,
    dynamic = TRUE
  )
  # add more vertices
  add.vertices(dynsim, 3)
  # simulate another step
  dynsim <- simulate(
    dynsim ~ Form(~edges) + Persist(~edges),
    coef = c(-9.326322, 2.197225),
    time.slices = 1,
    verbose = TRUE,
    dynamic = TRUE
  )
  # expect to see a pid defined, but will just be numeric values
  expect_equal(
    dynsim %n% "vertex.pid",
    "tergm_pid",
    label = "simulate.networkDynamic vertex.pid attribute"
  )
  expect_true(
    all(is.numeric(dynsim %v% "tergm_pid")),
    label = "simulate.networkDynamic numeric vertex pid"
  )


  # add more (at this point there should be some non-numeric pids created)
  add.vertices(dynsim, 3)
  expect_false(
    all(is.numeric(dynsim %v% "tergm_pid")),
    label = "add.vertices non-numeric vertex pids"
  )

  # simulate
  dynsim <- simulate(
    dynsim ~ Form(~edges) + Persist(~edges),
    coef = c(-9.326322, 2.197225),
    time.slices = 1,
    verbose = TRUE,
    dynamic = TRUE
  )
  # delete some vertices
  delete.vertices(dynsim, vid = 5:10)
  dynsim <- simulate(
    dynsim ~ Form(~edges) + Persist(~edges),
    coef = c(-9.326322, 2.197225),
    time.slices = 1,
    verbose = TRUE,
    dynamic = TRUE
  )
})

test_that("simulate.networkDynamic respects externally defined vertex PID", {
  set.seed(0)

  g0 <- network.initialize(20, dir = TRUE)
  g1 <- san(g0 ~ edges, target.stats = 20, verbose = TRUE)

  # test with an externally defined PID
  dynsim <- simulate(
    g1 ~ Form(~edges) + Persist(~edges),
    coef = c(-9.326322, 2.197225),
    time.slices = 1,
    verbose = TRUE,
    dynamic = TRUE
  )
  dynsim %n% "vertex.pid" <- "letters"
  dynsim %v% "letters" <- LETTERS[1:20]
  # simulate another step
  dynsim <- simulate(
    dynsim ~ Form(~edges) + Persist(~edges),
    coef = c(-9.326322, 2.197225),
    time.slices = 1,
    verbose = TRUE,
    dynamic = TRUE
  )
  # expect to see a pid defined, but will just be numeric values
  expect_equal(
    dynsim %n% "vertex.pid",
    "letters",
    label = "simulate.networkDynamic vertex pid preserved"
  )
  expect_true(
    all(is.na(dynsim %v% "tergm_pid")),
    label = "simulate.networkDynamic did not create tergm_pid when vertex.pid already existed"
  )
})

test_that("summary.statistics.networkDynamic works on network with active vertices", {
  # this was giving error in #958
  my.nD <- network.initialize(10, directed = FALSE)
  activate.vertices(my.nD, onset = 0, terminus = 10)
  expect_error(summary(my.nD ~ isolates, at = 1), NA)
})

test_that("simulate.networkDynamic handles Inf coefficients correctly", {
  # was giving error in #995
  n <- 10
  nw <- network.initialize(n, directed = FALSE)
  nw <- set.vertex.attribute(nw, "loc", rep(0:1, each = n / 2))

  fit <- ergm(
    nw ~ edges + offset(nodemix("loc", levels2 = -c(1, 3))),
    target.stats = 3,
    offset.coef = -Inf
  )
  coef.form <- coef(fit)
  coef.form[1] <- coef.form[1] - log(40 - 1)

  fit.sim <- simulate(fit, dynamic = FALSE)
  expect_error(
    simulate(
      fit.sim ~ Form(~ edges + offset(nodemix("loc", levels2 = -c(1, 3)))) +
        Persist(~offset(edges)),
      coef = c(coef.form, log(40 - 1)),
      time.slices = 10,
      monitor = "formation",
      nsim = 1,
      dynamic = TRUE
    ),
    NA
  )
})
