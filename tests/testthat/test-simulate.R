#  File tests/testthat/test-simulate.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2023 Statnet Commons
################################################################################

test_that("simulate_formula.network behaves reasonably", {
  nw <- network.initialize(100, directed = FALSE)

  ## if we don't specify 'dynamic=TRUE' we should go to ergm simulation by default
  expect_warning(rv <- simulate(nw ~ edges, coef = c(-3)), NA)
  expect_identical(class(rv), "network")
  ## should be no time or lasttoggle since the model isn't durational and we did ergm simulation
  expect_false("time" %in% list.network.attributes(rv))
  expect_false("lasttoggle" %in% list.network.attributes(rv))
  
  ## if we don't specify 'dynamic=TRUE' we should trigger an error: allow network series
  expect_error(rv <- simulate(nw ~ Form(~edges) + Persist(~edges), coef = c(-3, 1)), ".*This term requires.*(network series or network dynamic or last toggle information).*dynamic=TRUE.*")

  ## if we don't specify 'dynamic=TRUE' we should trigger an error: disallow network series
  expect_error(rv <- simulate(nw ~ edges + edges.ageinterval(3,4), coef = c(-3, 1)), ".*This term requires.*(network dynamic or last toggle information).*dynamic=TRUE.*")

  ## if we do specify 'dynamic=TRUE' we should go to tergm simulation
  expect_warning(rv <- simulate(nw ~ edges, coef = c(-3), output = "final", dynamic = TRUE), NA)
  expect_identical(class(rv), "network")
  ## should have time and lasttoggle since we did tergm simulation (even though the model is not durational)
  expect_identical(rv %n% "time", 1L)
  lt <- rv %n% "lasttoggle"
  expect_identical(cbind(as.edgelist(rv), 1L), lt[order(lt[,1], lt[,2]),,drop=FALSE])

  ## if we do specify 'dynamic=TRUE' we should go to tergm simulation
  expect_warning(rv <- simulate(nw ~ Form(~edges) + Persist(~edges), coef = c(-3, 1), output = "final", dynamic = TRUE), NA)
  expect_identical(class(rv), "network")
  ## should have time and lasttoggle since the model is durational (and also because we did tergm simulation)
  expect_identical(rv %n% "time", 1L)
  lt <- rv %n% "lasttoggle"
  expect_identical(cbind(as.edgelist(rv), 1L), lt[order(lt[,1], lt[,2]),,drop=FALSE])

  ## now do another time step with the previous return value
  rv2 <- simulate(rv ~ Form(~edges) + Persist(~edges), coef = c(-3, 1), output = "final", dynamic = TRUE)
  expect_identical(class(rv2), "network")
  expect_identical(rv2 %n% "time", 2L)
  rv1 <- rv
  mlt2 <- cbind(as.edgelist((rv2 - rv1) + (rv1 - rv2)), 2L)
  mlt1 <- cbind(as.edgelist(rv1 & rv2), 1L)
  mlt <- rbind(mlt1, mlt2)
  mlt <- mlt[order(mlt[,1], mlt[,2]),,drop=FALSE]
  lt <- rv2 %n% "lasttoggle"
  lt <- lt[order(lt[,1], lt[,2]),,drop=FALSE]
  expect_identical(mlt, lt)

  ## obtain also stats, changes, and ergm_state starting from rv
  rv_stats <- simulate(rv ~ Form(~edges) + Persist(~edges), coef = c(-3, 1), time.slices = 5, output = "stats", stats = TRUE, dynamic = TRUE)
  expect_true(NROW(rv_stats) == 5 && NCOL(rv_stats) == 2)
  expect_identical(colnames(rv_stats), c("Form~edges", "Persist~edges"))
  expect_true(all(rv_stats[,1] > rv_stats[,2]))

  rv_stats_mon <- simulate(rv ~ Form(~edges) + Persist(~edges), coef = c(-3, 1), output = "stats", monitor = ~edges, dynamic = TRUE)
  expect_true(NROW(rv_stats_mon) == 1 && NCOL(rv_stats_mon) == 1)
  expect_identical(colnames(rv_stats_mon), c("edges"))
  expect_true(rv_stats_mon[1,1] > 0)
  
  rv_changes <- simulate(rv ~ Form(~edges) + Persist(~edges), coef = c(-3, 0), output = "changes", dynamic = TRUE)
  expect_true(all(rv_changes[,1] == 2L))
  expect_true(all(c(0L, 1L) %in% rv_changes[,4]))

  rv_ergm_state <- simulate(rv ~ Form(~edges) + Persist(~edges), coef = c(-3, 0), time.slices = 3, output = "ergm_state", dynamic = TRUE)
  expect_identical(rv_ergm_state$nw0 %n% "time", 4L)
  expect_true(NROW(rv_ergm_state$el) > 0)
  expect_true(NROW(rv_ergm_state$nw0 %n% "lasttoggle") > NROW(rv_ergm_state$el))
  expect_true(all((rv_ergm_state$nw0 %n% "lasttoggle")[,3] %in% c(1L, 2L, 3L, 4L)))
  expect_true(all(c(1L, 2L, 3L, 4L) %in% (rv_ergm_state$nw0 %n% "lasttoggle")[,3]))

  ## now obtain networkDynamic return value, and restart sim from it, checking various components for basic sensibility
  rv_nwd <- simulate(nw ~ Form(~edges) + Persist(~edges), coef = c(-3, 0), time.slices = 3, dynamic = TRUE)
  expect_identical(get.change.times(rv_nwd), c(1, 2, 3))
  ## technically running simulate_formula.networkDynamic in this next call...
  rv_nw <- simulate(rv_nwd ~ Form(~edges) + Persist(~edges), coef = c(-3, 0), time.slices = 2, output = "final", dynamic = TRUE)
  expect_identical(rv_nw %n% "time", 5L)
  
  ## every edge in nw3 should have a lasttoggle time in nw5
  nw3 <- network.collapse(rv_nwd, at=3)
  nw5 <- rv_nw
  lt <- nw5 %n% "lasttoggle"
  lt_el <- lt[,1:2,drop=FALSE]
  nw_lt <- network(lt_el, directed = FALSE)
  expect_identical(network.edgecount(nw3 - nw_lt), 0L)
  
  ## any edges in nw3 and not in nw5 or in nw5 and not in nw3 should have lasttoggle times >= 4
  el45 <- as.edgelist((nw3 - nw5) + (nw5 - nw3))
  for(i in seq_len(NROW(el45))) {
    lt_time <- lt[lt[,1] == el45[i,1] & lt[,2] == el45[i,2],3]
    expect_true(identical(lt_time, 4L) || identical(lt_time, 5L))
  }
})

test_that("simulate_formula.networkDynamic behaves reasonably", {
  ## do some matching sims with sim_f.nw and sim_f.nwD and then compare times, edgelists, and lasttoggles
  set.seed(0)
  nw <- network.initialize(100, directed = FALSE)
  nw1 <- simulate(nw ~ edges, coef = c(-3), output = "final", dynamic = TRUE)
  nw2 <- simulate(nw1 ~ Form(~edges) + Persist(~edges), coef = c(-3, 1), output = "final", dynamic = TRUE)
  nw4 <- simulate(nw2 ~ Form(~edges) + Persist(~edges), coef = c(-3, 1), output = "final", time.slices = 2, dynamic = TRUE)
  expect_warning(nw5 <- simulate(nw2 ~ Form(~edges) + Persist(~edges), coef = c(-3, 1), output = "final", time.start = 3, time.slices = 2, dynamic = TRUE))
  
  set.seed(0)
  nw <- network.initialize(100, directed = FALSE)
  nwd1 <- simulate(nw ~ edges, coef = c(-3), dynamic = TRUE)
  nwd2 <- simulate(nwd1 ~ Form(~edges) + Persist(~edges), coef = c(-3, 1), dynamic = TRUE)
  nwd4 <- simulate(nwd2 ~ Form(~edges) + Persist(~edges), coef = c(-3, 1), time.slices = 2, dynamic = TRUE)
  expect_warning(nwd5 <- simulate(nwd2 ~ Form(~edges) + Persist(~edges), coef = c(-3, 1), time.start = 3, time.slices = 2, dynamic = TRUE))
  
  expect_identical(nw1 %n% "time", as.integer(.get.last.obs.time(nwd1)))
  expect_identical(nw2 %n% "time", as.integer(.get.last.obs.time(nwd2)))
  expect_identical(nw4 %n% "time", as.integer(.get.last.obs.time(nwd4)))
  expect_identical(nw5 %n% "time", as.integer(.get.last.obs.time(nwd5)))

  nwlt1 <- nw1 %n% "lasttoggle"
  nwlt2 <- nw2 %n% "lasttoggle"
  nwlt4 <- nw4 %n% "lasttoggle"
  nwlt5 <- nw5 %n% "lasttoggle"
  
  nwel1 <- as.edgelist(nw1)
  nwel2 <- as.edgelist(nw2)
  nwel4 <- as.edgelist(nw4)
  nwel5 <- as.edgelist(nw5)
  
  nw_from_nwd1 <- network.extract.with.lasttoggle(nwd1, at=1)
  nw_from_nwd2 <- network.extract.with.lasttoggle(nwd2, at=2)
  nw_from_nwd4 <- network.extract.with.lasttoggle(nwd4, at=4)
  nw_from_nwd5 <- network.extract.with.lasttoggle(nwd5, at=5)

  nwdlt1 <- nw_from_nwd1 %n% "lasttoggle"
  nwdlt2 <- nw_from_nwd2 %n% "lasttoggle"
  nwdlt4 <- nw_from_nwd4 %n% "lasttoggle"
  nwdlt5 <- nw_from_nwd5 %n% "lasttoggle"

  nwdel1 <- as.edgelist(nw_from_nwd1)
  nwdel2 <- as.edgelist(nw_from_nwd2)
  nwdel4 <- as.edgelist(nw_from_nwd4)
  nwdel5 <- as.edgelist(nw_from_nwd5)
  
  expect_identical(nwel1, nwdel1)  
  expect_identical(nwel2, nwdel2)  
  expect_identical(nwel4, nwdel4)  
  expect_identical(nwel5, nwdel5)  
  
  expect_identical(nwlt1[order(nwlt1[,1], nwlt1[,2]),,drop=FALSE], nwdlt1[order(nwdlt1[,1], nwdlt1[,2]),,drop=FALSE])  
  expect_identical(nwlt2[order(nwlt2[,1], nwlt2[,2]),,drop=FALSE], nwdlt2[order(nwdlt2[,1], nwdlt2[,2]),,drop=FALSE])  

  ## these next two should only include non-edges if they were toggled off on the most 
  ## recent time step, given the way network.extract.with.lasttoggle is implemented
  wtk4 <- logical(NROW(nwlt4))
  for(i in seq_along(wtk4)) {
    wtk4[i] <- any(nwel4[,1] == nwlt4[i,1] & nwel4[,2] == nwlt4[i,2]) || nwlt4[i,3] == 4
  }
  nwlt4 <- nwlt4[wtk4,,drop=FALSE]
  expect_identical(nwlt4[order(nwlt4[,1], nwlt4[,2]),,drop=FALSE], nwdlt4[order(nwdlt4[,1], nwdlt4[,2]),,drop=FALSE])  

  wtk5 <- logical(NROW(nwlt5))
  for(i in seq_along(wtk5)) {
    wtk5[i] <- any(nwel5[,1] == nwlt5[i,1] & nwel5[,2] == nwlt5[i,2]) || nwlt5[i,3] == 5
  }
  nwlt5 <- nwlt5[wtk5,,drop=FALSE]
  expect_identical(nwlt5[order(nwlt5[,1], nwlt5[,2]),,drop=FALSE], nwdlt5[order(nwdlt5[,1], nwdlt5[,2]),,drop=FALSE])  
})

test_that("simulate.network behaves reasonably", {
  nw <- network.initialize(100, directed = FALSE)

  ## obtain some results with the old interface
  set.seed(0)
  old_nw1 <- simulate(nw, formation = ~edges + concurrent, dissolution = ~edges, coef.form = c(-3, 0.5), coef.diss = 1, output = "final")
  old_nw2 <- simulate(old_nw1, formation = ~edges + concurrent, dissolution = ~edges, coef.form = c(-3, 0.5), coef.diss = 1, output = "final")
  old_nw4 <- simulate(old_nw2, formation = ~edges + concurrent, dissolution = ~edges, coef.form = c(-3, 0.5), coef.diss = 1, output = "final", time.slices = 2)
  expect_warning(old_nw5 <- simulate(old_nw2, formation = ~edges + concurrent, dissolution = ~edges, coef.form = c(-3, 0.5), coef.diss = 1, output = "final", time.start = 3, time.slices = 2))
  old_nw6 <- simulate(old_nw5, formation = ~edges + concurrent, dissolution = ~edges, coef.form = attr(old_nw5, "coef.form"), coef.diss = attr(old_nw5, "coef.diss"), output = "final")

  old_stats_nw <- simulate(old_nw6, formation = ~edges + concurrent, dissolution = ~edges, coef.form = c(-3, 0.5), coef.diss = 1, time.slices = 10, output = "final", stats.form = TRUE, stats.diss = TRUE, monitor = ~triangle + degree(1))
  old_stats <- simulate(old_nw6, formation = ~edges + concurrent, dissolution = ~edges, coef.form = c(-3, 0.5), coef.diss = 1, time.slices = 10, output = "stats", stats.form = TRUE, stats.diss = TRUE, monitor = ~triangle + degree(1))

  ## obtain what should be the same results with the new interface
  set.seed(0)
  new_nw1 <- simulate(nw ~ edges + concurrent, coef = c(-3, 0.5), output = "final", dynamic = TRUE)
  new_nw2 <- simulate(new_nw1 ~ Form(~edges + concurrent) + Persist(~edges), coef = c(-3, 0.5, 1), output = "final", dynamic = TRUE)
  new_nw4 <- simulate(new_nw2 ~ Form(~edges + concurrent) + Persist(~edges), coef = c(-3, 0.5, 1), output = "final", time.slices = 2, dynamic = TRUE)
  expect_warning(new_nw5 <- simulate(new_nw2 ~ Form(~edges + concurrent) + Persist(~edges), coef = c(-3, 0.5, 1), output = "final", time.start = 3, time.slices = 2, dynamic = TRUE))
  new_nw6 <- simulate(new_nw2 ~ Form(~edges + concurrent) + Persist(~edges), coef = attr(new_nw5, "coef"), output = "final", basis = new_nw5, dynamic = TRUE)

  new_stats_nw <- simulate(new_nw6 ~ Form(~edges + concurrent) + Persist(~edges), coef = c(-3, 0.5, 1), output = "final", time.slices = 10, stats = TRUE, monitor = ~triangle + degree(1), dynamic = TRUE)
  new_stats <- simulate(new_nw6 ~ Form(~edges + concurrent) + Persist(~edges), coef = c(-3, 0.5, 1), output = "stats", time.slices = 10, stats = TRUE, monitor = ~triangle + degree(1), dynamic = TRUE)
  
  ## test for equality (attributes will differ)
  expect_equal(old_nw1, new_nw1, ignore_attr = TRUE)
  expect_equal(old_nw2, new_nw2, ignore_attr = TRUE)
  expect_equal(old_nw4, new_nw4, ignore_attr = TRUE)
  expect_equal(old_nw5, new_nw5, ignore_attr = TRUE)
  expect_equal(old_nw6, new_nw6, ignore_attr = TRUE)

  expect_identical(attr(old_stats_nw, "stats"), attr(new_stats_nw, "stats"))
  expect_identical(attr(old_stats_nw, "stats.gen"), attr(new_stats_nw, "stats.gen"))
  expect_identical(cbind(attr(old_stats_nw, "stats.form"), attr(old_stats_nw, "stats.diss")), as.matrix(attr(new_stats_nw, "stats.gen")))

  expect_identical(old_stats$stats, new_stats$stats)
  expect_identical(old_stats$stats.gen, new_stats$stats.gen)
  expect_identical(cbind(old_stats$stats.form, old_stats$stats.diss), as.matrix(new_stats$stats.gen))
})

test_that("simulate.networkDynamic behaves reasonably", {
  nw <- network.initialize(100, directed = FALSE)

  ## obtain some results with the old interface
  set.seed(2)
  old_nwD1 <- simulate(nw, formation = ~edges + concurrent, dissolution = ~edges, coef.form = c(-3, 0.5), coef.diss = 1)
  old_nwD2 <- simulate(old_nwD1, formation = ~edges + concurrent, dissolution = ~edges, coef.form = c(-3, 0.5), coef.diss = 1)
  old_nwD4 <- simulate(old_nwD2, formation = ~edges + concurrent, dissolution = ~edges, coef.form = c(-3, 0.5), coef.diss = 1, time.slices = 2)
  expect_warning(old_nwD5 <- simulate(old_nwD2, formation = ~edges + concurrent, dissolution = ~edges, coef.form = c(-3, 0.5), coef.diss = 1, time.start = 3, time.slices = 2))
  old_nwD6 <- simulate(old_nwD5, formation = ~edges + concurrent, dissolution = ~edges, coef.form = attr(old_nwD5, "coef.form"), coef.diss = attr(old_nwD5, "coef.diss"))

  old_nwD7 <- simulate(old_nwD6, formation = ~edges + concurrent, dissolution = ~edges)

  old_stats_nwD <- simulate(old_nwD6, formation = ~edges + concurrent, dissolution = ~edges, coef.form = c(-3, 0.5), coef.diss = 1, time.slices = 10, stats.form = TRUE, stats.diss = TRUE, monitor = ~triangle + degree(1))
  old_stats <- simulate(old_nwD6, formation = ~edges + concurrent, dissolution = ~edges, coef.form = c(-3, 0.5), coef.diss = 1, time.slices = 10, output = "stats", stats.form = TRUE, stats.diss = TRUE, monitor = ~triangle + degree(1))

  old_nwD8 <- simulate(old_stats_nwD, formation = ~edges + concurrent, dissolution = ~edges, monitor = ~triangle + degree(1))

  old_nwD_constr <- simulate(nw, constraints = ~bd(maxout=1), formation = ~edges + concurrent, dissolution = ~edges, coef.form = c(-3, 0.5), coef.diss = 1, time.slices = 10)
  old_nwD_constr2 <- simulate(old_nwD_constr, constraints = ~bd(maxout=1), formation = ~edges + concurrent, dissolution = ~edges, time.slices=10)

  ## obtain what should be the same results with the new interface
  set.seed(2)
  new_nwD1 <- simulate(nw ~ edges + concurrent, coef = c(-3, 0.5), dynamic = TRUE)
  new_nwD2 <- simulate(new_nwD1 ~ Form(~edges + concurrent) + Persist(~edges), coef = c(-3, 0.5, 1), dynamic = TRUE)
  new_nwD4 <- simulate(new_nwD2 ~ Form(~edges + concurrent) + Persist(~edges), coef = c(-3, 0.5, 1), time.slices = 2, dynamic = TRUE)
  expect_warning(new_nwD5 <- simulate(new_nwD2 ~ Form(~edges + concurrent) + Persist(~edges), coef = c(-3, 0.5, 1), time.start = 3, time.slices = 2, dynamic = TRUE))
  new_nwD6 <- simulate(new_nwD2 ~ Form(~edges + concurrent) + Persist(~edges), coef = attr(new_nwD5, "coef"), basis = new_nwD5, dynamic = TRUE)

  new_nwD7 <- simulate(new_nwD6 ~ Form(~edges + concurrent) + Persist(~edges), dynamic = TRUE)

  new_stats_nwD <- simulate(new_nwD6 ~ Form(~edges + concurrent) + Persist(~edges), coef = c(-3, 0.5, 1), time.slices = 10, stats = TRUE, monitor = ~triangle + degree(1), dynamic = TRUE)
  new_stats <- simulate(new_nwD6 ~ Form(~edges + concurrent) + Persist(~edges), coef = c(-3, 0.5, 1), output = "stats", time.slices = 10, stats = TRUE, monitor = ~triangle + degree(1), dynamic = TRUE)

  new_nwD8 <- simulate(new_stats_nwD ~ Form(~edges + concurrent) + Persist(~edges), monitor = ~triangle + degree(1), dynamic = TRUE)
  
  new_nwD_constr <- simulate(nw ~ Form(~edges + concurrent) + Persist(~edges), coef = c(-3, 0.5, 1), constraints = ~bd(maxout=1), time.slices = 10, dynamic = TRUE)
  new_nwD_constr2 <- simulate(new_nwD_constr ~ Form(~edges + concurrent) + Persist(~edges), constraints = ~bd(maxout=1), time.slices = 10, dynamic = TRUE)

  ## test for equality (attributes will differ)
  expect_equal(old_nwD1, new_nwD1, ignore_attr = TRUE)
  expect_equal(old_nwD2, new_nwD2, ignore_attr = TRUE)
  expect_equal(old_nwD4, new_nwD4, ignore_attr = TRUE)
  expect_equal(old_nwD5, new_nwD5, ignore_attr = TRUE)
  expect_equal(old_nwD6, new_nwD6, ignore_attr = TRUE)
  expect_equal(old_nwD7, new_nwD7, ignore_attr = TRUE)
  expect_equal(old_nwD8, new_nwD8, ignore_attr = TRUE)

  expect_true(!is.null(attr(old_nwD8, "stats")))
  expect_true(!is.null(attr(new_nwD8, "stats")))

  expect_equal(attr(old_nwD8, "stats"), attr(new_nwD8, "stats"))

  expect_equal(attr(old_stats_nwD, "stats"), attr(new_stats_nwD, "stats"))
  expect_equal(attr(old_stats_nwD, "stats.gen"), attr(new_stats_nwD, "stats.gen"))
  expect_equal(cbind(attr(old_stats_nwD, "stats.form"), attr(old_stats_nwD, "stats.diss")), as.matrix(attr(new_stats_nwD, "stats.gen")))

  expect_equal(old_stats$stats, new_stats$stats)
  expect_equal(old_stats$stats.gen, new_stats$stats.gen)
  expect_equal(cbind(old_stats$stats.form, old_stats$stats.diss), as.matrix(new_stats$stats.gen))

  expect_equal(old_nwD_constr, new_nwD_constr, ignore_attr = TRUE)
  expect_equal(old_nwD_constr2, new_nwD_constr2, ignore_attr = TRUE)

  expect_true(summary(network.collapse(old_nwD_constr, at = 10) ~ edges) > 0)
  expect_true(summary(network.collapse(old_nwD_constr2, at = 20) ~ edges) > 0)
  expect_true(summary(network.collapse(new_nwD_constr, at = 10) ~ edges) > 0)
  expect_true(summary(network.collapse(new_nwD_constr2, at = 20) ~ edges) > 0)

  expect_equal(unname(summary(network.collapse(old_nwD_constr, at = 10) ~ concurrent)), 0)
  expect_equal(unname(summary(network.collapse(old_nwD_constr2, at = 20) ~ concurrent)), 0)
  expect_equal(unname(summary(network.collapse(new_nwD_constr, at = 10) ~ concurrent)), 0)
  expect_equal(unname(summary(network.collapse(new_nwD_constr2, at = 20) ~ concurrent)), 0)
})

test_that("Dynamic simulation error messages are correct", {
  nw <- network.initialize(100, directed = FALSE)

  expect_error(simulate(nw ~ Form(~edges) + Persist(~edges), coef = c(5, 5), dynamic = TRUE, control = control.simulate.formula.tergm(MCMC.maxedges=2)),
               "^Number of edges in a simulated network exceeds the maximum set by the 'MCMC.maxedges' control parameter.$")

  expect_error(simulate(nw ~ Form(~edges) + Persist(~edges), coef = c(5, 5), dynamic = TRUE, control = control.simulate.formula.tergm(MCMC.maxchanges=2)),
               "^Logging of changes in the network has been requested, and the storage capacity specified by 'MCMC.maxchanges' has been exceeded.$")
})

statnet.common::opttest({
test_that("simulate.tergm behaves reasonably", {
  nw <- network.initialize(100, directed = FALSE)

  nw1 <- simulate(nw ~ Form(~edges) + Persist(~edges), coef = c(-5, 1), dynamic = TRUE, output = "final")
  nw2 <- simulate(nw1 ~ Form(~edges) + Persist(~edges), coef = c(-5, 1), dynamic = TRUE, output = "final")
  nw3 <- simulate(nw2 ~ Form(~edges) + Persist(~edges), coef = c(-5, 1), dynamic = TRUE, output = "final")
  nw4 <- simulate(nw3 ~ Form(~edges) + Persist(~edges), coef = c(-5, 1), dynamic = TRUE, output = "final")

  nw100 <- simulate(nw ~ Form(~edges) + Persist(~edges), coef = c(-5, 1), dynamic = TRUE, output = "final", time.slices = 100)

  nwx <- network(100, directed = FALSE)
  
  nwL <- list(nw1, nw2, nw3, nw4)

  # pretty fast
  set.seed(0)
  CMLE <- tergm(nwL ~ Form(~edges) + Persist(~edges), estimate="CMLE")

  # this one can take a few minutes...
  set.seed(0)
  EGMME <- tergm(nw100 ~ Form(~edges) + Persist(~edges), targets = ~edges + mean.age, target.stats = c(120.9842, 3.718282), estimate = "EGMME")

  set.seed(0)
  rv_EGMME <- simulate(EGMME, output = "final")

  set.seed(0)
  rv_EGMME_x <- simulate(EGMME, nw.start = nwx, output = "final")
  
  set.seed(0)
  rv_CMLE_first <- simulate(CMLE, nw.start = "first", output = "final")  
  set.seed(0)
  rv_CMLE_last <- simulate(CMLE, nw.start = "last", output = "final")  

  set.seed(0)
  rv_CMLE_x <- simulate(CMLE, nw.start = nwx, output = "final")  

  set.seed(0)
  expect_error(simulate(CMLE, output = "final"))
  set.seed(0)
  expect_error(simulate(CMLE, nw.start = 0, output = "final"))
  set.seed(0)
  rv_CMLE_1 <- simulate(CMLE, nw.start = 1, output = "final")  
  set.seed(0)
  rv_CMLE_2 <- simulate(CMLE, nw.start = 2, output = "final")  
  set.seed(0)
  rv_CMLE_3 <- simulate(CMLE, nw.start = 3, output = "final")  
  set.seed(0)
  rv_CMLE_4 <- simulate(CMLE, nw.start = 4, output = "final")  
  set.seed(0)
  expect_error(simulate(CMLE, nw.start = 5, output = "final"))

  set.seed(0)
  rv_E_nw <- simulate(EGMME$network ~ Form(~edges) + Persist(~edges), coef = coef(EGMME), monitor = EGMME$targets, output = "final", dynamic = TRUE)

  set.seed(0)
  rv_E_x_nw <- simulate(nwx ~ Form(~edges) + Persist(~edges), coef = coef(EGMME), monitor = EGMME$targets, output = "final", dynamic = TRUE)

  set.seed(0)
  rv_C_first_nw <- simulate(nw1 ~ Form(~edges) + Persist(~edges), coef = coef(CMLE), output = "final", dynamic = TRUE)
  set.seed(0)
  rv_C_last_nw <- simulate(nw4 ~ Form(~edges) + Persist(~edges), coef = coef(CMLE), output = "final", dynamic = TRUE)

  set.seed(0)
  rv_C_x_nw <- simulate(nwx ~ Form(~edges) + Persist(~edges), coef = coef(CMLE), output = "final", dynamic = TRUE)

  set.seed(0)
  rv_C_1_nw <- simulate(nw1 ~ Form(~edges) + Persist(~edges), coef = coef(CMLE), output = "final", dynamic = TRUE)
  set.seed(0)
  rv_C_2_nw <- simulate(nw2 ~ Form(~edges) + Persist(~edges), coef = coef(CMLE), output = "final", dynamic = TRUE)
  set.seed(0)
  rv_C_3_nw <- simulate(nw3 ~ Form(~edges) + Persist(~edges), coef = coef(CMLE), output = "final", dynamic = TRUE)
  set.seed(0)
  rv_C_4_nw <- simulate(nw4 ~ Form(~edges) + Persist(~edges), coef = coef(CMLE), output = "final", dynamic = TRUE)

  ## subset network attributes to those present in direct simulations of the original networks to avoid
  ## considering unrelated attributes introduced during the CMLE fit
  names_to_keep <- list.network.attributes(rv_C_first_nw)

  rv_CMLE_first$gal <- rv_CMLE_first$gal[names_to_keep]
  rv_CMLE_last$gal <- rv_CMLE_last$gal[names_to_keep]
  rv_CMLE_x$gal <- rv_CMLE_x$gal[names_to_keep]
  rv_CMLE_1$gal <- rv_CMLE_1$gal[names_to_keep]
  rv_CMLE_2$gal <- rv_CMLE_2$gal[names_to_keep]
  rv_CMLE_3$gal <- rv_CMLE_3$gal[names_to_keep]
  rv_CMLE_4$gal <- rv_CMLE_4$gal[names_to_keep]

  rv_C_first_nw$gal <- rv_C_first_nw$gal[names_to_keep]
  rv_C_last_nw$gal <- rv_C_last_nw$gal[names_to_keep]
  rv_C_x_nw$gal <- rv_C_x_nw$gal[names_to_keep]
  rv_C_1_nw$gal <- rv_C_1_nw$gal[names_to_keep]
  rv_C_2_nw$gal <- rv_C_2_nw$gal[names_to_keep]
  rv_C_3_nw$gal <- rv_C_3_nw$gal[names_to_keep]
  rv_C_4_nw$gal <- rv_C_4_nw$gal[names_to_keep]
  
  expect_equal(rv_EGMME, rv_E_nw)
  expect_equal(rv_EGMME_x, rv_E_x_nw)
  
  expect_equal(rv_CMLE_first, rv_C_first_nw)  
  expect_equal(rv_CMLE_last, rv_C_last_nw)  

  expect_equal(rv_CMLE_x, rv_C_x_nw)  

  expect_equal(rv_CMLE_1, rv_C_1_nw)  
  expect_equal(rv_CMLE_2, rv_C_2_nw)  
  expect_equal(rv_CMLE_3, rv_C_3_nw)  
  expect_equal(rv_CMLE_4, rv_C_4_nw)  
})
}, "simulate.tergm")
