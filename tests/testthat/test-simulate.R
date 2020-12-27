#  File tests/testthat/test-simulate.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################

test_that("simulate_formula.network behaves reasonably", {
  nw <- network.initialize(100, directed = FALSE)

  ## if we don't specify 'dynamic=TRUE' we should go to ergm simulation by default
  expect_warning(rv <- simulate(nw ~ edges, coef = c(-3)), paste0("Attempting ", sQuote("ergm"), " simulation instead..."))
  expect_identical(class(rv), "network")
  ## should be no time or lasttoggle since the model isn't durational and we did ergm simulation
  expect_false("time" %in% list.network.attributes(rv))
  expect_false("lasttoggle" %in% list.network.attributes(rv))
  
  ## if we don't specify 'dynamic=TRUE' we should go to ergm simulation by default
  expect_warning(rv <- simulate(nw ~ Form(~edges) + Diss(~edges), coef = c(-3, 1)), paste0("Attempting ", sQuote("ergm"), " simulation instead..."))
  expect_identical(class(rv), "network")
  ## should have time and lasttoggle since the model is durational, but they are basically junk values in this case
  expect_true("time" %in% list.network.attributes(rv))
  expect_true("lasttoggle" %in% list.network.attributes(rv))

  ## if we do specify 'dynamic=TRUE' we should go to tergm simulation
  expect_warning(rv <- simulate(nw ~ edges, coef = c(-3), output = "final", dynamic = TRUE), NA)
  expect_identical(class(rv), "network")
  ## should have time and lasttoggle since we did tergm simulation (even though the model is not durational)
  expect_identical(rv %n% "time", 1L)
  lt <- rv %n% "lasttoggle"
  expect_identical(cbind(as.edgelist(rv), 1L), lt[order(lt[,1], lt[,2]),,drop=FALSE])

  ## if we do specify 'dynamic=TRUE' we should go to tergm simulation
  expect_warning(rv <- simulate(nw ~ Form(~edges) + Diss(~edges), coef = c(-3, 1), output = "final", dynamic = TRUE), NA)
  expect_identical(class(rv), "network")
  ## should have time and lasttoggle since the model is durational (and also because we did tergm simulation)
  expect_identical(rv %n% "time", 1L)
  lt <- rv %n% "lasttoggle"
  expect_identical(cbind(as.edgelist(rv), 1L), lt[order(lt[,1], lt[,2]),,drop=FALSE])

  ## now do another time step with the previous return value
  rv2 <- simulate(rv ~ Form(~edges) + Diss(~edges), coef = c(-3, 1), output = "final", dynamic = TRUE)
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
  rv_stats <- simulate(rv ~ Form(~edges) + Diss(~edges), coef = c(-3, 1), time.slices = 5, output = "stats", stats = TRUE, dynamic = TRUE)
  expect_true(NROW(rv_stats) == 5 && NCOL(rv_stats) == 2)
  expect_identical(colnames(rv_stats), c("Form(edges)", "Diss(edges)"))
  expect_true(all(rv_stats[,1] > rv_stats[,2]))

  rv_stats_mon <- simulate(rv ~ Form(~edges) + Diss(~edges), coef = c(-3, 1), output = "stats", monitor = ~edges, dynamic = TRUE)
  expect_true(NROW(rv_stats_mon) == 1 && NCOL(rv_stats_mon) == 1)
  expect_identical(colnames(rv_stats_mon), c("edges"))
  expect_true(rv_stats_mon[1,1] > 0)
  
  rv_changes <- simulate(rv ~ Form(~edges) + Diss(~edges), coef = c(-3, 0), output = "changes", dynamic = TRUE)
  expect_true(all(rv_changes[,1] == 2L))
  expect_true(all(c(0L, 1L) %in% rv_changes[,4]))

  rv_ergm_state <- simulate(rv ~ Form(~edges) + Diss(~edges), coef = c(-3, 0), time.slices = 3, output = "ergm_state", dynamic = TRUE)
  expect_identical(rv_ergm_state$nw0 %n% "time", 4L)
  expect_true(NROW(rv_ergm_state$el) > 0)
  expect_true(NROW(rv_ergm_state$nw0 %n% "lasttoggle") > NROW(rv_ergm_state$el))
  expect_true(all((rv_ergm_state$nw0 %n% "lasttoggle")[,3] %in% c(1L, 2L, 3L, 4L)))
  expect_true(all(c(1L, 2L, 3L, 4L) %in% (rv_ergm_state$nw0 %n% "lasttoggle")[,3]))

  ## now obtain networkDynamic return value, and restart sim from it, checking various components for basic sensibility
  rv_nwd <- simulate(nw ~ Form(~edges) + Diss(~edges), coef = c(-3, 0), time.slices = 3, dynamic = TRUE)
  expect_identical(get.change.times(rv_nwd), c(1, 2, 3))
  ## technically running simulate_formula.networkDynamic in this next call...
  rv_nw <- simulate(rv_nwd ~ Form(~edges) + Diss(~edges), coef = c(-3, 0), time.slices = 2, output = "final", dynamic = TRUE)
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
  nw2 <- simulate(nw1 ~ Form(~edges) + Diss(~edges), coef = c(-3, 1), output = "final", dynamic = TRUE)
  nw4 <- simulate(nw2 ~ Form(~edges) + Diss(~edges), coef = c(-3, 1), output = "final", time.slices = 2, dynamic = TRUE)
  nw5 <- simulate(nw2 ~ Form(~edges) + Diss(~edges), coef = c(-3, 1), output = "final", time.start = 3, time.slices = 2, dynamic = TRUE)
  
  set.seed(0)
  nw <- network.initialize(100, directed = FALSE)
  nwd1 <- simulate(nw ~ edges, coef = c(-3), dynamic = TRUE)
  nwd2 <- simulate(nwd1 ~ Form(~edges) + Diss(~edges), coef = c(-3, 1), dynamic = TRUE)
  nwd4 <- simulate(nwd2 ~ Form(~edges) + Diss(~edges), coef = c(-3, 1), time.slices = 2, dynamic = TRUE)
  nwd5 <- simulate(nwd2 ~ Form(~edges) + Diss(~edges), coef = c(-3, 1), time.start = 3, time.slices = 2, dynamic = TRUE)
  
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
  
  nw_from_nwd1 <- network.extract.with.lasttoggle(nwd1, at=1, TRUE)
  nw_from_nwd2 <- network.extract.with.lasttoggle(nwd2, at=2, TRUE)
  nw_from_nwd4 <- network.extract.with.lasttoggle(nwd4, at=4, TRUE)
  nw_from_nwd5 <- network.extract.with.lasttoggle(nwd5, at=5, TRUE)

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
