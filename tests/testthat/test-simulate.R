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
