#  File tests/testthat/test-stergm-init.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################

test_that("stergm correctly handles initial coefficient specifications", {
  set.seed(2)
  nw <- network.initialize(100, directed = FALSE)
  nw %v% "attr" <- rep(1:5, length.out = 100)
  nw %v% "cov" <- runif(100)
  nw <- simulate(nw ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "final", time.slices = 50)
  nw1 <- simulate(nw ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "final", time.slices = 1)
  nw2 <- simulate(nw1 ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "final", time.slices = 1)
  nw3 <- simulate(nw2 ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "final", time.slices = 1)
  nwlist <- list(nw1, nw2, nw3)
  
  ## valid specification tests
  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.form = c(1,1),
               offset.coef.diss = c(-1),
               estimate = "CMLE",
               control = control.stergm())

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,-1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.form = c(1,1),
               estimate = "CMLE",
               control = control.stergm(init.diss = c(NA, -1, NA, NA, NA, NA, NA)))

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,-1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.diss = c(-1),
               estimate = "CMLE",
               control = control.stergm(init.form = c(NA, NA, 1, NA, 1)))

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,-1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               estimate = "CMLE",
               control = control.stergm(init.form = c(NA, NA, 1, NA, 1),
                                        init.diss = c(NA, -1, NA, NA, NA, NA, NA)))

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,-1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.form = c(1,1),
               estimate = "CMLE",
               control = control.stergm(CMLE.diss.ergm = control.ergm(init = c(NA, -1, NA, NA, NA, NA, NA))))

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,-1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.diss = c(-1),
               estimate = "CMLE",
               control = control.stergm(CMLE.form.ergm = control.ergm(init = c(NA, NA, 1, NA, 1))))

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,-1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               estimate = "CMLE",
               control = control.stergm(CMLE.form.ergm = control.ergm(init = c(NA, NA, 1, NA, 1)),
                                        CMLE.diss.ergm = control.ergm(init = c(NA, -1, NA, NA, NA, NA, NA))))

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,-1), check.attributes = FALSE)

  ## override tests
  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.form = c(1,1),
               offset.coef.diss = c(2),
               estimate = "CMLE",
               control = control.stergm(init.diss = c(NA, -1, NA, NA, NA, NA, NA)))

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,2), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.form = c(5,-3),
               offset.coef.diss = c(-1),
               estimate = "CMLE",
               control = control.stergm(init.form = c(NA, NA, 1, NA, 1)))

  expect_equal(coef(rv)[c(3,5,7)], c(5,-3,-1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               estimate = "CMLE",
               offset.coef.form = c(5,-3),
               offset.coef.diss = c(1),               
               control = control.stergm(init.form = c(NA, NA, 1, NA, 1),
                                        init.diss = c(NA, -1, NA, NA, NA, NA, NA)))

  expect_equal(coef(rv)[c(3,5,7)], c(5,-3,1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + nodecov("cov") + nodematch("attr", diff = TRUE),
               offset.coef.form = c(4,-7),
               estimate = "CMLE",
               control = control.stergm(init.form = c(NA, NA, 1, NA, 1)))

  expect_equal(coef(rv)[c(3,5)], c(4,-7), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.form = c(1,1),
               offset.coef.diss = c(2),
               estimate = "CMLE",
               control = control.stergm(CMLE.diss.ergm = control.ergm(init = c(NA, -1, NA, NA, NA, NA, NA))))

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,2), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.form = c(5,-3),
               offset.coef.diss = c(-1),
               estimate = "CMLE",
               control = control.stergm(CMLE.form.ergm = control.ergm(init = c(NA, NA, 1, NA, 1))))

  expect_equal(coef(rv)[c(3,5,7)], c(5,-3,-1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               estimate = "CMLE",
               offset.coef.form = c(5,-3),
               offset.coef.diss = c(1),               
               control = control.stergm(CMLE.form.ergm = control.ergm(init = c(NA, NA, 1, NA, 1)),
                                        CMLE.diss.ergm = control.ergm(init = c(NA, -1, NA, NA, NA, NA, NA))))

  expect_equal(coef(rv)[c(3,5,7)], c(5,-3,1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + nodecov("cov") + nodematch("attr", diff = TRUE),
               offset.coef.form = c(4,-7),
               estimate = "CMLE",
               control = control.stergm(CMLE.form.ergm = control.ergm(init = c(NA, NA, 1, NA, 1))))

  expect_equal(coef(rv)[c(3,5)], c(4,-7), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.form = c(1,1),
               estimate = "CMLE",
               control = control.stergm(init.diss = c(NA, -2, NA, NA, NA, NA, NA),
                                        CMLE.diss.ergm = control.ergm(init = c(NA, -1, NA, NA, NA, NA, NA))))

  expect_equal(coef(rv)[c(3,5,7)], c(1,1,-1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.diss = c(-1),
               estimate = "CMLE",
               control = control.stergm(init.form = c(NA, NA, 1, NA, 1),
                                        CMLE.form.ergm = control.ergm(init = c(NA, NA, 5, NA, -3))))

  expect_equal(coef(rv)[c(3,5,7)], c(5,-3,-1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               estimate = "CMLE",
               control = control.stergm(init.form = c(NA, NA, 7, NA, 8),
                                        init.diss = c(NA, 9, NA, NA, NA, NA, NA),
                                        CMLE.form.ergm = control.ergm(init = c(NA, NA, 5, NA, -3)),
                                        CMLE.diss.ergm = control.ergm(init = c(NA, 1, NA, NA, NA, NA, NA))))

  expect_equal(coef(rv)[c(3,5,7)], c(5,-3,1), check.attributes = FALSE)

  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + nodecov("cov") + nodematch("attr", diff = TRUE),
               estimate = "CMLE",
               control = control.stergm(init.form = c(NA, NA, 1, NA, 1),
                                        CMLE.form.ergm = control.ergm(init = c(NA, NA, 4, NA, -7))))

  expect_equal(coef(rv)[c(3,5)], c(4,-7), check.attributes = FALSE)
  
  ## test errors
  expect_error(rv <- stergm(nwlist,
                            formation = ~edges + offset(nodefactor("attr"), c(2,4)),
                            dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
                            offset.coef.form = c(1),
                            offset.coef.diss = c(-1),
                            estimate = "CMLE",
                            control = control.stergm()), 
                            "Invalid offset parameter vector offset.coef: wrong number of parameters: expected 3 got 2.")
               
  expect_error(rv <- stergm(nwlist,
                            formation = ~edges + offset(nodefactor("attr"), c(2,4)),
                            dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
                            offset.coef.form = c(1,1),
                            estimate = "CMLE",
                            control = control.stergm()), 
                            "Invalid offset parameter vector offset.coef: wrong number of parameters: expected 3 got 2.")
               
  expect_error(rv <- stergm(nwlist,
                            formation = ~edges + offset(nodefactor("attr"), c(2,4)),
                            dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
                            offset.coef.form = c(1,1),
                            estimate = "CMLE",
                            control = control.stergm(init.diss = c(NA, NA, NA, NA, NA, NA, NA))), 
                            "Dissolution model contains offsets whose coefficients have not been specified.")
               
  expect_error(rv <- stergm(nwlist,
                            formation = ~edges + offset(nodefactor("attr"), c(2,4)),
                            dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
                            offset.coef.diss = c(-1),
                            estimate = "CMLE",
                            control = control.stergm(init.form = c(NA, NA, NA, NA, 1))), 
                            "Formation model contains offsets whose coefficients have not been specified.")
               
  expect_error(rv <- stergm(nwlist,
                            formation = ~edges + offset(nodefactor("attr"), c(2,4)),
                            dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
                            estimate = "CMLE",
                            control = control.stergm(init.form = c(NA, NA, 1, NA, 1),
                                                     init.diss = c(NA, -1, NA, NA, NA, NA))), 
                            "Incorrect length of init.diss passed to stergm\\(\\); expected 7, got 6.")
    
  ## test with tergm() redefined to simply return what it's passed
  actualtergm <- tergm::tergm
  unlockBinding("tergm", environment(stergm))
  environment(stergm)$tergm <- function(...) list(...)
  lockBinding("tergm", environment(stergm))
  
  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               estimate = "CMLE",
               control = control.stergm(CMLE.form.ergm = control.ergm(init = c(NA, 2, 1, NA, 1)),
                                        CMLE.diss.ergm = control.ergm(init = c(NA, -1, NA, NA, 3, NA, NA))))

  expect_equal(rv$offset.coef, c(1,1,-1), check.attributes = FALSE)
  expect_equal(rv$control$init, c(NA,2,1,NA,1,NA,-1,NA,NA,3,NA,NA), check.attributes = FALSE)
  
  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               offset.coef.form = c(5,-3),
               offset.coef.diss = c(-1),
               estimate = "CMLE",
               control = control.stergm(CMLE.form.ergm = control.ergm(init = c(NA, NA, 1, 2, 1))))

  expect_equal(rv$offset.coef, c(5,-3,-1), check.attributes = FALSE)
  expect_equal(rv$control$init, c(NA,NA,5,2,-3,NA,-1,NA,NA,NA,NA,NA), check.attributes = FALSE)
  
  rv <- stergm(nwlist,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               estimate = "CMLE",
               control = control.stergm(init.form = c(NA, 1, 7, 2, 8),
                                        init.diss = c(NA, 9, NA, 5, NA, 6, NA),
                                        CMLE.form.ergm = control.ergm(init = c(1, NA, 5, 3, -3)),
                                        CMLE.diss.ergm = control.ergm(init = c(NA, 1, NA, 99, NA, NA, 8))))

  expect_equal(rv$offset.coef, c(5,-3,1), check.attributes = FALSE)
  expect_equal(rv$control$init, c(1, NA, 5, 3, -3, NA, 1, NA, 99, NA, NA, 8), check.attributes = FALSE)
  
  ## test EGMME
  
  ## valid specification tests
  rv <- stergm(nw1,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               targets = ~nodefactor("attr") + nodematch("attr", diff = TRUE),
               offset.coef.form = c(1,1),
               offset.coef.diss = c(-1),
               estimate = "EGMME",
               control = control.stergm())

  expect_equal(rv$offset.coef, c(1,1,-1), check.attributes = FALSE)
  expect_equal(rv$control$init, NULL, check.attributes = FALSE)

  rv <- stergm(nw1,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               targets = ~nodefactor("attr") + nodematch("attr", diff = TRUE),
               offset.coef.form = c(1,1),
               estimate = "EGMME",
               control = control.stergm(init.diss = c(NA, -1, NA, NA, NA, NA, NA)))

  expect_equal(rv$offset.coef, c(1,1,-1), check.attributes = FALSE)
  expect_equal(rv$control$init, c(NA,NA,1,NA,1,NA,-1,NA,NA,NA,NA,NA), check.attributes = FALSE)

  rv <- stergm(nw1,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               targets = ~nodefactor("attr") + nodematch("attr", diff = TRUE),
               offset.coef.diss = c(-1),
               estimate = "EGMME",
               control = control.stergm(init.form = c(NA, 5, 1, NA, 1)))

  expect_equal(rv$offset.coef, c(1,1,-1), check.attributes = FALSE)
  expect_equal(rv$control$init, c(NA,5,1,NA,1,NA,-1,NA,NA,NA,NA,NA), check.attributes = FALSE)

  rv <- stergm(nw1,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               targets = ~nodefactor("attr") + nodematch("attr", diff = TRUE),
               estimate = "EGMME",
               control = control.stergm(init.form = c(NA, NA, 1, 11, 1),
                                        init.diss = c(NA, -1, NA, 12, NA, NA, NA)))

  expect_equal(rv$offset.coef, c(1,1,-1), check.attributes = FALSE)
  expect_equal(rv$control$init, c(NA,NA,1,11,1,NA,-1,NA,12,NA,NA,NA), check.attributes = FALSE)

  ## override tests
  rv <- stergm(nw1,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               targets = ~nodefactor("attr") + nodematch("attr", diff = TRUE),
               offset.coef.form = c(1,1),
               offset.coef.diss = c(2),
               estimate = "EGMME",
               control = control.stergm(init.diss = c(NA, -1, NA, NA, NA, NA, NA)))

  expect_equal(rv$offset.coef, c(1,1,2), check.attributes = FALSE)
  expect_equal(rv$control$init, c(NA,NA,1,NA,1,NA,2,NA,NA,NA,NA,NA), check.attributes = FALSE)

  rv <- stergm(nw1,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               targets = ~nodefactor("attr") + nodematch("attr", diff = TRUE),
               offset.coef.form = c(5,-3),
               offset.coef.diss = c(-1),
               estimate = "EGMME",
               control = control.stergm(init.form = c(NA, 2, 1, 3, 1)))

  expect_equal(rv$offset.coef, c(5,-3,-1), check.attributes = FALSE)
  expect_equal(rv$control$init, c(NA,2,5,3,-3,NA,-1,NA,NA,NA,NA,NA), check.attributes = FALSE)

  rv <- stergm(nw1,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
               targets = ~nodefactor("attr") + nodematch("attr", diff = TRUE),
               estimate = "EGMME",
               offset.coef.form = c(5,-3),
               offset.coef.diss = c(1),               
               control = control.stergm(init.form = c(NA, NA, 1, NA, 1),
                                        init.diss = c(NA, -1, NA, NA, NA, NA, 8)))

  expect_equal(rv$offset.coef, c(5,-3,1), check.attributes = FALSE)
  expect_equal(rv$control$init, c(NA,NA,5,NA,-3,NA,1,NA,NA,NA,NA,8), check.attributes = FALSE)

  rv <- stergm(nw1,
               formation = ~edges + offset(nodefactor("attr"), c(2,4)),
               dissolution = ~edges + nodecov("cov") + nodematch("attr", diff = TRUE),
               targets = ~edges + nodefactor("attr") + nodematch("attr", diff = TRUE),
               offset.coef.form = c(4,-7),
               estimate = "EGMME",
               control = control.stergm(init.form = c(NA, NA, 1, NA, 1)))

  expect_equal(rv$offset.coef, c(4,-7), check.attributes = FALSE)
  expect_equal(rv$control$init, c(NA,NA,4,NA,-7,NA,NA,NA,NA,NA,NA,NA), check.attributes = FALSE)
  
  ## test errors               
  expect_error(rv <- stergm(nw1,
                            formation = ~edges + offset(nodefactor("attr"), c(2,4)),
                            dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
                            targets = ~nodefactor("attr") + nodematch("attr", diff = TRUE),
                            offset.coef.form = c(1,1),
                            estimate = "EGMME",
                            control = control.stergm(init.diss = c(NA, NA, NA, NA, NA, NA, NA))), 
                            "Dissolution model contains offsets whose coefficients have not been specified.")
               
  expect_error(rv <- stergm(nw1,
                            formation = ~edges + offset(nodefactor("attr"), c(2,4)),
                            dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
                            targets = ~nodefactor("attr") + nodematch("attr", diff = TRUE),
                            offset.coef.diss = c(-1),
                            estimate = "EGMME",
                            control = control.stergm(init.form = c(NA, NA, NA, NA, 1))), 
                            "Formation model contains offsets whose coefficients have not been specified.")
               
  expect_error(rv <- stergm(nw1,
                            formation = ~edges + offset(nodefactor("attr"), c(2,4)),
                            dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
                            targets = ~nodefactor("attr") + nodematch("attr", diff = TRUE),
                            estimate = "EGMME",
                            control = control.stergm(init.form = c(NA, NA, 1, NA, 1),
                                                     init.diss = c(NA, -1, NA, NA, NA, NA))), 
                            "Incorrect length of init.diss passed to stergm\\(\\); expected 7, got 6.")  

  # restore tergm definition
  unlockBinding("tergm", environment(stergm))
  environment(stergm)$tergm <- actualtergm
  lockBinding("tergm", environment(stergm))

  ## test a couple of errors that are caught lower down
  expect_error(rv <- stergm(nw1,
                            formation = ~edges + offset(nodefactor("attr"), c(2,4)),
                            dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
                            targets = ~nodefactor("attr") + nodematch("attr", diff = TRUE),
                            offset.coef.form = c(1),
                            offset.coef.diss = c(-1),
                            estimate = "EGMME",
                            control = control.stergm()), 
                            "Invalid offset parameter vector offset.coef: wrong number of parameters: expected 3 got 2.")
               
  expect_error(rv <- stergm(nw1,
                            formation = ~edges + offset(nodefactor("attr"), c(2,4)),
                            dissolution = ~edges + offset(nodecov("cov")) + nodematch("attr", diff = TRUE),
                            targets = ~nodefactor("attr") + nodematch("attr", diff = TRUE),
                            offset.coef.form = c(1,1),
                            estimate = "EGMME",
                            control = control.stergm()), 
                            "Invalid offset parameter vector offset.coef: wrong number of parameters: expected 3 got 2.")
  
})
