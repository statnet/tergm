#  File tests/testthat/test-extract-formulae.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################
#  File tests/testthat/test-extract-formulae.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

test_that(".extract.fd.formulae behaves reasonably", {
  
  F1 <- ~Form(~edges) + Persist(~edges)
  
  F2 <- ~Form(~edges + gwesp(0,fixed=TRUE)) + Persist(~edges)
  
  F3 <- ~Form(~edges + gwesp(0,fixed=TRUE)) + Persist(~edges) + edges - triangle
  
  X <- ~edges
  
  F4 <- ~Form(X) + Persist(X)
  
  F5 <- ~Form(~offset(edges) + gwesp(0,fixed=TRUE)) + Persist(~offset(edges))

  Y <- ~offset(edges)
  
  F6 <- ~Form(Y) - Persist(~edges)
  
  F7 <- ~Form(Y) - Persist(Y)
  
  F8 <- ~-Form(~edges - triangle) + Persist(~edges) + Persist(~triangle) + Persist(~-offset(triangle))
  
  F9 <- ~-Form(~offset(edges) - offset(triangle) + offset(gwesp(0,fixed=TRUE))) + Persist(~edges)
  
  F10 <- ~-Form(~-edges + triangle) + Form(~-gwesp(0,fixed=TRUE)) - Persist(~-edges) + Persist(~triangle) + Persist(~-gwesp(0,fixed=TRUE))
  
  R1 <- .extract.fd.formulae(F1)
  expect_equal(~edges, R1$form)
  expect_equal(~edges, R1$diss)
  expect_equal(~., R1$nonsep)
  expect_equal(~edges+edges, R1$all)
  
  R2 <- .extract.fd.formulae(F2)
  expect_equal(~edges + gwesp(0,fixed=TRUE), R2$form)
  expect_equal(~edges, R2$diss)
  expect_equal(~., R2$nonsep)
  expect_equal(~edges+gwesp(0,fixed=TRUE)+edges, R2$all)
  
  R3 <- .extract.fd.formulae(F3)
  expect_equal(~edges + gwesp(0,fixed=TRUE), R3$form)
  expect_equal(~edges, R3$diss)
  expect_equal(~edges - triangle, R3$nonsep)
  expect_equal(~edges+gwesp(0,fixed=TRUE)+edges+edges-triangle, R3$all)
  
  R4 <- .extract.fd.formulae(F4)
  expect_equal(~edges, R4$form)
  expect_equal(~edges, R4$diss)
  expect_equal(~., R4$nonsep)
  expect_equal(~edges+edges, R4$all)
  
  R5 <- .extract.fd.formulae(F5)
  expect_equal(~offset(edges) + gwesp(0,fixed=TRUE), R5$form)
  expect_equal(~offset(edges), R5$diss)
  expect_equal(~., R5$nonsep)
  expect_equal(~offset(edges) + gwesp(0,fixed=TRUE) + offset(edges), R5$all)
  
  R6 <- .extract.fd.formulae(F6)
  expect_equal(~offset(edges), R6$form)
  expect_equal(~-edges, R6$diss)
  expect_equal(~., R6$nonsep)
  expect_equal(~offset(edges) - edges, R6$all)
  
  R7 <- .extract.fd.formulae(F7)
  expect_equal(~offset(edges), R7$form)
  expect_equal(~-offset(edges), R7$diss)
  expect_equal(~., R7$nonsep)
  expect_equal(~offset(edges) - offset(edges), R7$all)
  
  R8 <- .extract.fd.formulae(F8)
  expect_equal(~-edges+triangle, R8$form)
  expect_equal(~edges+triangle-offset(triangle), R8$diss)
  expect_equal(~., R8$nonsep)
  expect_equal(~-edges+triangle+edges+triangle-offset(triangle), R8$all)
  
  R9 <- .extract.fd.formulae(F9)
  expect_equal(~-offset(edges)+offset(triangle)-offset(gwesp(0,fixed=TRUE)), R9$form)
  expect_equal(~edges, R9$diss)
  expect_equal(~., R9$nonsep)
  expect_equal(~-offset(edges)+offset(triangle)-offset(gwesp(0,fixed=TRUE))+edges, R9$all)
  
  R10 <- .extract.fd.formulae(F10)
  expect_equal(~edges-triangle-gwesp(0,fixed=TRUE), R10$form)
  expect_equal(~edges+triangle-gwesp(0,fixed=TRUE), R10$diss)
  expect_equal(~., R10$nonsep)
  expect_equal(~edges-triangle-gwesp(0,fixed=TRUE)+edges+triangle-gwesp(0,fixed=TRUE), R10$all)
  
  F11 <- ~Form(~edges) + Persist(~edges) + edges
  environment(F11) <- new.env()
  
  R11 <- .extract.fd.formulae(F11)
  
  expect_identical(environment(F11), environment(R11$form))
  expect_identical(environment(F11), environment(R11$diss))
  expect_identical(environment(F11), environment(R11$nonsep))
  expect_error(expect_identical(environment(F11), environment(R11$all)))
  
  a <- ~edges
  b <- ~triangle
  environment(a) <- new.env()
  environment(b) <- new.env()
  
  F12 <- ~Form(a) + Persist(b) + edges
  environment(F12) <- new.env()
  environment(F12)$a <- a
  environment(F12)$b <- b
  
  R12 <- .extract.fd.formulae(F12)
  
  expect_error(expect_identical(environment(F12), environment(R12$form)))
  expect_error(expect_identical(environment(F12), environment(R12$diss)))
  expect_identical(environment(F12), environment(R12$nonsep))
  expect_error(expect_identical(environment(F12), environment(R12$all)))
  expect_identical(environment(a), environment(R12$form))
  expect_identical(environment(b), environment(R12$diss))
})
