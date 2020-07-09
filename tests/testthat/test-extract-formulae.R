#  File tests/testthat/test-extract-formulae.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

context("test-extract-formulae.R")

test_that(".extract.fd.formulae behaves reasonably", {
  
  F1 <- ~Form(~edges) + Diss(~edges)
  
  F2 <- ~Form(~edges + gwesp(0,fixed=TRUE)) + Diss(~edges)
  
  F3 <- ~Form(~edges + gwesp(0,fixed=TRUE)) + Diss(~edges) + edges
  
  X <- ~edges
  
  F4 <- ~Form(X) + Diss(X)
  
  F5 <- ~Form(~offset(edges) + gwesp(0,fixed=TRUE)) + offset(Diss(~edges))
  
  F6 <- ~offset(Form(X)) - Diss(~edges)
  
  Y <- ~offset(edges)
  
  F7 <- ~offset(Form(Y)) - Diss(Y)
  
  F8 <- ~-Form(~edges - triangle) + Diss(~edges) + Diss(~triangle) + offset(Diss(~-triangle))
  
  F9 <- ~-offset(Form(~edges - triangle + offset(gwesp(0,fixed=TRUE)))) + Diss(~edges)
  
  F10 <- ~-Form(~-edges + triangle) + Form(~-gwesp(0,fixed=TRUE)) - Diss(~-edges) + Diss(~triangle) + Diss(~-gwesp(0,fixed=TRUE))
  
  R1 <- .extract.fd.formulae(F1)
  expect_identical(~edges, R1$form)
  expect_identical(~edges, R1$diss)
  expect_identical(~., R1$nonsep)
  expect_identical(~edges+edges, R1$all)
  
  R2 <- .extract.fd.formulae(F2)
  expect_identical(~edges + gwesp(0,fixed=TRUE), R2$form)
  expect_identical(~edges, R2$diss)
  expect_identical(~., R2$nonsep)
  expect_identical(~edges+gwesp(0,fixed=TRUE)+edges, R2$all)
  
  R3 <- .extract.fd.formulae(F3)
  expect_identical(~edges + gwesp(0,fixed=TRUE), R3$form)
  expect_identical(~edges, R3$diss)
  expect_identical(~edges, R3$nonsep)
  expect_identical(~edges+gwesp(0,fixed=TRUE)+edges+edges, R3$all)
  
  R4 <- .extract.fd.formulae(F4)
  expect_identical(~edges, R4$form)
  expect_identical(~edges, R4$diss)
  expect_identical(~., R4$nonsep)
  expect_identical(~edges+edges, R4$all)
  
  R5 <- .extract.fd.formulae(F5)
  expect_identical(~offset(edges) + gwesp(0,fixed=TRUE), R5$form)
  expect_identical(~offset(edges), R5$diss)
  expect_identical(~., R5$nonsep)
  expect_identical(~offset(edges) + gwesp(0,fixed=TRUE) + offset(edges), R5$all)
  
  R6 <- .extract.fd.formulae(F6)
  expect_identical(~offset(edges), R6$form)
  expect_identical(~-edges, R6$diss)
  expect_identical(~., R6$nonsep)
  expect_identical(~offset(edges) - edges, R6$all)
  
  R7 <- .extract.fd.formulae(F7)
  expect_identical(~offset(edges), R7$form)
  expect_identical(~-offset(edges), R7$diss)
  expect_identical(~., R7$nonsep)
  expect_identical(~offset(edges) - offset(edges), R7$all)
  
  R8 <- .extract.fd.formulae(F8)
  expect_identical(~-edges+triangle, R8$form)
  expect_identical(~edges+triangle-offset(triangle), R8$diss)
  expect_identical(~., R8$nonsep)
  expect_identical(~-edges+triangle+edges+triangle-offset(triangle), R8$all)
  
  R9 <- .extract.fd.formulae(F9)
  expect_identical(~-offset(edges)+offset(triangle)-offset(gwesp(0,fixed=TRUE)), R9$form)
  expect_identical(~edges, R9$diss)
  expect_identical(~., R9$nonsep)
  expect_identical(~-offset(edges)+offset(triangle)-offset(gwesp(0,fixed=TRUE))+edges, R9$all)
  
  R10 <- .extract.fd.formulae(F10)
  expect_identical(~edges-triangle-gwesp(0,fixed=TRUE), R10$form)
  expect_identical(~edges+triangle-gwesp(0,fixed=TRUE), R10$diss)
  expect_identical(~., R10$nonsep)
  expect_identical(~edges-triangle-gwesp(0,fixed=TRUE)+edges+triangle-gwesp(0,fixed=TRUE), R10$all)
})
