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
  
  F1 <- ~Form(~edges) + Diss(~edges)
  
  F2 <- ~Form(~edges + gwesp(0,fixed=TRUE)) + Diss(~edges)
  
  F3 <- ~Form(~edges + gwesp(0,fixed=TRUE)) + Diss(~edges) + edges - triangle
  
  X <- ~edges
  
  F4 <- ~Form(X) + Diss(X)
  
  F5 <- ~Form(~offset(edges) + gwesp(0,fixed=TRUE)) + Diss(~offset(edges))

  Y <- ~offset(edges)
  
  F6 <- ~Form(Y) - Diss(~edges)
  
  F7 <- ~Form(Y) - Diss(Y)
  
  F8 <- ~-Form(~edges - triangle) + Diss(~edges) + Diss(~triangle) + Diss(~-offset(triangle))
  
  F9 <- ~-Form(~offset(edges) - offset(triangle) + offset(gwesp(0,fixed=TRUE))) + Diss(~edges)
  
  F10 <- ~-Form(~-edges + triangle) + Form(~-gwesp(0,fixed=TRUE)) - Diss(~-edges) + Diss(~triangle) + Diss(~-gwesp(0,fixed=TRUE))
  
  R1 <- .extract.fd.formulae(F1)
  expect_equal(~Sum(~edges,label=I), R1$form)
  expect_equal(~Sum(~edges,label=I), R1$diss)
  expect_equal(~., R1$nonsep)
  expect_equal(~Sum(~edges,label=I)+Sum(~edges,label=I), R1$all)
  
  R2 <- .extract.fd.formulae(F2)
  expect_equal(~Sum(~edges + gwesp(0,fixed=TRUE),label=I), R2$form)
  expect_equal(~Sum(~edges,label=I), R2$diss)
  expect_equal(~., R2$nonsep)
  expect_equal(~Sum(~edges+gwesp(0,fixed=TRUE),label=I)+Sum(~edges,label=I), R2$all)
  
  R3 <- .extract.fd.formulae(F3)
  expect_equal(~Sum(~edges + gwesp(0,fixed=TRUE),label=I), R3$form)
  expect_equal(~Sum(~edges,label=I), R3$diss)
  expect_equal(~edges - triangle, R3$nonsep)
  expect_equal(~Sum(~edges+gwesp(0,fixed=TRUE),label=I)+Sum(~edges,label=I)+edges-triangle, R3$all)
  
  R4 <- .extract.fd.formulae(F4)
  expect_equal(~Sum(X,label=I), R4$form)
  expect_equal(~Sum(X,label=I), R4$diss)
  expect_equal(~., R4$nonsep)
  expect_equal(~Sum(X,label=I)+Sum(X,label=I), R4$all)
  
  R5 <- .extract.fd.formulae(F5)
  expect_equal(~Sum(~offset(edges) + gwesp(0,fixed=TRUE),label=I), R5$form)
  expect_equal(~Sum(~offset(edges),label=I), R5$diss)
  expect_equal(~., R5$nonsep)
  expect_equal(~Sum(~offset(edges) + gwesp(0,fixed=TRUE),label=I) + Sum(~offset(edges),label=I), R5$all)
  
  R6 <- .extract.fd.formulae(F6)
  expect_equal(~Sum(Y,label=I), R6$form)
  expect_equal(~-Sum(~edges,label=I), R6$diss)
  expect_equal(~., R6$nonsep)
  expect_equal(~Sum(Y,label=I) - Sum(~edges,label=I), R6$all)
  
  R7 <- .extract.fd.formulae(F7)
  expect_equal(~Sum(Y,label=I), R7$form)
  expect_equal(~-Sum(Y,label=I), R7$diss)
  expect_equal(~., R7$nonsep)
  expect_equal(~Sum(Y,label=I) - Sum(Y,label=I), R7$all)
  
  R8 <- .extract.fd.formulae(F8)
  expect_equal(~-Sum(~edges-triangle,label=I), R8$form)
  expect_equal(~Sum(~edges,label=I)+Sum(~triangle,label=I)+Sum(~-offset(triangle),label=I), R8$diss)
  expect_equal(~., R8$nonsep)
  expect_equal(~-Sum(~edges-triangle,label=I)+Sum(~edges,label=I)+Sum(~triangle,label=I)+Sum(~-offset(triangle),label=I), R8$all)
  
  R9 <- .extract.fd.formulae(F9)
  expect_equal(~-Sum(~offset(edges)-offset(triangle)+offset(gwesp(0,fixed=TRUE)), label=I), R9$form)
  expect_equal(~Sum(~edges,label=I), R9$diss)
  expect_equal(~., R9$nonsep)
  expect_equal(~-Sum(~offset(edges)-offset(triangle)+offset(gwesp(0,fixed=TRUE)),label=I)+Sum(~edges,label=I), R9$all)
  
  R10 <- .extract.fd.formulae(F10)
  expect_equal(~-Sum(~-edges+triangle,label=I)+Sum(~-gwesp(0,fixed=TRUE),label=I), R10$form)
  expect_equal(~-Sum(~-edges,label=I)+Sum(~triangle,label=I)+Sum(~-gwesp(0,fixed=TRUE),label=I), R10$diss)
  expect_equal(~., R10$nonsep)
  expect_equal(~-Sum(~-edges+triangle,label=I)+Sum(~-gwesp(0,fixed=TRUE),label=I)-Sum(~-edges,label=I)+Sum(~triangle,label=I)+Sum(~-gwesp(0,fixed=TRUE),label=I), R10$all)
  
  F11 <- ~Form(~edges) + Diss(~edges) + edges
  environment(F11) <- new.env()
  
  R11 <- .extract.fd.formulae(F11)
  
  expect_identical(environment(F11), environment(R11$form))
  expect_identical(environment(F11), environment(R11$diss))
  expect_identical(environment(F11), environment(R11$nonsep))
  expect_identical(environment(F11), environment(R11$all))
  
  a <- ~edges
  b <- ~triangle
  environment(a) <- new.env()
  environment(b) <- new.env()
  
  F12 <- ~Form(a) + Diss(b) + edges
  environment(F12) <- new.env()
  environment(F12)$a <- a
  environment(F12)$b <- b
  
  R12 <- .extract.fd.formulae(F12)
  
  expect_identical(environment(F12), environment(R12$form))
  expect_identical(environment(F12), environment(R12$diss))
  expect_identical(environment(F12), environment(R12$nonsep))
  expect_identical(environment(F12), environment(R12$all))
  expect_identical(environment(a), environment(eval(R12$form[[2]][[2]], env = environment(R12$form))))
  expect_identical(environment(b), environment(eval(R12$diss[[2]][[2]], env = environment(R12$form))))
  
  
  
  ## diss tests
  F13 <- ~Form(~edges) + Diss(~edges)
  R13 <- .extract.fd.formulae(F13)
  expect_equal(~Sum(~edges,label=I), R13$form)
  expect_equal(~Sum(~edges,label=I), R13$diss)
  expect_equal(~., R13$nonsep)
  expect_equal(~Sum(~edges,label=I)+Sum(~edges,label=I), R13$all)
  
  ## test diss environment
  x <- ~edges
  environment(x) <- new.env()
  F14 <- ~Form(~edges) + Diss(x)
  R14 <- .extract.fd.formulae(F14)
  expect_equal(~Sum(~edges,label=I), R14$form)
  expect_equal(~Sum(x,label=I), R14$diss)
  expect_equal(~., R14$nonsep)
  expect_equal(~Sum(~edges,label=I)+Sum(x,label=I), R14$all)
  expect_identical(environment(x), environment(eval(R14$diss[[2]][[2]], env = environment(R14$diss))))

  ## combine diss
  F15 <- ~Form(~edges) + Diss(~edges) + Diss(~triangle)
  R15 <- .extract.fd.formulae(F15)
  expect_equal(~Sum(~edges,label=I), R15$form)
  expect_equal(~Sum(~edges,label=I)+Sum(~triangle,label=I), R15$diss)
  expect_equal(~., R15$nonsep)
  expect_equal(~Sum(~edges,label=I)+Sum(~edges,label=I)+Sum(~triangle,label=I), R15$all)
    
  ## make some ergm_model calls that would fail if environments were not handled correctly (in particular for diss)
  make_formula <- function() {
    aaaaa <- 52
    ~nodefactor(~0, levels = I(aaaaa))
  }
  make_formula_2 <- function() {
    bbbbb <- 2
    ~nodematch(~3:6, levels = bbbbb, diff = TRUE)
  }
  ff <- make_formula()
  ff2 <- make_formula_2()
  F16 <- ~Form(~edges) + Diss(ff2) + Diss(ff)
  R16 <- .extract.fd.formulae(F16)
  expect_equal(~Sum(~edges,label=I), R16$form)
  expect_equal(~Sum(ff2,label=I)+Sum(ff,label=I), R16$diss)
  expect_equal(~., R16$nonsep)
  expect_equal(~Sum(~edges,label=I)+Sum(ff2,label=I)+Sum(ff,label=I), R16$all)
  
  nw <- network.initialize(8, dir = FALSE)
  m <- ergm_model(R16$diss, nw = nw)
  expect_identical(param_names(m), c("nodematch.3:6.4", "nodefactor.0.52"))
  m2 <- ergm_model(R16$all, nw = nw)  
  expect_identical(param_names(m2), c("edges", "nodematch.3:6.4", "nodefactor.0.52"))
})
