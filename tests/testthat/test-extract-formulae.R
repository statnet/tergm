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
  expect_equal(~Sum(~edges,label=I), R1$form)
  expect_equal(~Sum(~edges,label=I), R1$pers)
  expect_equal(~., R1$nonsep)
  expect_equal(~Sum(~edges,label=I)+Sum(~edges,label=I), R1$all)
  
  R2 <- .extract.fd.formulae(F2)
  expect_equal(~Sum(~edges + gwesp(0,fixed=TRUE),label=I), R2$form)
  expect_equal(~Sum(~edges,label=I), R2$pers)
  expect_equal(~., R2$nonsep)
  expect_equal(~Sum(~edges+gwesp(0,fixed=TRUE),label=I)+Sum(~edges,label=I), R2$all)
  
  R3 <- .extract.fd.formulae(F3)
  expect_equal(~Sum(~edges + gwesp(0,fixed=TRUE),label=I), R3$form)
  expect_equal(~Sum(~edges,label=I), R3$pers)
  expect_equal(~edges - triangle, R3$nonsep)
  expect_equal(~Sum(~edges+gwesp(0,fixed=TRUE),label=I)+Sum(~edges,label=I)+edges-triangle, R3$all)
  
  R4 <- .extract.fd.formulae(F4)
  expect_equal(~Sum(~edges,label=I), R4$form)
  expect_equal(~Sum(~edges,label=I), R4$pers)
  expect_equal(~., R4$nonsep)
  expect_equal(~Sum(~edges,label=I)+Sum(~edges,label=I), R4$all)
  
  R5 <- .extract.fd.formulae(F5)
  expect_equal(~Sum(~offset(edges) + gwesp(0,fixed=TRUE),label=I), R5$form)
  expect_equal(~Sum(~offset(edges),label=I), R5$pers)
  expect_equal(~., R5$nonsep)
  expect_equal(~Sum(~offset(edges) + gwesp(0,fixed=TRUE),label=I) + Sum(~offset(edges),label=I), R5$all)
  
  R6 <- .extract.fd.formulae(F6)
  expect_equal(~Sum(~offset(edges),label=I), R6$form)
  expect_equal(~-Sum(~edges,label=I), R6$pers)
  expect_equal(~., R6$nonsep)
  expect_equal(~Sum(~offset(edges),label=I) - Sum(~edges,label=I), R6$all)
  
  R7 <- .extract.fd.formulae(F7)
  expect_equal(~Sum(~offset(edges),label=I), R7$form)
  expect_equal(~-Sum(~offset(edges),label=I), R7$pers)
  expect_equal(~., R7$nonsep)
  expect_equal(~Sum(~offset(edges),label=I) - Sum(~offset(edges),label=I), R7$all)
  
  R8 <- .extract.fd.formulae(F8)
  expect_equal(~-Sum(~edges-triangle,label=I), R8$form)
  expect_equal(~Sum(~edges,label=I)+Sum(~triangle,label=I)+Sum(~-offset(triangle),label=I), R8$pers)
  expect_equal(~., R8$nonsep)
  expect_equal(~-Sum(~edges-triangle,label=I)+Sum(~edges,label=I)+Sum(~triangle,label=I)+Sum(~-offset(triangle),label=I), R8$all)
  
  R9 <- .extract.fd.formulae(F9)
  expect_equal(~-Sum(~offset(edges)-offset(triangle)+offset(gwesp(0,fixed=TRUE)), label=I), R9$form)
  expect_equal(~Sum(~edges,label=I), R9$pers)
  expect_equal(~., R9$nonsep)
  expect_equal(~-Sum(~offset(edges)-offset(triangle)+offset(gwesp(0,fixed=TRUE)),label=I)+Sum(~edges,label=I), R9$all)
  
  R10 <- .extract.fd.formulae(F10)
  expect_equal(~-Sum(~-edges+triangle,label=I)+Sum(~-gwesp(0,fixed=TRUE),label=I), R10$form)
  expect_equal(~-Sum(~-edges,label=I)+Sum(~triangle,label=I)+Sum(~-gwesp(0,fixed=TRUE),label=I), R10$pers)
  expect_equal(~., R10$nonsep)
  expect_equal(~-Sum(~-edges+triangle,label=I)+Sum(~-gwesp(0,fixed=TRUE),label=I)-Sum(~-edges,label=I)+Sum(~triangle,label=I)+Sum(~-gwesp(0,fixed=TRUE),label=I), R10$all)
  
  F11 <- ~Form(~edges) + Persist(~edges) + edges
  environment(F11) <- new.env()
  
  R11 <- .extract.fd.formulae(F11)
  
  expect_identical(environment(F11), environment(R11$form))
  expect_identical(environment(F11), environment(R11$pers))
  expect_identical(environment(F11), environment(R11$nonsep))
  expect_identical(environment(F11), environment(R11$all))
  
  a <- ~edges
  b <- ~triangle
  environment(a) <- new.env()
  environment(b) <- new.env()
  
  F12 <- ~Form(a) + Persist(b) + edges
  environment(F12) <- new.env()
  environment(F12)$a <- a
  environment(F12)$b <- b
  
  R12 <- .extract.fd.formulae(F12)
  
  expect_identical(environment(F12), environment(R12$form))
  expect_identical(environment(F12), environment(R12$pers))
  expect_identical(environment(F12), environment(R12$nonsep))
  expect_identical(environment(F12), environment(R12$all))
  expect_identical(environment(a), environment(eval(R12$form[[2]][[2]], env = environment(R12$form))))
  expect_identical(environment(b), environment(eval(R12$pers[[2]][[2]], env = environment(R12$form))))
    
  F13 <- ~Form(~edges) + Persist(~edges)
  R13 <- .extract.fd.formulae(F13)
  expect_equal(~Sum(~edges,label=I), R13$form)
  expect_equal(~Sum(~edges,label=I), R13$pers)
  expect_equal(~., R13$nonsep)
  expect_equal(~Sum(~edges,label=I)+Sum(~edges,label=I), R13$all)
  
  x <- ~edges
  environment(x) <- new.env()
  F14 <- ~Form(~edges) + Persist(x)
  R14 <- .extract.fd.formulae(F14)
  expect_equal(~Sum(~edges,label=I), R14$form)
  expect_equal(~Sum(~edges,label=I), R14$pers)
  expect_equal(~., R14$nonsep)
  expect_equal(~Sum(~edges,label=I)+Sum(~edges,label=I), R14$all)
  expect_identical(environment(x), environment(eval(R14$pers[[2]][[2]], env = environment(R14$pers))))

  F15 <- ~Form(~edges) + Persist(~edges) + Persist(~triangle)
  R15 <- .extract.fd.formulae(F15)
  expect_equal(~Sum(~edges,label=I), R15$form)
  expect_equal(~Sum(~edges,label=I)+Sum(~triangle,label=I), R15$pers)
  expect_equal(~., R15$nonsep)
  expect_equal(~Sum(~edges,label=I)+Sum(~edges,label=I)+Sum(~triangle,label=I), R15$all)
    
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
  F16 <- ~Form(~edges) + Persist(ff2) + Persist(ff)
  R16 <- .extract.fd.formulae(F16)
  expect_equal(~Sum(~edges,label=I), R16$form)
  expect_equal(~Sum(~nodematch(~3:6, levels = bbbbb, diff = TRUE),label=I)+Sum(~nodefactor(~0, levels = I(aaaaa)),label=I), R16$pers)
  expect_equal(~., R16$nonsep)
  expect_equal(~Sum(~edges,label=I)+Sum(~nodematch(~3:6, levels = bbbbb, diff = TRUE),label=I)+Sum(~nodefactor(~0, levels = I(aaaaa)),label=I), R16$all)
  
  nw <- network.initialize(8, dir = FALSE)
  m <- ergm_model(R16$pers, nw = nw)
  expect_identical(param_names(m), c("nodematch.3:6.4", "nodefactor.0.52"))
  m2 <- ergm_model(R16$all, nw = nw)  
  expect_identical(param_names(m2), c("edges", "nodematch.3:6.4", "nodefactor.0.52"))







  ## diss tests
  F17 <- ~Form(~edges) + Persist(~edges) + Diss(~triangle)
  R17 <- .extract.fd.formulae(F17)
  expect_equal(~Sum(~edges,label=I), R17$form)
  expect_equal(~Sum(~edges,label=I)+Sum(-1~triangle,label=I), R17$pers)
  expect_equal(~., R17$nonsep)
  expect_equal(~Sum(~edges,label=I)+Sum(~edges,label=I)+Sum(-1~triangle,label=I), R17$all)
  
  ## test diss environment
  x <- ~edges
  environment(x) <- new.env()
  y <- ~offset(triangle)
  environment(y) <- new.env()  
  F18 <- ~Form(~edges) + Diss(y) + Persist(x)
  R18 <- .extract.fd.formulae(F18)
  expect_equal(~Sum(~edges,label=I), R18$form)
  expect_equal(~Sum(-1~offset(triangle),label=I)+Sum(~edges,label=I), R18$pers)
  expect_equal(~., R18$nonsep)
  expect_equal(~Sum(~edges,label=I)+Sum(-1~offset(triangle),label=I)+Sum(~edges,label=I), R18$all)
  expect_identical(environment(x), environment(eval(R18$pers[[2]][[3]][[2]], env = environment(R18$pers))))
  expect_identical(environment(y), environment(eval(R18$pers[[2]][[2]][[2]], env = environment(R18$pers))))

  ## combine diss
  F19 <- ~Diss(~triangle) + Form(~edges) + gwesp(0, fixed = TRUE) + Persist(~edges) + Diss(~concurrent)
  R19 <- .extract.fd.formulae(F19)
  expect_equal(~Sum(~edges,label=I), R19$form)
  expect_equal(~Sum(-1~triangle,label=I)+Sum(~edges,label=I)+Sum(-1~concurrent,label=I), R19$pers)
  expect_equal(~gwesp(0,fixed=TRUE), R19$nonsep)
  expect_equal(~Sum(-1~triangle,label=I) + Sum(~edges,label=I) + gwesp(0, fixed = TRUE) + Sum(~edges,label=I) + Sum(-1~concurrent,label=I), R19$all)
    
  ## make some ergm_model calls that would fail if environments were not handled correctly (in particular for diss)
  n43 <- "abs"
  make_formula_3 <- function() {
    a1 <- 87
    ~nodefactor(~0:7, levels = I(a1))
  }
  make_formula_4 <- function() {
    b6 <- 4
    ~nodematch(~11:18, levels = b6, diff = TRUE)
  }
  make_formula_5 <- function() {
    n43 <- "wealth"
    ~nodecov(n43)
  }
  make_formula_6 <- function() {
    n43 <- "priorates"
    ~nodecov(n43)
  }
  ff3 <- make_formula_3()
  ff4 <- make_formula_4()
  ff5 <- make_formula_5()
  ff6 <- make_formula_6()
  
  n43 <- "diff"
  
  
  F20 <- ~Persist(ff4) + Diss(ff5) + nodecov(n43) + Form(ff3) + Persist(ff6)
  R20 <- .extract.fd.formulae(F20)
  expect_equal(~Sum(~nodefactor(~0:7, levels = I(a1)),label=I), R20$form)
  expect_equal(~Sum(~nodematch(~11:18, levels = b6, diff = TRUE),label=I)+Sum(-1~nodecov(n43),label=I)+Sum(~nodecov(n43),label=I), R20$pers)
  expect_equal(~nodecov(n43), R20$nonsep)
  expect_equal(~Sum(~nodematch(~11:18, levels = b6, diff = TRUE),label=I)+Sum(-1~nodecov(n43),label=I)+nodecov(n43)+Sum(~nodefactor(~0:7, levels = I(a1)),label=I)+Sum(~nodecov(n43),label=I), R20$all)
  
  nw <- network.initialize(8, dir = FALSE)
  nw %v% "abs" <- letters[1:8]
  nw %v% "diff" <- runif(8)
  nw %v% "wealth" <- runif(8)
  nw %v% "priorates" <- sample(1:10, 8, TRUE)

  m <- ergm_model(R20$form, nw = nw)
  expect_identical(param_names(m), c("nodefactor.0:7.87"))  
  m <- ergm_model(R20$pers, nw = nw)
  expect_identical(param_names(m), c("nodematch.11:18.14", "nodecov.wealth", "nodecov.priorates"))
  m <- ergm_model(R20$nonsep, nw = nw)
  expect_identical(param_names(m), c("nodecov.diff"))
  m <- ergm_model(R20$all, nw = nw)
  expect_identical(param_names(m), c("nodematch.11:18.14", "nodecov.wealth", "nodecov.diff", "nodefactor.0:7.87", "nodecov.priorates"))

})
