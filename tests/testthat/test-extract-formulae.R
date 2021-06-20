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
  expect_equal(~.P(~edges), R1$form)
  expect_equal(~.P(~edges), R1$pers)
  expect_equal(~., R1$nonsep)
  expect_equal(~.P(~edges)+.P(~edges), R1$all)
  
  R2 <- .extract.fd.formulae(F2)
  expect_equal(~.P(~edges + gwesp(0,fixed=TRUE)), R2$form)
  expect_equal(~.P(~edges), R2$pers)
  expect_equal(~., R2$nonsep)
  expect_equal(~.P(~edges+gwesp(0,fixed=TRUE))+.P(~edges), R2$all)
  
  R3 <- .extract.fd.formulae(F3)
  expect_equal(~.P(~edges + gwesp(0,fixed=TRUE)), R3$form)
  expect_equal(~.P(~edges), R3$pers)
  expect_equal(~edges - triangle, R3$nonsep)
  expect_equal(~.P(~edges+gwesp(0,fixed=TRUE))+.P(~edges)+edges-triangle, R3$all)
  
  R4 <- .extract.fd.formulae(F4)
  expect_equal(~.P(X), R4$form)
  expect_equal(~.P(X), R4$pers)
  expect_equal(~., R4$nonsep)
  expect_equal(~.P(X)+.P(X), R4$all)
  
  R5 <- .extract.fd.formulae(F5)
  expect_equal(~.P(~offset(edges) + gwesp(0,fixed=TRUE)), R5$form)
  expect_equal(~.P(~offset(edges)), R5$pers)
  expect_equal(~., R5$nonsep)
  expect_equal(~.P(~offset(edges) + gwesp(0,fixed=TRUE)) + .P(~offset(edges)), R5$all)
  
  R6 <- .extract.fd.formulae(F6)
  expect_equal(~.P(Y), R6$form)
  expect_equal(~-.P(~edges), R6$pers)
  expect_equal(~., R6$nonsep)
  expect_equal(~.P(Y) - .P(~edges), R6$all)
  
  R7 <- .extract.fd.formulae(F7)
  expect_equal(~.P(Y), R7$form)
  expect_equal(~-.P(Y), R7$pers)
  expect_equal(~., R7$nonsep)
  expect_equal(~.P(Y) - .P(Y), R7$all)
  
  R8 <- .extract.fd.formulae(F8)
  expect_equal(~-.P(~edges-triangle), R8$form)
  expect_equal(~.P(~edges)+.P(~triangle)+.P(~-offset(triangle)), R8$pers)
  expect_equal(~., R8$nonsep)
  expect_equal(~-.P(~edges-triangle)+.P(~edges)+.P(~triangle)+.P(~-offset(triangle)), R8$all)
  
  R9 <- .extract.fd.formulae(F9)
  expect_equal(~-.P(~offset(edges)-offset(triangle)+offset(gwesp(0,fixed=TRUE))), R9$form)
  expect_equal(~.P(~edges), R9$pers)
  expect_equal(~., R9$nonsep)
  expect_equal(~-.P(~offset(edges)-offset(triangle)+offset(gwesp(0,fixed=TRUE)))+.P(~edges), R9$all)
  
  R10 <- .extract.fd.formulae(F10)
  expect_equal(~-.P(~-edges+triangle)+.P(~-gwesp(0,fixed=TRUE)), R10$form)
  expect_equal(~-.P(~-edges)+.P(~triangle)+.P(~-gwesp(0,fixed=TRUE)), R10$pers)
  expect_equal(~., R10$nonsep)
  expect_equal(~-.P(~-edges+triangle)+.P(~-gwesp(0,fixed=TRUE))-.P(~-edges)+.P(~triangle)+.P(~-gwesp(0,fixed=TRUE)), R10$all)
  
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
  expect_equal(~.P(~edges), R13$form)
  expect_equal(~.P(~edges), R13$pers)
  expect_equal(~., R13$nonsep)
  expect_equal(~.P(~edges)+.P(~edges), R13$all)
  
  x <- ~edges
  environment(x) <- new.env()
  F14 <- ~Form(~edges) + Persist(x)
  R14 <- .extract.fd.formulae(F14)
  expect_equal(~.P(~edges), R14$form)
  expect_equal(~.P(x), R14$pers)
  expect_equal(~., R14$nonsep)
  expect_equal(~.P(~edges)+.P(x), R14$all)
  expect_identical(environment(x), environment(eval(R14$pers[[2]][[2]], env = environment(R14$pers))))

  F15 <- ~Form(~edges) + Persist(~edges) + Persist(~triangle)
  R15 <- .extract.fd.formulae(F15)
  expect_equal(~.P(~edges), R15$form)
  expect_equal(~.P(~edges)+.P(~triangle), R15$pers)
  expect_equal(~., R15$nonsep)
  expect_equal(~.P(~edges)+.P(~edges)+.P(~triangle), R15$all)
    
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
  expect_equal(~.P(~edges), R16$form)
  expect_equal(~.P(ff2)+.P(ff), R16$pers)
  expect_equal(~., R16$nonsep)
  expect_equal(~.P(~edges)+.P(ff2)+.P(ff), R16$all)
  
  nw <- network.initialize(8, dir = FALSE)
  m <- ergm_model(R16$pers, nw = nw)
  expect_identical(param_names(m), c("nodematch.3:6.4", "nodefactor.0.52"))
  m2 <- ergm_model(R16$all, nw = nw)  
  expect_identical(param_names(m2), c("edges", "nodematch.3:6.4", "nodefactor.0.52"))

  ## diss tests
  F17 <- ~Form(~edges) + Persist(~edges) + Diss(~triangle)
  R17 <- .extract.fd.formulae(F17)
  expect_equal(~.P(~edges), R17$form)
  expect_equal(~.P(~edges)+.M(~triangle), R17$pers)
  expect_equal(~., R17$nonsep)
  expect_equal(~.P(~edges)+.P(~edges)+.M(~triangle), R17$all)
  
  ## test diss environment
  x <- ~edges
  environment(x) <- new.env()
  y <- ~offset(triangle)
  environment(y) <- new.env()  
  F18 <- ~Form(~edges) + Diss(y) + Persist(x)
  R18 <- .extract.fd.formulae(F18)
  expect_equal(~.P(~edges), R18$form)
  expect_equal(~.M(y)+.P(x), R18$pers)
  expect_equal(~., R18$nonsep)
  expect_equal(~.P(~edges)+.M(y)+.P(x), R18$all)
  expect_identical(environment(x), environment(eval(R18$pers[[2]][[3]][[2]], env = environment(R18$pers))))
  expect_identical(environment(y), environment(eval(R18$pers[[2]][[2]][[2]], env = environment(R18$pers))))

  ## combine diss
  F19 <- ~Diss(~triangle) + Form(~edges) + gwesp(0, fixed = TRUE) + Persist(~edges) + Diss(~concurrent)
  R19 <- .extract.fd.formulae(F19)
  expect_equal(~.P(~edges), R19$form)
  expect_equal(~.M(~triangle)+.P(~edges)+.M(~concurrent), R19$pers)
  expect_equal(~gwesp(0,fixed=TRUE), R19$nonsep)
  expect_equal(~.M(~triangle) + .P(~edges) + gwesp(0, fixed = TRUE) + .P(~edges) + .M(~concurrent), R19$all)
    
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
  expect_equal(~.P(ff3), R20$form)
  expect_equal(~.P(ff4)+.M(ff5)+.P(ff6), R20$pers)
  expect_equal(~nodecov(n43), R20$nonsep)
  expect_equal(~.P(ff4)+.M(ff5)+nodecov(n43)+.P(ff3)+.P(ff6), R20$all)
  
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

  ## test offsets
  F <- ~Form(~edges) + Form(~offset(triangle)) + offset(Form(~gwesp)) + offset(Diss(~cycle(4))) + edges + offset(concurrent) + Persist(~isolates)
  R <- .extract.fd.formulae(F)
  expect_equal(~.P(~edges)+.P(~offset(triangle))+offset(.P(~gwesp)),R$form)
  expect_equal(~offset(.M(~cycle(4)))+.P(~isolates),R$pers)
  expect_equal(~edges+offset(concurrent),R$nonsep)
  expect_equal(~.P(~edges)+.P(~offset(triangle))+offset(.P(~gwesp))+offset(.M(~cycle(4)))+edges+offset(concurrent)+.P(~isolates),R$all)
  
  ## test partial offsets
  F <- ~Form(~edges) + Form(~offset(triangle,3:4)) + offset(Form(~gwesp),c(TRUE,FALSE)) + offset(Diss(~cycle(4)),6) + edges + offset(concurrent) + Persist(~isolates)
  R <- .extract.fd.formulae(F)
  expect_equal(~.P(~edges)+.P(~offset(triangle,3:4))+offset(.P(~gwesp),c(TRUE,FALSE)),R$form)
  expect_equal(~offset(.M(~cycle(4)),6)+.P(~isolates),R$pers)
  expect_equal(~edges+offset(concurrent),R$nonsep)
  expect_equal(~.P(~edges)+.P(~offset(triangle,3:4))+offset(.P(~gwesp),c(TRUE,FALSE))+offset(.M(~cycle(4)),6)+edges+offset(concurrent)+.P(~isolates),R$all)
  
  ## test sticking in a formula-class formula with a non-default environment  
  F <- ~Form(.) + edges + offset(Diss(.)) + Persist(~concurrent)
  u <- ~nodecov(uniquename)
  environment(u) <- new.env()
  environment(u)$uniquename <- "particularattrname"
  v <- ~isolates
  environment(v) <- new.env()
  F[[2]][[2]][[3]][[2]][[2]] <- v
  F[[2]][[2]][[2]][[2]][[2]] <- u
  R <- .extract.fd.formulae(F)
  expect_identical(environment(u), environment(R$form[[2]][[2]]))
  expect_identical(environment(v), environment(R$pers[[2]][[2]][[2]][[2]]))
  expect_equal(~.P(~nodecov(uniquename)), R$form)
  expect_equal(~offset(.M(~isolates))+.P(~concurrent), R$pers)
  expect_equal(~edges, R$nonsep)
  expect_equal(~.P(~nodecov(uniquename)) + edges + offset(.M(~isolates))+.P(~concurrent), R$all)

  nw <- network.initialize(10, dir = FALSE)
  expect_error(m <- ergm_model(R$form, nw = nw), "is/are not valid nodal attribute\\(s\\)")
  nw %v% "particularattrname" <- sample(1:5, 10, TRUE)
  expect_error(m <- ergm_model(R$form, nw = nw), NA)
  expect_identical(param_names(m, canonical = FALSE), "nodecov.particularattrname")
  
  ## test leading signs
  F <- ~edges - Form(~edges)
  R <- .extract.fd.formulae(F)
  expect_equal(~-.P(~edges), R$form)
  expect_equal(~., R$pers)
  expect_equal(~edges, R$nonsep)
  expect_equal(~edges-.P(~edges), R$all)
  
  F <- ~-edges + Form(~triangle) - Form(~edges) - Diss(~concurrent)
  R <- .extract.fd.formulae(F)
  expect_equal(~.P(~triangle)-.P(~edges), R$form)
  expect_equal(~-.M(~concurrent), R$pers)
  expect_equal(~-edges, R$nonsep)
  expect_equal(~-edges+.P(~triangle)-.P(~edges)-.M(~concurrent), R$all)
  
  F <- ~-offset(Form(~edges)) - Persist(~offset(concurrent)) + offset(Diss(~-isolates)) - offset(triangle) + edges
  R <- .extract.fd.formulae(F)
  expect_equal(~-offset(.P(~edges)), R$form)
  expect_equal(~-.P(~offset(concurrent))+offset(.M(~-isolates)), R$pers)
  expect_equal(~-offset(triangle)+edges, R$nonsep)
  expect_equal(~-offset(.P(~edges))-.P(~offset(concurrent))+offset(.M(~-isolates))-offset(triangle)+edges, R$all)

  ## test interactions
  F <- ~Form(~edges + triangle:gwesp(0,fixed=TRUE)) + edges*concurrent + Diss(~mean.age) + triangle + isolates:degree(1)
  R <- .extract.fd.formulae(F)
  expect_equal(~.P(~edges + triangle:gwesp(0,fixed=TRUE)), R$form)
  expect_equal(~.M(~mean.age), R$pers)
  expect_equal(~edges*concurrent + triangle + isolates:degree(1), R$nonsep)
  expect_equal(~.P(~edges + triangle:gwesp(0,fixed=TRUE)) + edges*concurrent + .M(~mean.age) + triangle + isolates:degree(1), R$all)

  ## test that interactions are always treated as non-separable (maybe not ideal, but at least documentable)
  F <- ~edges:Form(~edges) + Form(~edges)*edges + Diss(~edges)*Form(~edges) + Persist(~edges):Diss(~edges) +
        offset(edges):Form(~edges) + offset(edges:Form(~edges)) + edges:offset(Form(~edges))
  R <- .extract.fd.formulae(F)
  expect_equal(~., R$form)
  expect_equal(~., R$pers)
  expect_equal(F, R$nonsep)
  expect_equal(F, R$all)
})

test_that("terms .P/.M behave as plus/minus identity", {
  nw0 <- network.initialize(100, dir = FALSE)
  nw <- simulate(nw0 ~ edges, coef = c(-4), dynamic = TRUE, output = "final")

  ## non-durational, non-curved
  ff0 <- ~edges + triangle + gwesp(0, fixed = TRUE) + isolates + concurrent
  s1 <- summary(ff0, basis = nw)
  s2 <- summary(~.P(ff0), basis = nw)
  s3 <- summary(~.M(ff0), basis = nw)
  expect_identical(s1, s2)
  expect_identical(s1, -s3)
  
  ## non-durational, curved
  ff1 <- ~edges + triangle + gwesp(0, fixed = TRUE) + gwesp(fixed = FALSE, cutoff = 5)
  s1 <- summary(ff1, basis = nw)
  s2 <- summary(~.P(ff1), basis = nw)
  s3 <- summary(~.M(ff1), basis = nw)
  expect_identical(s1, s2)
  expect_identical(s1, -s3)
  
  set.seed(0)
  nw1 <- simulate(nw ~ Form(~edges + triangle) + Diss(~edges + gwesp(0, fixed = TRUE)), coef = c(-5, 0.01, -4, 0.01), time.slices = 10, dynamic = TRUE, output = "final")
  set.seed(0)
  nw2 <- simulate(nw ~ .P(~Form(~edges + triangle) + Diss(~edges + gwesp(0, fixed = TRUE))), coef = c(-5, 0.01, -4, 0.01), time.slices = 10, dynamic = TRUE, output = "final")
  set.seed(0)
  nw3 <- simulate(nw ~ .M(~Form(~edges + triangle) + Diss(~edges + gwesp(0, fixed = TRUE))), coef = c(5, -0.01, 4, -0.01), time.slices = 10, dynamic = TRUE, output = "final")
  attr(nw1, "coef") <- NULL
  attr(nw2, "coef") <- NULL
  attr(nw3, "coef") <- NULL
  expect_identical(nw1, nw2)
  expect_identical(nw1, nw3)
  nw <- nw1

  ## durational, non-curved
  ff2 <- ~edges + triangle + gwesp(0, fixed = TRUE) + isolates + concurrent + edges.ageinterval(1) + edges.ageinterval(4) + edges.ageinterval(15) + mean.age + Form(~edges + triangle + mean.age) + Persist(~edges + concurrent + isolates + edges.ageinterval(2) + edges.ageinterval(6))
  s1 <- summary(ff2, basis = nw)
  s2 <- summary(~.P(ff2), basis = nw)
  s3 <- summary(~.M(ff2), basis = nw)
  expect_identical(s1, s2)
  expect_identical(s1, -s3)
  
  ## durational, curved
  ff3 <- ~edges + triangle + gwesp(0, fixed = TRUE) + gwesp(fixed = FALSE, cutoff = 5) + mean.age + Form(~edges + triangle + mean.age) + Diss(~gwesp(fixed = FALSE, cutoff = 5))
  s1 <- summary(ff3, basis = nw)
  s2 <- summary(~.P(ff3), basis = nw)
  s3 <- summary(~.M(ff3), basis = nw)
  expect_identical(s1, s2)
  expect_identical(s1, -s3)
  
  set.seed(0)
  sim01 <- simulate(nw ~ Form(~edges) + Persist(~edges), 
                    coef = c(-6, 2), 
                    time.slices = 10, 
                    dynamic = TRUE,
                    monitor = ff0,                    
                    output = "stats")

  set.seed(0)
  sim02 <- simulate(nw ~ Form(~edges) + Persist(~edges), 
                    coef = c(-6, 2), 
                    time.slices = 10, 
                    dynamic = TRUE,
                    monitor = ~.P(ff0),                    
                    output = "stats")
                    
  set.seed(0)
  sim03 <- simulate(nw ~ Form(~edges) + Persist(~edges), 
                    coef = c(-6, 2), 
                    time.slices = 10, 
                    dynamic = TRUE,
                    monitor = ~.M(ff0),                    
                    output = "stats")

  expect_identical(sim01, sim02)
  expect_identical(sim01, -sim03)

  set.seed(1)
  sim11 <- simulate(nw ~ Form(~edges) + Persist(~edges), 
                    coef = c(-6, 2), 
                    time.slices = 10, 
                    dynamic = TRUE,
                    monitor = ff1,                    
                    output = "stats")

  set.seed(1)
  sim12 <- simulate(nw ~ Form(~edges) + Persist(~edges), 
                    coef = c(-6, 2), 
                    time.slices = 10, 
                    dynamic = TRUE,
                    monitor = ~.P(ff1),                    
                    output = "stats")
                    
  set.seed(1)
  sim13 <- simulate(nw ~ Form(~edges) + Persist(~edges), 
                    coef = c(-6, 2), 
                    time.slices = 10, 
                    dynamic = TRUE,
                    monitor = ~.M(ff1),                    
                    output = "stats")

  expect_identical(sim11, sim12)
  expect_identical(sim11, -sim13)

  set.seed(2)
  sim21 <- simulate(nw ~ Form(~edges) + Persist(~edges), 
                    coef = c(-6, 2), 
                    time.slices = 10, 
                    dynamic = TRUE,
                    monitor = ff2,                    
                    output = "stats")

  set.seed(2)
  sim22 <- simulate(nw ~ Form(~edges) + Persist(~edges), 
                    coef = c(-6, 2), 
                    time.slices = 10, 
                    dynamic = TRUE,
                    monitor = ~.P(ff2),                    
                    output = "stats")
                    
  set.seed(2)
  sim23 <- simulate(nw ~ Form(~edges) + Persist(~edges), 
                    coef = c(-6, 2), 
                    time.slices = 10, 
                    dynamic = TRUE,
                    monitor = ~.M(ff2),                    
                    output = "stats")

  expect_identical(sim21, sim22)
  expect_identical(sim21, -sim23)

  set.seed(3)
  sim31 <- simulate(nw ~ Form(~edges) + Persist(~edges), 
                    coef = c(-6, 2), 
                    time.slices = 10, 
                    dynamic = TRUE,
                    monitor = ff3,                    
                    output = "stats")

  set.seed(3)
  sim32 <- simulate(nw ~ Form(~edges) + Persist(~edges), 
                    coef = c(-6, 2), 
                    time.slices = 10, 
                    dynamic = TRUE,
                    monitor = ~.P(ff3),                    
                    output = "stats")
                    
  set.seed(3)
  sim33 <- simulate(nw ~ Form(~edges) + Persist(~edges), 
                    coef = c(-6, 2), 
                    time.slices = 10, 
                    dynamic = TRUE,
                    monitor = ~.M(ff3),                    
                    output = "stats")

  expect_identical(sim31, sim32)
  expect_identical(sim31, -sim33)
  
  
  m01 <- ergm_model(ff0, nw = nw)
  m02 <- ergm_model(~.P(ff0), nw = nw)
  m03 <- ergm_model(~.M(ff0), nw = nw)
  expect_true(!is.curved(m01))
  expect_true(!is.curved(m02))
  expect_true(!is.curved(m03))
  expect_true(!is.durational(m01))
  expect_true(!is.durational(m02))
  expect_true(!is.durational(m03))
  expect_identical(param_names(m01, canonical = FALSE), param_names(m02, canonical = FALSE))
  expect_identical(param_names(m01, canonical = TRUE), param_names(m02, canonical = TRUE))
  expect_identical(param_names(m01, canonical = FALSE), param_names(m03, canonical = FALSE))
  expect_identical(param_names(m01, canonical = TRUE), param_names(m03, canonical = TRUE))
  
  m11 <- ergm_model(ff1, nw = nw)
  m12 <- ergm_model(~.P(ff1), nw = nw)
  m13 <- ergm_model(~.M(ff1), nw = nw)
  expect_true(is.curved(m11))
  expect_true(is.curved(m12))
  expect_true(is.curved(m13))
  expect_true(!is.durational(m11))
  expect_true(!is.durational(m12))
  expect_true(!is.durational(m13))
  expect_identical(param_names(m11, canonical = FALSE), param_names(m12, canonical = FALSE))
  expect_identical(param_names(m11, canonical = TRUE), param_names(m12, canonical = TRUE))
  expect_identical(param_names(m11, canonical = FALSE), param_names(m13, canonical = FALSE))
  expect_identical(param_names(m11, canonical = TRUE), param_names(m13, canonical = TRUE))

  m21 <- ergm_model(ff2, nw = nw)
  m22 <- ergm_model(~.P(ff2), nw = nw)
  m23 <- ergm_model(~.M(ff2), nw = nw)
  expect_true(!is.curved(m21))
  expect_true(!is.curved(m22))
  expect_true(!is.curved(m23))
  expect_true(is.durational(m21))
  expect_true(is.durational(m22))
  expect_true(is.durational(m23))
  expect_identical(param_names(m21, canonical = FALSE), param_names(m22, canonical = FALSE))
  expect_identical(param_names(m21, canonical = TRUE), param_names(m22, canonical = TRUE))
  expect_identical(param_names(m21, canonical = FALSE), param_names(m23, canonical = FALSE))
  expect_identical(param_names(m21, canonical = TRUE), param_names(m23, canonical = TRUE))

  m31 <- ergm_model(ff3, nw = nw)
  m32 <- ergm_model(~.P(ff3), nw = nw)
  m33 <- ergm_model(~.M(ff3), nw = nw)
  expect_true(is.curved(m31))
  expect_true(is.curved(m32))
  expect_true(is.curved(m33))
  expect_true(is.durational(m31))
  expect_true(is.durational(m32))
  expect_true(is.durational(m33))
  expect_identical(param_names(m31, canonical = FALSE), param_names(m32, canonical = FALSE))
  expect_identical(param_names(m31, canonical = TRUE), param_names(m32, canonical = TRUE))
  expect_identical(param_names(m31, canonical = FALSE), param_names(m33, canonical = FALSE))
  expect_identical(param_names(m31, canonical = TRUE), param_names(m33, canonical = TRUE))
})
