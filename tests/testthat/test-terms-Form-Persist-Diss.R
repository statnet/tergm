#  File tests/testthat/test-terms-Form-Persist-Diss.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2024 Statnet Commons
################################################################################

test_that("terms Form, Diss, Persist behave reasonably in dynamic contexts", {
  nw0 <- network.initialize(100, dir = FALSE)
  set.seed(0)
  nw1 <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "final", time.slices = 1)
  set.seed(0)
  nw <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "final", time.slices = 10)
  set.seed(0)
  nwd <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "networkDynamic", time.slices = 10)

  nw10 <- network.collapse(nwd, at = 10)
  nw9 <- network.collapse(nwd, at = 9)
  
  nwunion <- nw9 + nw10
  nwunion %n% "time" <- 10
  nwunion %n% "lasttoggle" <- cbind(as.edgelist(nwunion), 10)
  
  nwinter <- nw9 & nw10
  nwinter %n% "time" <- 10
  nwinter %n% "lasttoggle" <- cbind(as.edgelist(nwinter), 10)
  
  nwunion01 <- nw0 + nw1
  nwunion01 %n% "time" <- 1
  nwunion01 %n% "lasttoggle" <- cbind(as.edgelist(nwunion01), 1)
  
  nwinter01 <- nw0 # no edges
  
  ## non-durational, non-curved
  ff0 <- ~edges + triangle + gwesp(0, fixed = TRUE) + isolates + concurrent
  
  # Form
  s1 <- summary(ff0, basis = nwunion)
  s2 <- summary(~Form(ff0), basis = nw)
  names(s1) <- paste0("Form~", names(s1))
  expect_equal(s1, s2)
  
  s1 <- summary(ff0, basis = nwunion01)
  s2 <- summary(~Form(ff0), basis = nw1)
  names(s1) <- paste0("Form~", names(s1))
  expect_equal(s1, s2)
  
  set.seed(0)
  s3 <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "stats", monitor = ~Form(ff0), time.slices = 1)
  expect_equal(s1, s3[1,])

  # Persist
  s1 <- summary(ff0, basis = nwinter)
  s2 <- summary(~Persist(ff0), basis = nw)
  names(s1) <- paste0("Persist~", names(s1))
  expect_equal(s1, s2)
  
  s1 <- summary(ff0, basis = nwinter01)
  s2 <- summary(~Persist(ff0), basis = nw1)
  names(s1) <- paste0("Persist~", names(s1))
  expect_equal(s1, s2)
  
  set.seed(0)
  s3 <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "stats", monitor = ~Persist(ff0), time.slices = 1)
  expect_equal(s1, s3[1,])

  # Diss
  s1 <- summary(ff0, basis = nwinter)
  s2 <- summary(~Diss(ff0), basis = nw)
  names(s1) <- paste0("Diss~", names(s1))
  expect_equal(-s1, s2)
  
  s1 <- summary(ff0, basis = nwinter01)
  s2 <- summary(~Diss(ff0), basis = nw1)
  names(s1) <- paste0("Diss~", names(s1))
  expect_equal(-s1, s2)
  
  set.seed(0)
  s3 <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "stats", monitor = ~Diss(ff0), time.slices = 1)
  expect_equal(-s1, s3[1,])

  ## non-durational, curved
  ff1 <- ~edges + triangle + gwesp(0, fixed = TRUE) + gwesp(fixed = FALSE, cutoff = 100) + isolates

  # Form
  s1 <- summary(ff1, basis = nwunion)
  s2 <- summary(~Form(ff1), basis = nw)
  names(s1) <- paste0("Form~", names(s1))
  expect_equal(s1, s2)
  
  s1 <- summary(ff1, basis = nwunion01)
  s2 <- summary(~Form(ff1), basis = nw1)
  names(s1) <- paste0("Form~", names(s1))
  expect_equal(s1, s2)
  
  set.seed(0)
  s3 <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "stats", monitor = ~Form(ff1), time.slices = 1)
  expect_equal(s1, s3[1,])

  # Persist
  s1 <- summary(ff1, basis = nwinter)
  s2 <- summary(~Persist(ff1), basis = nw)
  names(s1) <- paste0("Persist~", names(s1))
  expect_equal(s1, s2)
  
  s1 <- summary(ff1, basis = nwinter01)
  s2 <- summary(~Persist(ff1), basis = nw1)
  names(s1) <- paste0("Persist~", names(s1))
  expect_equal(s1, s2)
  
  set.seed(0)
  s3 <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "stats", monitor = ~Persist(ff1), time.slices = 1)
  expect_equal(s1, s3[1,])

  # Diss
  s1 <- summary(ff1, basis = nwinter)
  s2 <- summary(~Diss(ff1), basis = nw)
  names(s1) <- paste0("Diss~", names(s1))
  expect_equal(-s1, s2)
  
  s1 <- summary(ff1, basis = nwinter01)
  s2 <- summary(~Diss(ff1), basis = nw1)
  names(s1) <- paste0("Diss~", names(s1))
  expect_equal(-s1, s2)
  
  set.seed(0)
  s3 <- simulate(nw0 ~ Form(~edges) + Persist(~edges), coef = c(-5, 2), dynamic = TRUE, output = "stats", monitor = ~Diss(ff1), time.slices = 1)
  expect_equal(-s1, s3[1,])

  set.seed(0)
  sim01 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-5, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ff0,
                    output = "networkDynamic")

  set.seed(0)
  sim02 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-5, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~Form(ff0),
                    output = "stats")
                    
  set.seed(0)
  sim03 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-5, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~Persist(ff0),
                    output = "stats")

  set.seed(0)
  sim04 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-5, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~Diss(ff0),
                    output = "stats")
                    
                    
  for(time_step in 10:19) {
    nw_initial <- network.collapse(sim01, at = time_step)
    nw_final <- network.collapse(sim01, at = time_step + 1)

    nw_union <- nw_initial + nw_final
    nw_inter <- nw_initial & nw_final
        
    summ_stats <- summary(ff0, basis = nw_union)
    names(summ_stats) <- paste0("Form~", names(summ_stats))
    expect_equal(summ_stats, sim02[time_step - 10 + 1,])

    summ_stats <- summary(ff0, basis = nw_inter)
    names(summ_stats) <- paste0("Persist~", names(summ_stats))
    expect_equal(summ_stats, sim03[time_step - 10 + 1,])

    summ_stats <- summary(ff0, basis = nw_inter)
    names(summ_stats) <- paste0("Diss~", names(summ_stats))
    expect_equal(-summ_stats, sim04[time_step - 10 + 1,])
  }
                    
  set.seed(1)
  sim11 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-5, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ff1,
                    output = "networkDynamic")

  set.seed(1)
  sim12 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-5, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~Form(ff1),
                    output = "stats")
                    
  set.seed(1)
  sim13 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-5, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~Persist(ff1),
                    output = "stats")

  set.seed(1)
  sim14 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                    coef = c(-5, 2),
                    time.slices = 10,
                    dynamic = TRUE,
                    monitor = ~Diss(ff1),
                    output = "stats")
                    
                    
  for(time_step in 10:19) {
    nw_initial <- network.collapse(sim11, at = time_step)
    nw_final <- network.collapse(sim11, at = time_step + 1)

    nw_union <- nw_initial + nw_final
    nw_inter <- nw_initial & nw_final
        
    summ_stats <- summary(ff1, basis = nw_union)
    names(summ_stats) <- paste0("Form~", names(summ_stats))
    expect_equal(summ_stats, sim12[time_step - 10 + 1,])

    summ_stats <- summary(ff1, basis = nw_inter)
    names(summ_stats) <- paste0("Persist~", names(summ_stats))
    expect_equal(summ_stats, sim13[time_step - 10 + 1,])

    summ_stats <- summary(ff1, basis = nw_inter)
    names(summ_stats) <- paste0("Diss~", names(summ_stats))
    expect_equal(-summ_stats, sim14[time_step - 10 + 1,])
  }  
})
