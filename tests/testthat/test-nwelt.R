#  File tests/testthat/test-nwelt.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################

test_that("network.extract.with.lasttoggle behaves reasonably", {

  base_net <- network.initialize(10, dir = FALSE)
  
  # add an edge with no toggle information;
  # it should not show up in any of the
  # extracted lasttoggle matrices
  base_net[8,9] <- 1
  
  # these edges should have lasttoggle information
  # when it is extracted at appropriate timepoints
  edge_toggles <- matrix(c(1L, 1L, 2L,
                           1L, 1L, 3L,
                           1L, 2L, 3L,
                           2L, 3L, 4L,
                           2L, 3L, 5L,
                           2L, 4L, 5L,
                           3L, 1L, 3L,
                           3L, 3L, 4L,
                           3L, 1L, 4L), ncol = 3, byrow = TRUE)
                          
  nwd <- networkDynamic(base.net = base_net, edge.toggles = edge_toggles)
  
  lt_0 <- matrix(0L, ncol = 3, nrow = 0, byrow = TRUE)
  
  lt_1 <- matrix(c(1L, 2L, 1L, 
                   1L, 3L, 1L, 
                   2L, 3L, 1L), ncol = 3, byrow = TRUE)
  
  lt_2 <- matrix(c(1L, 2L, 1L, 
                   1L, 3L, 1L, 
                   2L, 3L, 1L, 
                   3L, 4L, 2L, 
                   3L, 5L, 2L, 
                   4L, 5L, 2L), ncol = 3, byrow = TRUE)
  
  lt_3 <- matrix(c(1L, 2L, 1L, 
                   2L, 3L, 1L, 
                   3L, 5L, 2L, 
                   4L, 5L, 2L, 
                   1L, 3L, 3L, 
                   3L, 4L, 3L, 
                   1L, 4L, 3L), ncol = 3, byrow = TRUE)
  
  lt_4 <- matrix(c(1L, 2L, 1L, 
                   2L, 3L, 1L, 
                   3L, 5L, 2L, 
                   4L, 5L, 2L, 
                   1L, 4L, 3L), ncol = 3, byrow = TRUE)
  
  lt_list <- list(lt_0, lt_1, lt_2, lt_3, lt_4)
  
  for(timeslice in 0:4) {
    nwe <- network.extract.with.lasttoggle(nwd, timeslice)
    expect_identical(nwe %n% "time", timeslice)
    
    extracted_lt <- nwe %n% "lasttoggle"
    manual_lt <- lt_list[[timeslice + 1]]
    
    extracted_lt <- extracted_lt[order(extracted_lt[,1], extracted_lt[,2]),,drop=FALSE]
    manual_lt <- manual_lt[order(manual_lt[,1], manual_lt[,2]),,drop=FALSE]
    
    expect_identical(extracted_lt, manual_lt)
  }
  
  
  
  ## now test deactivating a vertex, ensuring other vertex ids get remapped appropriately
  nwd2 <- deactivate.vertices(nwd, onset=4,terminus=Inf,v=c(1), deactivate.edges=TRUE)
  
  ## we have deleted edges incident on vertex 1, and then decremented other vertex indices by 1
  ## so e.g. the vertex "1L" in the matrix below is actually the second vertex in the original 10 node network
  lt_4_mod <- matrix(c(1L, 2L, 1L, 
                       2L, 4L, 2L, 
                       3L, 4L, 2L), ncol = 3, byrow = TRUE)
  
  # should be same as before up and including time 3
  lt_list_mod <- lt_list
  lt_list_mod[[5]] <- lt_4_mod
  
  for(timeslice in 0:4) {
    nwe <- network.extract.with.lasttoggle(nwd2, timeslice)
    expect_identical(nwe %n% "time", timeslice)
    
    extracted_lt <- nwe %n% "lasttoggle"
    manual_lt <- lt_list_mod[[timeslice + 1]]
    
    extracted_lt <- extracted_lt[order(extracted_lt[,1], extracted_lt[,2]),,drop=FALSE]
    manual_lt <- manual_lt[order(manual_lt[,1], manual_lt[,2]),,drop=FALSE]
    
    expect_identical(extracted_lt, manual_lt)
  }
  
  # now add vertex 1 back in; since we deactivated its edges above, they should
  # not show up in lasttoggle extractions, even after this reactivate.vertices
  # call, but other vertex ids should shift by +1, back to their original values
  nwd3 <- activate.vertices(nwd2, v=c(1), onset=7, terminus=Inf)
  
  lt_10_mod_3 <- matrix(c(2L, 3L, 1L, 
                         3L, 5L, 2L, 
                         4L, 5L, 2L), ncol = 3, byrow = TRUE)
                         
  nwe <- network.extract.with.lasttoggle(nwd3, at=10)
  
  expect_identical(nwe %n% "time", 10)
  expect_identical(nwe %n% "lasttoggle", lt_10_mod_3)
  
  
  ## now try deactivating/reactivating vertex 1 while keeping its edges
  ## (in a couple of different ways)
  
  ## get a fresh one
  nwd <- networkDynamic(base.net = base_net, edge.toggles = edge_toggles)
  
  nwd4 <- deactivate.vertices(nwd, onset=4,terminus=7,v=c(1), deactivate.edges=FALSE)
  
  nwe <- network.extract.with.lasttoggle(nwd4, at=10)
  
  expect_identical(nwe %n% "time", 10)
  expect_identical(nwe %n% "lasttoggle", lt_4)
  
  
  ## get a fresh one
  nwd <- networkDynamic(base.net = base_net, edge.toggles = edge_toggles)
  
  nwd5 <- deactivate.vertices(nwd, onset=4,terminus=Inf,v=c(1), deactivate.edges=FALSE)
  nwd6 <- activate.vertices(nwd, onset=7,terminus=Inf,v=c(1))
  
  nwe <- network.extract.with.lasttoggle(nwd6, at=10)
  
  expect_identical(nwe %n% "time", 10)
  expect_identical(nwe %n% "lasttoggle", lt_4)
})

test_that("network.extract.with.lasttoggle handles both edges and non-edges appropriately", {
  nw <- network.initialize(5, dir=FALSE)
  
  edge_toggles <- matrix(c(1L, 1L, 3L,
                           1L, 1L, 4L,
                           2L, 4L, 5L,
                           3L, 1L, 3L,
                           4L, 1L, 4L), ncol = 3, byrow = TRUE)
                           
  nwd <- networkDynamic(base.net = nw, edge.toggles = edge_toggles)
  deactivate.vertices(nwd, onset = 2, terminus = Inf, v = c(2), deactivate.edges = TRUE)
  deactivate.vertices(nwd, onset = 4, terminus = Inf, v = c(4), deactivate.edges = TRUE)
  
  nw0 <- network.extract.with.lasttoggle(nwd, at = 0)
  nw1 <- network.extract.with.lasttoggle(nwd, at = 1)
  nw2 <- network.extract.with.lasttoggle(nwd, at = 2)
  nw3 <- network.extract.with.lasttoggle(nwd, at = 3)
  nw4 <- network.extract.with.lasttoggle(nwd, at = 4)
  nw5 <- network.extract.with.lasttoggle(nwd, at = 5)
  
  expect_true(network.edgecount(nw0) == 0)
  expect_true(network.edgecount(nw1) == 2 && nw1[1,3] && nw1[1,4])
  expect_true(network.edgecount(nw2) == 3 && nw2[1,2] && nw2[1,3] && nw2[3,4])
  expect_true(network.edgecount(nw3) == 2 && nw3[1,3] && nw3[3,4])
  expect_true(network.edgecount(nw4) == 0)
  expect_true(network.edgecount(nw5) == 0)
  
  expect_identical(nw0 %n% "lasttoggle", matrix(0L, nrow = 0L, ncol = 3L))
  expect_identical(nw1 %n% "lasttoggle", matrix(c(1L, 3L, 1L, 1L, 4L, 1L), nrow = 2L, ncol = 3L, byrow = TRUE))
  expect_identical(nw2 %n% "lasttoggle", matrix(c(1L, 2L, 1L, 1L, 3L, 1L, 3L, 4L, 2L), nrow = 3L, ncol = 3L, byrow = TRUE))
  expect_identical(nw3 %n% "lasttoggle", matrix(c(1L, 3L, 1L, 3L, 4L, 2L, 1L, 2L, 3L), nrow = 3L, ncol = 3L, byrow = TRUE))
  expect_identical(nw4 %n% "lasttoggle", matrix(0L, nrow = 0L, ncol = 3L))
  expect_identical(nw5 %n% "lasttoggle", matrix(0L, nrow = 0L, ncol = 3L))
})
