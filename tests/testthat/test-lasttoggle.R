#  File tests/testthat/test-lasttoggle.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2022 Statnet Commons
################################################################################
# These are tests of basic functionality of lasttoggle. They will eventually be merged into others.

# TODO: Figure out which of these are redundant.

test_that("lasttoggle tests pass", {

# Construct a random network
nw <- structure(list(mel = list(list(inl = 1L, outl = 4L, atl = list( na = FALSE)), list(inl = 1L, outl = 5L, atl = list(na = FALSE)), list(inl = 1L, outl = 9L, atl = list(na = FALSE)), list(inl = 2L, outl = 4L, atl = list(na = FALSE)), list(inl = 2L, outl = 6L, atl = list(na = FALSE)), list(inl = 2L, outl = 7L, atl = list( na = FALSE)), list(inl = 2L, outl = 9L, atl = list(na = FALSE)), list(inl = 3L, outl = 6L, atl = list(na = FALSE)), list(inl = 3L, outl = 9L, atl = list(na = FALSE)), list(inl = 4L, outl = 5L, atl = list(na = FALSE)), list(inl = 5L, outl = 6L, atl = list( na = FALSE)), list(inl = 5L, outl = 7L, atl = list(na = FALSE)), list(inl = 5L, outl = 10L, atl = list(na = FALSE)), list( inl = 6L, outl = 7L, atl = list(na = FALSE)), list(inl = 6L, outl = 8L, atl = list(na = FALSE)), list(inl = 9L, outl = 10L, atl = list(na = FALSE))), gal = list(n = 10, mnext = 17L, directed = FALSE, hyper = FALSE, loops = FALSE, multiple = FALSE, bipartite = FALSE, lasttoggle = structure(c(1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5, 5, 5, 6, 6, 9, 3, 4, 5, 9, 4, 6, 7, 9, 6, 9, 5, 6, 7, 10, 7, 8, 10, 7, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), .Dim = c(17L, 3L)), time = 1), val = list( list(na = FALSE, vertex.names = "1"), list(na = FALSE, vertex.names = "2"), list(na = FALSE, vertex.names = "3"), list(na = FALSE, vertex.names = "4"), list(na = FALSE, vertex.names = "5"), list(na = FALSE, vertex.names = "6"), list(na = FALSE, vertex.names = "7"), list(na = FALSE, vertex.names = "8"), list(na = FALSE, vertex.names = "9"), list(na = FALSE, vertex.names = "10")), iel = list(3:1, 7:4, 9:8, 10L, 13:11, 15:14, integer(0), integer(0), 16L, integer(0)), oel = list(integer(0), integer(0), integer(0), c(1L, 4L), c(2L, 10L), c(5L, 8L, 11L), c(6L, 12L, 14L), 15L, c(3L, 7L, 9L), c(13L, 16L))), class = "network")

el <- cbind(as.edgelist(nw, output="tibble"), edge=TRUE)
lt <- nw %n% "lasttoggle"
colnames(lt) <- c(".tail", ".head", "lt")
merge(el, lt, all=TRUE) # I.e., lasttoggle without an edge, and 1 edge with an old lasttoggle.

nw %n% "time" <- 1
# Here, there are 16 edges in the network, all but 1 added "today" (so Persist(~edges) is 1), and 1 removed "today" (so Form(~edges) is 17):
expect_equal(summary(nw~edges+Form(~edges)+Persist(~edges)), c(16,17,1), ignore_attr=TRUE)

nw %n% "time" <- 2
# No changes at time 2, so all should be 16.
expect_equal(summary(nw~edges+Form(~edges)+Persist(~edges)), c(16,16,16), ignore_attr=TRUE)

# Now, let's start at Time 2, add edge (1,2) at time 4, remove it at time 5, and then delete edge (1,4) at time 7 and re-add it at time 8.
o <- tergm.godfather(nw~edges+Form(~edges)+Persist(~edges), toggles=rbind(c(4,1,2),c(5,1,2),c(7,1,4),c(8,1,4)),start=2, end=10, end.network=TRUE)
attr(o, "stats") # Statistics are appropriate: note how both formation and dissolution "lag":
expect_equal(attr(o,"stats"), structure(c(16,17,16,16,15,16,16,16,
                                          16,17,17,16,16,16,16,16,
                                          16,16,16,16,15,15,16,16),
                                        .Dim=c(8L,3L)), ignore_attr=TRUE)

# Finally, the lasttoggle structure should have been updated, with (2,9) last toggled at 0, (1,4) at 8, and (1,2) at 5:
(lt <- o%n%"lasttoggle")
expect_equal(lt[lt[,1]==2 & lt[,2]==9,3], 0)
expect_equal(lt[lt[,1]==1 & lt[,2]==4,3], 8)
expect_equal(lt[lt[,1]==1 & lt[,2]==2,3], 5)
expect_equal(sum(lt[,3]==1), 15)

# Test that auxiliary networks are being tracked properly.
el <- matrix(c(1, 2,
               1, 3,
               1, 4,
               2, 3), ncol=2, nrow=4, byrow=TRUE)

nw2 <- network(el, type="edgelist", dir=FALSE)
nw2 %n% "lasttoggle" <- cbind(el, c(1,1,1,2))
nw2 %n% "time" <- 2
expect_equal(summary(nw2 ~ Form(~edges + triangle) + Persist(~edges + triangle)), c(4,1,3,0), ignore_attr=TRUE)
nw2 %n% "time" <- 3
expect_equal(summary(nw2 ~ Form(~edges + triangle) + Persist(~edges + triangle)), c(4,1,4,1), ignore_attr=TRUE)



## test initialization (via dyad-dependent terms in formation model) and time advancement
nw <- network.initialize(100,dir=F)
nw <- san(nw ~ edges + offset(concurrent), offset.coef = c(-Inf), target.stats=40)

expect_true(summary(nw ~ edges) == 40 && summary(nw ~ concurrent) == 0)

nw %n% "time" <- 0
nw %n% "lasttoggle" <- cbind(as.edgelist(nw), 0)

nw <- simulate(nw ~ Form(~edges + concurrent) + Persist(~edges), coef = c(1,-Inf,1), dynamic = TRUE, output="final")
rv <- summary(nw ~ edges + concurrent)
expect_equal(summary(nw ~ concurrent), 0, ignore_attr=TRUE)
expect_equal(nw %n% "time", 1)
nw <- simulate(nw ~ Form(~edges + concurrent) + Persist(~edges), coef = c(1,-Inf,1), dynamic = TRUE, output="final")
rv <- summary(nw ~ edges + concurrent)
expect_equal(summary(nw ~ concurrent), 0, ignore_attr=TRUE)
expect_equal(nw %n% "time", 2)
nw <- simulate(nw ~ Form(~edges + concurrent) + Persist(~edges), coef = c(1,-Inf,1), dynamic = TRUE, output="final")
rv <- summary(nw ~ edges + concurrent)
expect_equal(summary(nw ~ concurrent), 0, ignore_attr=TRUE)
expect_equal(nw %n% "time", 3)
nw <- simulate(nw ~ Form(~edges + concurrent) + Persist(~edges), coef = c(1,-Inf,1), dynamic = TRUE, output="final")
rv <- summary(nw ~ edges + concurrent)
expect_equal(summary(nw ~ concurrent), 0, ignore_attr=TRUE)
expect_equal(nw %n% "time", 4)
nw <- simulate(nw ~ Form(~edges + concurrent) + Persist(~edges), coef = c(1,-Inf,1), dynamic = TRUE, output="final")
rv <- summary(nw ~ edges + concurrent)
expect_equal(summary(nw ~ concurrent), 0, ignore_attr=TRUE)
expect_equal(nw %n% "time", 5)


nw <- network.initialize(100,dir=F)
nw <- san(nw ~ edges + offset(triangle), offset.coef = c(-Inf), target.stats=40)

expect_true(summary(nw ~ edges) == 40 && summary(nw ~ triangle) == 0)

nw %n% "time" <- 0
nw %n% "lasttoggle" <- cbind(as.edgelist(nw), 0)

nw <- simulate(nw ~ Form(~edges + triangle) + Persist(~edges), coef = c(-1,-Inf,1), dynamic = TRUE, output="final")
rv <- summary(nw ~ edges + triangle)
expect_equal(summary(nw ~ triangle), 0, ignore_attr=TRUE)
expect_equal(nw %n% "time", 1)
nw <- simulate(nw ~ Form(~edges + triangle) + Persist(~edges), coef = c(-1,-Inf,1), dynamic = TRUE, output="final")
rv <- summary(nw ~ edges + triangle)
expect_equal(summary(nw ~ triangle), 0, ignore_attr=TRUE)
expect_equal(nw %n% "time", 2)
nw <- simulate(nw ~ Form(~edges + triangle) + Persist(~edges), coef = c(-1,-Inf,1), dynamic = TRUE, output="final")
rv <- summary(nw ~ edges + triangle)
expect_equal(summary(nw ~ triangle), 0, ignore_attr=TRUE)
expect_equal(nw %n% "time", 3)
nw <- simulate(nw ~ Form(~edges + triangle) + Persist(~edges), coef = c(-1,-Inf,1), dynamic = TRUE, output="final")
rv <- summary(nw ~ edges + triangle)
expect_equal(summary(nw ~ triangle), 0, ignore_attr=TRUE)
expect_equal(nw %n% "time", 4)
nw <- simulate(nw ~ Form(~edges + triangle) + Persist(~edges), coef = c(-1,-Inf,1), dynamic = TRUE, output="final")
rv <- summary(nw ~ edges + triangle)
expect_equal(summary(nw ~ triangle), 0, ignore_attr=TRUE)
expect_equal(nw %n% "time", 5)


## durational operator tests, including emptynwstats
nw <- network.initialize(100, dir = FALSE)

stat_names <- c("Form~edges", "Form~isolates", "Persist~edges", "Persist~isolates", "Diss~edges", "Diss~isolates", "Change~edges", "Change~isolates")

# y0, y1 empty.

expected <- c(0, 100, 0, 100, 0, -100, 0, 100)
names(expected) <- stat_names
actual <- summary(nw ~ Form(~edges+isolates) + Persist(~edges+isolates) + Diss(~edges+isolates) + Change(~edges+isolates), dynamic=TRUE)
expect_identical(actual, expected)

# y0, y1 both have (1,2)
nw[1,2] <- TRUE

expected <- c(1, 98, 1, 98, -1, -98, 0, 100)
names(expected) <- stat_names
actual <- summary(nw ~ Form(~edges+isolates) + Persist(~edges+isolates) + Diss(~edges+isolates) + Change(~edges+isolates), dynamic=TRUE)
expect_identical(actual, expected)

# y1 has (1,2), y0 has (1,2) and (3,4)
nw %n% "time" <- 1
nw %n% "lasttoggle" <- cbind(3,4,1)

expected <- c(2, 96, 1, 98, -1, -98, 1, 98)
names(expected) <- stat_names
actual <- summary(nw ~ Form(~edges+isolates) + Persist(~edges+isolates) + Diss(~edges+isolates) + Change(~edges+isolates), dynamic=TRUE)
expect_identical(actual, expected)

# y1 has (1,2), y0 has (1,2), (2,3), and (3,4)
nw %n% "lasttoggle" <- rbind(nw %n% "lasttoggle", cbind(2,3,1))

expected <- c(3, 96, 1, 98, -1, -98, 2, 97)
names(expected) <- stat_names
actual <- summary(nw ~ Form(~edges+isolates) + Persist(~edges+isolates) + Diss(~edges+isolates) + Change(~edges+isolates), dynamic=TRUE)
expect_identical(actual, expected)
})
