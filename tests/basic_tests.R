# These are tests of basic functionality. They will eventually be merged into others.

# Construct a random network
library(tergm)
library(tibble)

nw <- structure(list(mel = list(list(inl = 1L, outl = 4L, atl = list( na = FALSE)), list(inl = 1L, outl = 5L, atl = list(na = FALSE)), list(inl = 1L, outl = 9L, atl = list(na = FALSE)), list(inl = 2L, outl = 4L, atl = list(na = FALSE)), list(inl = 2L, outl = 6L, atl = list(na = FALSE)), list(inl = 2L, outl = 7L, atl = list( na = FALSE)), list(inl = 2L, outl = 9L, atl = list(na = FALSE)), list(inl = 3L, outl = 6L, atl = list(na = FALSE)), list(inl = 3L, outl = 9L, atl = list(na = FALSE)), list(inl = 4L, outl = 5L, atl = list(na = FALSE)), list(inl = 5L, outl = 6L, atl = list( na = FALSE)), list(inl = 5L, outl = 7L, atl = list(na = FALSE)), list(inl = 5L, outl = 10L, atl = list(na = FALSE)), list( inl = 6L, outl = 7L, atl = list(na = FALSE)), list(inl = 6L, outl = 8L, atl = list(na = FALSE)), list(inl = 9L, outl = 10L, atl = list(na = FALSE))), gal = list(n = 10, mnext = 17L, directed = FALSE, hyper = FALSE, loops = FALSE, multiple = FALSE, bipartite = FALSE, lasttoggle = structure(c(1, 1, 1, 2, 2, 2, 2, 3, 3, 4, 5, 5, 5, 6, 6, 9, 3, 4, 5, 9, 4, 6, 7, 9, 6, 9, 5, 6, 7, 10, 7, 8, 10, 7, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1), .Dim = c(17L, 3L)), time = 1), val = list( list(na = FALSE, vertex.names = "1"), list(na = FALSE, vertex.names = "2"), list(na = FALSE, vertex.names = "3"), list(na = FALSE, vertex.names = "4"), list(na = FALSE, vertex.names = "5"), list(na = FALSE, vertex.names = "6"), list(na = FALSE, vertex.names = "7"), list(na = FALSE, vertex.names = "8"), list(na = FALSE, vertex.names = "9"), list(na = FALSE, vertex.names = "10")), iel = list(3:1, 7:4, 9:8, 10L, 13:11, 15:14, integer(0), integer(0), 16L, integer(0)), oel = list(integer(0), integer(0), integer(0), c(1L, 4L), c(2L, 10L), c(5L, 8L, 11L), c(6L, 12L, 14L), 15L, c(3L, 7L, 9L), c(13L, 16L))), class = "network")

el <- cbind(as.edgelist(nw, output="tibble"), edge=TRUE)
lt <- nw %n% "lasttoggle"
colnames(lt) <- c(".tail", ".head", "lt")
merge(el, lt, all=TRUE) # I.e., lasttoggle without an edge, and 1 edge with an old lasttoggle.

nw %n% "time" <- 1
# Here, there are 16 edges in the network, all but 1 added "today" (so DissE(~edges) is 1), and 1 removed "today" (so FormE(~edges) is 17):
stopifnot(all(summary(nw~edges+FormE(~edges)+DissE(~edges))==c(16,17,1))) 

nw %n% "time" <- 2
# No changes at time 2, so all should be 16.
stopifnot(all(summary(nw~edges+FormE(~edges)+DissE(~edges))==c(16,16,16))) 

# Now, let's start at Time 2, add edge (1,2) at time 4, remove it at time 5, and then delete edge (1,4) at time 7 and re-add it at time 8.
o <- tergm.godfather(nw~edges+FormE(~edges)+DissE(~edges), toggles=rbind(c(4,1,2),c(5,1,2),c(7,1,4),c(8,1,4)),start=2, end=10, end.network=TRUE)
attr(o, "stats") # Statistics are appropriate: note how both formation and dissolution "lag":
stopifnot(all(attr(o,"stats")==structure(c(16,17,16,16,15,16,16,16,
                                           16,17,17,16,16,16,16,16,
                                           16,16,16,16,15,15,16,16),
                                         .Dim=c(8L,3L))))
# Finally, the lasttoggle structure should have been updated, with (2,9) last toggled at 0, (1,4) at 8, and (1,2) at 5:
(lt <- o%n%"lasttoggle")
stopifnot(lt[lt[,1]==2 & lt[,2]==9,3] == 0)
stopifnot(lt[lt[,1]==1 & lt[,2]==4,3] == 8)
stopifnot(lt[lt[,1]==2 & lt[,2]==2,3] == 5)
stopifnot(sum(lt[,3]==1)==15)
