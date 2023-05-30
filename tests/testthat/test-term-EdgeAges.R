#  File tests/testthat/test-term-EdgeAges.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2023 Statnet Commons
################################################################################

test_that("the EdgeAges term behaves consistently with some existing durational terms", {
  net_size <- 100
  bip_size <- 40
  edges_target <- 100
  logit <- function(p) log(p/(1-p))
  density <- 1/50
  D <- 10
  seed <- 0

  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }

      nw0 <- network.initialize(net_size, directed = directed, bipartite = bipartite)
      nw0 %v% "attr" <- rep(1:3, length.out = net_size)
      nw1 <- san(nw0 ~ edges, target.stats = c(edges_target))
      nw1 %n% "lasttoggle" <- cbind(as.edgelist(nw1), -as.integer(100*runif(network.edgecount(nw1))))
      nw1 %n% "time" <- 0

      for(nw in list(nw0, nw1)) {
        mat <- matrix(runif(net_size*net_size), net_size, net_size)
        if(!directed && !bipartite) mat <- (mat + t(mat))/2

        set.seed(seed)
        s1 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                       coef = c(logit(density) - log(D), log(D - 1)),
                       time.slices = 10,
                       dynamic = TRUE,
                       output = "stats",
                       monitor = ~edge.ages + edgecov.ages(mat) + nodemix("attr", levels2 = TRUE) + nodemix.mean.age("attr", levels2 = TRUE))
        set.seed(seed)
        s2 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                       coef = c(logit(density) - log(D), log(D - 1)),
                       time.slices = 10,
                       dynamic = TRUE,
                       output = "stats",
                       monitor = ~EdgeAges(~edges + edgecov(mat) + nodematch("attr", diff = TRUE) + nodemix("attr", levels2 = TRUE)))

        s1 <- unname(s1)
        s2 <- unname(s2)

        s1edges <- s1[,1,drop=FALSE]
        s1edgecov <- s1[,2,drop=FALSE]
        s1 <- s1[,-(1:2),drop=FALSE]
        s1mix <- s1[,seq_len(NCOL(s1)/2)]*s1[,-seq_len(NCOL(s1)/2)]
        diag_indices <- if(directed || bipartite) c(1,5,9) else c(1,3,6)
        s1match <- s1mix[,diag_indices,drop=FALSE]

        s2edges <- s2[,1,drop=FALSE]
        s2edgecov <- s2[,2,drop=FALSE]
        s2 <- s2[,-(1:2),drop=FALSE]
        s2mix <- s2[,-(1:3),drop=FALSE]
        s2match <- s2[,1:3,drop=FALSE]

        expect_equal(s1edges, s2edges)
        expect_equal(s1edgecov, s2edgecov)
        expect_equal(s1mix, s2mix)
        expect_equal(s1match, s2match)
      }
    }
  }
})

test_that("the EdgeAges term behaves appropriately for general submodels", {
  net_size <- 100
  bip_size <- 40
  edges_target <- 100
  logit <- function(p) log(p/(1-p))
  density <- 1/50
  D <- 10
  seed <- 0

  ff <- ~edges +
         nodematch("attr") +
         nodefactor("attr") +
         nodemix("attr") +
         nodecov("qattr") +
         edgecov(mat)

  ffi <- ~edge.ages:(edges +
                     nodematch("attr") +
                     nodefactor("attr") +
                     nodemix("attr") +
                     nodecov("qattr") +
                     edgecov(mat))

  for(directed in list(FALSE, TRUE)) {
    for(bipartite in list(FALSE, bip_size)) {
      if(directed && bipartite) {
        next
      }

      nw0 <- network.initialize(net_size, directed = directed, bipartite = bipartite)
      nw0 %v% "attr" <- rep(1:3, length.out = net_size)
      nw0 %v% "qattr" <- runif(net_size)
      nw1 <- san(nw0 ~ edges, target.stats = c(edges_target))
      nw1 %n% "lasttoggle" <- cbind(as.edgelist(nw1), -as.integer(100*runif(network.edgecount(nw1))))
      nw1 %n% "time" <- 0

      for(nw in list(nw0, nw1)) {
        mat <- matrix(runif(net_size*net_size), net_size, net_size)
        if(!directed && !bipartite) mat <- (mat + t(mat))/2

        set.seed(seed)
        s1 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                       coef = c(logit(density) - log(D), log(D - 1)),
                       time.slices = 10,
                       dynamic = TRUE,
                       output = "networkDynamic")
        set.seed(seed)
        s2 <- simulate(nw ~ Form(~edges) + Persist(~edges),
                       coef = c(logit(density) - log(D), log(D - 1)),
                       time.slices = 10,
                       dynamic = TRUE,
                       output = "stats",
                       monitor = ~EdgeAges(ff))

        s1 <- summary(ffi, basis = s1, at = 1:10)
        s1 <- unname(s1)
        s2 <- unname(as.matrix(s2))
        expect_equal(s1, s2)
      }
    }
  }
})
