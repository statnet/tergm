#  File tests/testthat/test-networkLite.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2024 Statnet Commons
################################################################################

## run tests conditionally on availability of the networkLite package
if(require("networkLite")) {

  test_that("network and networkLite simulate and summarize formulas equally in tergm", {

    net_size <- 100
    bip_size <- 40

    ffdir <- ~nodemix(~a) + absdiff(~b) + odegrange(2) + idegrange(2) +
              gwesp(cutoff = 100) + mean.age + edge.ages + nodemix.mean.age(~a) +
              gwnsp(0.3, fixed = TRUE)
    ffundir <- ~nodemix(~a) + absdiff(~b) + concurrent + gwesp(cutoff = 100) +
                mean.age + edge.ages + nodemix.mean.age(~a) +
                gwnsp(0.3, fixed = TRUE)

    for(directed in list(FALSE, TRUE)) {
      for(bipartite in list(FALSE, bip_size)) {
        if(directed && bipartite) {
          next
        }
        if (directed) {
          ff <- ffdir
        } else {
          ff <- ffundir
        }

        set.seed(0)
        nw <- network.initialize(net_size, directed = directed, bipartite = bipartite)
        nw %v% "a" <- rep(letters[1:5], length.out = net_size)
        nw %v% "b" <- runif(net_size)

        nwL <- as.networkLite(nw)

        coef <- c(-4, 1, 1.5, 0.5, -1, 0.5, 3)

        set.seed(0)
        nw_1 <- simulate(nw ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) +
                              Persist(~edges),
                         coef = coef, output = "final", dynamic = TRUE)
        set.seed(0)
        nwL_1 <- simulate(nwL ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) +
                                Persist(~edges),
                          coef = coef, output = "final", dynamic = TRUE)
        expect_s3_class(nwL_1, "networkLite")

        expect_equal(as.edgelist(nw_1), as.edgelist(nwL_1))
        expect_identical(nw_1 %n% "lasttoggle", nwL_1 %n% "lasttoggle")
        expect_identical(nw_1 %n% "time", nwL_1 %n% "time")
        expect_identical(summary(ff, basis = nw_1),
                         summary(ff, basis = nwL_1))


        set.seed(0)
        nw_2 <- simulate(nw_1 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) +
                                Persist(~edges),
                         coef = coef, output = "final", dynamic = TRUE)
        set.seed(0)
        nwL_2 <- simulate(nwL_1 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) +
                                  Persist(~edges),
                          coef = coef, output = "final", dynamic = TRUE)
        expect_s3_class(nwL_2, "networkLite")

        expect_equal(as.edgelist(nw_2), as.edgelist(nwL_2))
        expect_identical(nw_2 %n% "lasttoggle", nwL_2 %n% "lasttoggle")
        expect_identical(nw_2 %n% "time", nwL_2 %n% "time")
        expect_identical(summary(ff, basis = nw_2),
                         summary(ff, basis = nwL_2))

        set.seed(0)
        nw_3 <- simulate(nw_2 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) +
                                Persist(~edges),
                         coef = coef, output = "final", dynamic = TRUE)
        set.seed(0)
        nwL_3 <- simulate(nwL_2 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) +
                                  Persist(~edges),
                          coef = coef, output = "final", dynamic = TRUE)
        expect_s3_class(nwL_3, "networkLite")

        expect_equal(as.edgelist(nw_3), as.edgelist(nwL_3))
        expect_identical(nw_3 %n% "lasttoggle", nwL_3 %n% "lasttoggle")
        expect_identical(nw_3 %n% "time", nwL_3 %n% "time")
        expect_identical(summary(ff, basis = nw_3),
                         summary(ff, basis = nwL_3))

        set.seed(0)
        nw_4 <- simulate(nw_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) +
                                Persist(~edges),
                         coef = coef, dynamic = TRUE)
        set.seed(0)
        nwL_4 <- simulate(nwL_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) +
                                  Persist(~edges),
                          coef = coef, dynamic = TRUE)

        # comparison of networkDynamics
        expect_equal(nw_4, nwL_4)


        ## for completeness, also get stats and changes as output
        set.seed(0)
        s <- simulate(nw_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) +
                             Persist(~edges),
                      coef = coef, dynamic = TRUE, output = "stats", stats = TRUE,
                      monitor = if(directed) ~edges + idegree(0:10) + odegree(0:10) +
                                              mean.age + Form(~odegree(0:2))
                                else ~edges + degree(0:10) + mean.age +
                                      Form(~degree(0:2)))
        set.seed(0)
        sL <- simulate(nwL_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) +
                               Persist(~edges),
                       coef = coef, dynamic = TRUE, output = "stats", stats = TRUE,
                       monitor = if(directed) ~edges + idegree(0:10) + odegree(0:10) +
                                               mean.age + Form(~odegree(0:2))
                                 else ~edges + degree(0:10) + mean.age +
                                       Form(~degree(0:2)))

        # comparison of stats
        expect_equal(s, sL)

        set.seed(0)
        c <- simulate(nw_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) +
                             Persist(~edges),
                      coef = coef, dynamic = TRUE, output = "changes")
        set.seed(0)
        cL <- simulate(nwL_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) +
                               Persist(~edges),
                       coef = coef, dynamic = TRUE, output = "changes")

        # comparison of changes
        expect_equal(c, cL)

        # again, without lasttoggle
        nw_3 %n% "lasttoggle" <- NULL
        nwL_3 %n% "lasttoggle" <- NULL

        set.seed(0)
        nw_4 <- simulate(nw_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) +
                                Persist(~edges),
                         coef = coef, dynamic = TRUE)
        set.seed(0)
        nwL_4 <- simulate(nwL_3 ~ Form(~edges + nodefactor("a") + nodecov(~b^2 + b)) +
                                  Persist(~edges),
                          coef = coef, dynamic = TRUE)

        # comparison of networkDynamics
        expect_equal(nw_4, nwL_4)
      }
    }
  })

  test_that("conversions between network, networkLite, and networkDynamic behave as expected with non-atomic attributes", {

    logit <- function(p) log(p/(1-p))
    nw <- network.initialize(100, directed = FALSE)

    ## set some arbitrary non-atomic vertex attribute
    set.vertex.attribute(nw, "vertex_attr",
                         lapply(seq_len(network.size(nw)),
                                function(i) { runif(rbinom(1,10,0.5)) } ))

    nw <- san(nw ~ edges, target.stats = c(100))

    ## set some arbitrary non-atomic edge attribute
    set.edge.attribute(nw, "edge_attr",
                       lapply(seq_len(network.edgecount(nw)),
                              function(i) { list(runif(rbinom(1,10,0.5))) } ))

    ## edge activities will be non-atomic
    nwD <- simulate(nw ~ edges, coef = c(logit(1/10)), dynamic = TRUE, time.slices = 100)

    for (dynamic in list(FALSE, TRUE)) {
      if (dynamic == FALSE) {
        nw_base <- nw
        nwL <- as.networkLite(nw_base)
        nw_rebase <- to_network_networkLite(nwL)
      } else {
        nw_base <- nwD
        nwL <- as.networkLite(nw_base)
        nw_rebase <- as.networkDynamic(nwL)
      }

      expect_identical(as.edgelist(nw_base), as.edgelist(nwL))
      expect_identical(as.edgelist(nwL), as.edgelist(nw_rebase))

      expect_identical(list.vertex.attributes(nw_base), list.vertex.attributes(nwL))
      expect_identical(list.vertex.attributes(nwL), list.vertex.attributes(nw_rebase))

      expect_identical(list.edge.attributes(nw_base), list.edge.attributes(nwL))
      expect_identical(list.edge.attributes(nwL), list.edge.attributes(nw_rebase))

      expect_identical(setdiff(list.network.attributes(nw_base), "mnext"), list.network.attributes(nwL))
      expect_identical(list.network.attributes(nwL), setdiff(list.network.attributes(nw_rebase), "mnext"))

      for (attrname in list.vertex.attributes(nwL)) {
        for (unlist in list(FALSE, TRUE)) {
          expect_identical(get.vertex.attribute(nw_base, attrname, unlist = unlist),
                           get.vertex.attribute(nwL, attrname, unlist = unlist))
          expect_identical(get.vertex.attribute(nwL, attrname, unlist = unlist),
                           get.vertex.attribute(nw_rebase, attrname, unlist = unlist))
        }
      }

      ## need to consistently order edges before comparing edge attributes
      el <- as.edgelist(nwL)
      eidsD <- unlist(get.dyads.eids(nw_base, el[,1], el[,2]))
      eidsLD <- unlist(get.dyads.eids(nw_rebase, el[,1], el[,2]))

      for (attrname in list.edge.attributes(nwL)) {
        eaD <- get.edge.attribute(nw_base, attrname, null.na = FALSE, unlist = FALSE)[eidsD]
        eaL <- get.edge.attribute(nwL, attrname, null.na = FALSE, unlist = FALSE)
        eaLD <- get.edge.attribute(nw_rebase, attrname, null.na = FALSE, unlist = FALSE)[eidsLD]

        expect_identical(eaD, eaL)
        expect_identical(eaL, eaLD)
      }

      for (attrname in list.network.attributes(nwL)) {
        expect_identical(get.network.attribute(nw_base, attrname),
                         get.network.attribute(nwL, attrname))
        expect_identical(get.network.attribute(nwL, attrname),
                         get.network.attribute(nw_rebase, attrname))
      }
    }
  })

}
