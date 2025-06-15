#  File tests/testthat/test-edgelist_with_lasttoggle.R in package tergm, part
#  of the Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################

test_that("edgelist_with_lasttoggle behaves reasonably", {
  for(net_type in c("U", "D", "B")) {
    if(net_type == "U") {
      nw0 <- network.initialize(100, dir = FALSE)
    } else if(net_type == "D") {
      nw0 <- network.initialize(100, dir = TRUE)
    } else {
      nw0 <- network.initialize(100, dir = FALSE, bip = 40)    
    }
    
    nw <- san(nw0 ~ edges, target.stats = 500)
    
    el <- as.edgelist(nw)
    lt <- cbind(el, -rpois(NROW(el), 10))
    
    nw %n% "lasttoggle" <- lt
    
    expect_equal(lt, edgelist_with_lasttoggle(nw))
    
    ## scramble order
    nw %n% "lasttoggle" <- lt[sample(seq_len(NROW(lt))),,drop=FALSE]
    
    expect_equal(lt, edgelist_with_lasttoggle(nw))
    
    ## add some
    nw2 <- san(nw0 ~ edges, target.stats = 500)
    nwdiff <- nw2 - nw
    ltdiff <- cbind(as.edgelist(nwdiff), -rpois(network.edgecount(nwdiff), 10))
    ltadd <- rbind(lt, ltdiff)
    ltadd <- ltadd[sample(seq_len(NROW(ltadd))),,drop=FALSE]
    nw %n% "lasttoggle" <- ltadd
    expect_equal(lt, edgelist_with_lasttoggle(nw))
    
    ## remove some
    rtk <- sample(seq_len(NROW(el)), as.integer(NROW(el)/2), FALSE)
    ltdrop <- lt[rtk,,drop=FALSE]
    ltdrop <- ltdrop[sample(seq_len(NROW(ltdrop))),,drop=FALSE]
    nw %n% "lasttoggle" <- ltdrop
    ltmod <- lt
    ltmod[-rtk,3] <- as.integer(-.Machine$integer.max/2)
    expect_equal(ltmod, edgelist_with_lasttoggle(nw))
    
    ## add and remove some
    ltdropadd <- rbind(ltdrop, ltdiff)
    ltdropadd <- ltdropadd[sample(seq_len(NROW(ltdropadd))),,drop=FALSE]
    nw %n% "lasttoggle" <- ltdropadd
    expect_equal(ltmod, edgelist_with_lasttoggle(nw))
    
    ## test corner cases
    expect_equal(matrix(0L, 0, 3), edgelist_with_lasttoggle(nw0))
    nw0 %n% "lasttoggle" <- lt
    expect_equal(matrix(0L, 0, 3), edgelist_with_lasttoggle(nw0))

    nw %n% "lasttoggle" <- NULL
    expect_equal(cbind(as.edgelist(nw), as.integer(-.Machine$integer.max/2)), edgelist_with_lasttoggle(nw))
  }
})
