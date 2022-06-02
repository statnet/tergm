#  File tests/testthat/test-discord-proposals.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2022 Statnet Commons
################################################################################

NITER <- 2

test_that("discordStratTNT behaves reasonably", {

  net_size <- 500L

  nw <- network.initialize(net_size, dir = FALSE)

  vattr <- sample(c("A","B","C"), net_size, TRUE)
  
  nw %v% "vattr" <- vattr
  
  pmat <- 1 - matrix(c(1,0,0,0,1,0,0,0,0),3,3)
    
  nw_sim <- nw
  
  for(i in 1:NITER) {
    nw_sim <- simulate(nw_sim ~ edges, 
                       coef = c(-5), 
                       time.slices = 3,
                       constraints = ~strat(attr = "vattr", pmat = pmat),
                       dynamic = TRUE,
                       output = "final")
    summ_stats <- summary(nw_sim ~ nodemix("vattr", levels2=TRUE))
    expect_true(summ_stats["mix.vattr.A.A"] == 0)
    expect_true(summ_stats["mix.vattr.B.B"] == 0)
    expect_true(summ_stats["mix.vattr.A.B"] > 0)
    expect_true(summ_stats["mix.vattr.A.C"] > 0)
    expect_true(summ_stats["mix.vattr.B.C"] > 0)
    expect_true(summ_stats["mix.vattr.C.C"] > 0)    
  }  
})


test_that("discordBDTNT behaves reasonably", {
  for(deg_bound in c(2)) {
    net_size <- 500L
  
    nw <- network.initialize(net_size, dir = FALSE)
  
    vattr <- sample(c("A","B","C"), net_size, TRUE)
    
    nw %v% "vattr" <- vattr
    
    levels2 <- matrix(c(1,0,0,0,1,0,0,0,0),3,3)
    levels2 <- levels2 > 0
    
    nw_sim <- nw
    
    for(i in 1:NITER) {
      nw_sim <- simulate(nw_sim ~ edges, 
                         coef = c(-4),
                         time.slices = 3,                       
                         constraints = ~bd(maxout = deg_bound) + blocks(attr = "vattr", levels2 = levels2),
                         dynamic = TRUE,
                         output = "final")
      summ_stats <- summary(nw_sim ~ nodemix("vattr", levels2=TRUE) + degrange(deg_bound + 1))
      expect_true(summ_stats[paste0("deg", deg_bound + 1, "+")] == 0)
      expect_true(summ_stats["mix.vattr.A.A"] == 0)
      expect_true(summ_stats["mix.vattr.B.B"] == 0)
      expect_true(summ_stats["mix.vattr.A.B"] > 0)
      expect_true(summ_stats["mix.vattr.A.C"] > 0)
      expect_true(summ_stats["mix.vattr.B.C"] > 0)
      expect_true(summ_stats["mix.vattr.C.C"] > 0)    
    }
  }  
})


test_that("discordBDStratTNT behaves reasonably", {
  for(deg_bound in c(1,3)) {
    net_size <- 2000L
  
    nw <- network.initialize(net_size, dir = FALSE)
  
    vattr <- sample(c("A","B","C"), net_size, TRUE)
    sex <- sample(c("X","Y","Z"), net_size, TRUE)
    
    nw %v% "vattr" <- vattr
    nw %v% "sex" <-  sex
    
    levels2 <- matrix(c(1,0,1,0,0,0,1,0,0),3,3)
    levels2 <- levels2 > 0
    
    pmat <- 1 - matrix(c(1,0,0,0,1,0,0,0,0),3,3)
        
    nw_sim <- nw
    
    for(i in 1:NITER) {
      nw_sim <- simulate(nw_sim ~ edges, 
                         coef = c(-5), 
                         time.slices = 3,
                         constraints = ~bd(maxout = deg_bound) + blocks(attr = "sex", levels2 = levels2) + strat(attr = "vattr", pmat = pmat),
                         dynamic = TRUE,
                         output = "final")
      summ_stats <- summary(nw_sim ~ nodemix("vattr", levels2=TRUE) + nodemix("sex", levels2=TRUE) + degrange(deg_bound + 1))
      expect_true(summ_stats["mix.vattr.A.A"] == 0)
      expect_true(summ_stats["mix.vattr.B.B"] == 0)
      expect_true(summ_stats["mix.vattr.A.B"] > 0)
      expect_true(summ_stats["mix.vattr.A.C"] > 0)
      expect_true(summ_stats["mix.vattr.B.C"] > 0)
      expect_true(summ_stats["mix.vattr.C.C"] > 0)    
  
      expect_true(summ_stats["mix.sex.X.X"] == 0)
      expect_true(summ_stats["mix.sex.X.Z"] == 0)
      expect_true(summ_stats["mix.sex.X.Y"] > 0)
      expect_true(summ_stats["mix.sex.Y.Y"] > 0)
      expect_true(summ_stats["mix.sex.Y.Z"] > 0)
      expect_true(summ_stats["mix.sex.Z.Z"] > 0)
  
      expect_true(summ_stats[paste0("deg", deg_bound + 1, "+")] == 0)    
    }  
  }




  ### older set of tests:
  
  net_size <- 500L

  nw <- network.initialize(net_size, dir = FALSE)

  sex <- c(rep(c("M","F"), 240), rep("O", 20))
  uncommon <- c(rep(c("A","B"), 10), rep(c("C","D"), 10), rep("E", 460))
  
  nw %v% "sex" <- sex
  nw %v% "uncommon" <- uncommon
  
  sex_levels <- sort(unique(sex))
  uncommon_levels <- sort(unique(uncommon))
  
  levels2 <- matrix(FALSE,length(sex_levels),length(sex_levels))
  colnames(levels2) <- rownames(levels2) <- sex_levels
  
  pmat <- matrix(0,length(uncommon_levels),length(uncommon_levels))
  colnames(pmat) <- rownames(pmat) <- uncommon_levels
  
  levels2["F","F"] <- TRUE
  levels2["M","M"] <- TRUE
  levels2[,"O"] <- TRUE
  levels2["O",] <- TRUE
  
  
  pmat["A","B"] <- 1
  pmat["B","A"] <- 1
  pmat["A","C"] <- 1
  pmat["C","A"] <- 1
  
  control <- control.simulate.formula.tergm(MCMC.prop.weights = "discordBDStratTNT", 
                                            MCMC.prop.args = list(BD_attr = "sex",
                                                                  levels2 = levels2,
                                                                  Strat_attr = "uncommon",
                                                                  pmat = pmat))
  
  out_one <- simulate(nw ~ edges, 
                      coef = c(2000), 
                      constraints = ~blocks(attr = "sex", levels2 = levels2) + strat(attr = "uncommon", pmat = pmat),
                      dynamic = TRUE,
                      output = "final")
  
  expect_true(summary(out_one ~ nodemix("sex", levels2=TRUE))["mix.sex.F.M"] == 100L)
  expect_true(summary(out_one ~ nodemix("uncommon", levels2=TRUE))["mix.uncommon.A.B"] == 100L)
  expect_true(network.edgecount(out_one) == 100L)
  
  out_two <- simulate(nw ~ edges, 
                      coef = c(2000), 
                      constraints = ~bd(maxout = 2) + blocks(attr = "sex", levels2 = levels2) + strat(attr = "uncommon", pmat = pmat),
                      dynamic = TRUE,
                      output = "final")
  
  expect_true(summary(out_two ~ nodemix("sex", levels2=TRUE))["mix.sex.F.M"] == network.edgecount(out_two))
  expect_true(summary(out_two ~ nodemix("uncommon", levels2=TRUE))["mix.uncommon.A.B"] == network.edgecount(out_two))
  expect_true(19 <= network.edgecount(out_two) && network.edgecount(out_two) <= 20L)
  
  ## now allow ties to O but continue to give them zero proposal weight
  levels2["O",] <- FALSE
  levels2[,"O"] <- FALSE

  out_three <- simulate(nw ~ edges, 
                        coef = c(2000), 
                        constraints = ~bd(maxout = 2) + blocks(attr = "sex", levels2 = levels2) + strat(attr = "uncommon", pmat = pmat),
                        dynamic = TRUE,
                        output = "final")
  
  expect_true(summary(out_three ~ nodemix("sex", levels2=TRUE))["mix.sex.F.M"] == network.edgecount(out_three))
  expect_true(summary(out_three ~ nodemix("uncommon", levels2=TRUE))["mix.uncommon.A.B"] == network.edgecount(out_three))
  expect_true(19 <= network.edgecount(out_three) && network.edgecount(out_three) <= 20L)

  ## now allow M-M ties and remove degree bound
  
  levels2["M","M"] <- FALSE
  
  out_four <- simulate(nw ~ edges, 
                       coef = c(2000), 
                       constraints = ~blocks(attr = "sex", levels2 = levels2) + strat(attr = "uncommon", pmat = pmat),
                       dynamic = TRUE,
                       output = "final")
  
  expect_true(summary(out_four ~ nodemix("sex", levels2=TRUE))["mix.sex.F.M"] == 100L)
  expect_true(summary(out_four ~ nodemix("uncommon", levels2=TRUE))["mix.uncommon.A.B"] == 100L)
  expect_true(summary(out_four ~ nodemix("sex", levels2=TRUE))["mix.sex.M.M"] == 100L)
  expect_true(summary(out_four ~ nodemix("uncommon", levels2=TRUE))["mix.uncommon.A.C"] == 100L)
  expect_true(network.edgecount(out_four) == 200L)  

  ## now to test both initialization code and edge-removal code, start at the empty network and go up and down in 
  ## edge count repeatedly (with the same weights and constraints each time), checking that we satisfy the necessary
  ## conditions at each stage
  
  pmat["A","C"] <- 0
  pmat["C","A"] <- 0
  pmat["C","C"] <- 1
  pmat["B","B"] <- 1
  pmat["D","D"] <- 1
  pmat["B","D"] <- 1
  pmat["D","B"] <- 1

  out_five <- simulate(nw ~ edges, 
                       coef = c(2000), 
                       constraints = ~bd(maxout = 1) + blocks(attr = "sex", levels2 = levels2) + strat(attr = "uncommon", pmat = pmat),
                       dynamic = TRUE,
                       output = "final")

## A-B (M-F) and C-C (M-M) are only ties we should see
  expect_true(sum(summary(out_five ~ nodemix("uncommon", levels2=TRUE))[c("mix.uncommon.A.B", "mix.uncommon.C.C")]) == network.edgecount(out_five))  
  expect_true(network.edgecount(out_five) == 15L)  


  out_six <- simulate(out_five ~ edges, 
                      coef = c(0), 
                      constraints = ~bd(maxout = 1) + blocks(attr = "sex", levels2 = levels2) + strat(attr = "uncommon", pmat = pmat),
                      dynamic = TRUE,
                      output = "final")

  expect_true(sum(summary(out_six ~ nodemix("uncommon", levels2=TRUE))[c("mix.uncommon.A.B", "mix.uncommon.C.C")]) == network.edgecount(out_six))  
  
  out_seven <- simulate(out_six ~ edges, 
                        coef = c(0), 
                        constraints = ~bd(maxout = 1) + blocks(attr = "sex", levels2 = levels2) + strat(attr = "uncommon", pmat = pmat),
                        dynamic = TRUE,
                        output = "final")

  expect_true(sum(summary(out_seven ~ nodemix("uncommon", levels2=TRUE))[c("mix.uncommon.A.B", "mix.uncommon.C.C")]) == network.edgecount(out_seven))  

  out_eight <- simulate(out_seven ~ edges, 
                        coef = c(12), 
                        constraints = ~bd(maxout = 1) + blocks(attr = "sex", levels2 = levels2) + strat(attr = "uncommon", pmat = pmat),
                        dynamic = TRUE,
                        output = "final")

  expect_true(sum(summary(out_eight ~ nodemix("uncommon", levels2=TRUE))[c("mix.uncommon.A.B", "mix.uncommon.C.C")]) == network.edgecount(out_eight))  

  out_nine <- simulate(out_eight ~ edges, 
                       coef = c(0.5), 
                       constraints = ~bd(maxout = 1) + blocks(attr = "sex", levels2 = levels2) + strat(attr = "uncommon", pmat = pmat),
                       dynamic = TRUE,
                       output = "final")

  expect_true(sum(summary(out_nine ~ nodemix("uncommon", levels2=TRUE))[c("mix.uncommon.A.B", "mix.uncommon.C.C")]) == network.edgecount(out_nine))  

  out_ten <- simulate(out_nine ~ edges, 
                      coef = c(2000), 
                      constraints = ~bd(maxout = 1) + blocks(attr = "sex", levels2 = levels2) + strat(attr = "uncommon", pmat = pmat),
                      dynamic = TRUE,
                      output = "final")

  expect_true(sum(summary(out_ten ~ nodemix("uncommon", levels2=TRUE))[c("mix.uncommon.A.B", "mix.uncommon.C.C")]) == network.edgecount(out_ten))  
  expect_true(network.edgecount(out_ten) == 15L)  

  out_eleven <- simulate(out_ten ~ edges, 
                         coef = c(-2000), 
                         constraints = ~bd(maxout = 1) + blocks(attr = "sex", levels2 = levels2) + strat(attr = "uncommon", pmat = pmat),
                         dynamic = TRUE,
                         output = "final")

  expect_true(network.edgecount(out_eleven) == 0L)  

})

test_that("discordBDStratTNT simulates reasonably with heterogeneous degree bounds", {
  for(deg_bound in c(1)) {
    net_size <- 2000L
  
    nw <- network.initialize(net_size, dir = FALSE)
  
    vattr <- sample(c("A","B","C"), net_size, TRUE)
    sex <- sample(c(1,2,3), net_size, TRUE)
    
    attribs <- matrix(FALSE, nrow = net_size, ncol = 3)
    attribs[cbind(seq_len(net_size), sex)] <- TRUE    
    
    nw %v% "vattr" <- vattr
    nw %v% "sex" <-  sex
    nw %v% "blocks_attr" <- sample(1:6, net_size, TRUE)
    
    blocks_levels_2 <- matrix(FALSE, 6, 6)
    blocks_levels_2[cbind(c(1,2,2,4), c(5,2,3,4))] <- TRUE
    blocks_levels_2 <- blocks_levels_2 | t(blocks_levels_2)
    
    levels2 <- matrix(c(1,0,1,0,0,0,1,0,0),3,3)
    levels2 <- levels2 > 0
    
    pmat <- 1 - matrix(c(1,0,0,0,1,0,0,0,0),3,3)
          
    nw_sim <- nw

    maxout <- matrix(0, nrow = net_size, ncol = 3)

    for(row_index in 1:3) {
      for(col_index in 1:3) {
        if(!levels2[row_index, col_index]) {
          maxout[sex == row_index, col_index] <- deg_bound
        }
      }
    }
    maxout <- maxout + round(5*(runif(length(maxout)) - 1/2))
    maxout[maxout < 0] <- 0
    
    for(i in 1:NITER) {
      nw_sim <- simulate(nw_sim ~ Form(~edges) + Persist(~edges),
                         coef = c(-7,2),
                         dynamic = TRUE,
                         time.slices = 3,                         
                         constraints = ~bd(attribs = attribs, maxout = maxout) + blocks(~blocks_attr, levels2 = blocks_levels_2) + strat(attr = "vattr", pmat = pmat),
                         output = "final")
      
      summ_stats_vattr <- summary(nw_sim ~ nodemix("vattr",levels2=TRUE))
      expect_true(all(summ_stats_vattr[c(1,3)] == 0))
      expect_true(all(summ_stats_vattr[-c(1,3)] > 0))
  
      summ_stats_blocks_attr <- summary(nw_sim ~ nodemix("blocks_attr",levels2=TRUE))
      expect_true(all(summ_stats_blocks_attr[c(3,5,10,11)] == 0))
      expect_true(all(summ_stats_blocks_attr[-c(3,5,10,11)] > 0))
      
      el <- as.edgelist(nw_sim)
      degs <- table(from = factor(c(el), levels = seq_len(net_size)), to = factor(sex[c(el[,c(2,1)])], levels = seq_len(3)))
      expect_true(all(degs <= maxout))
    }  
  }
})

test_that("discordBDStratTNT simulates reasonably with bipartite heterogeneous degree bounds", {
  for(deg_bound in c(1)) {
    net_size <- 2000L
    bip <- 700L
    
    nw <- network.initialize(net_size, dir = FALSE, bip = bip)
  
    vattr <- c(sample(c("A","B","C","D"), bip, TRUE), sample(c("X","Y","Z"), net_size - bip, TRUE))
    sex <- c(sample(c(1,2,3,4,5), bip, TRUE), sample(c(6,7,8,9,10,11), net_size - bip, TRUE))
    
    attribs <- matrix(FALSE, nrow = net_size, ncol = length(unique(sex)))
    attribs[cbind(seq_len(net_size), sex)] <- TRUE    
    
    nw %v% "vattr" <- vattr
    nw %v% "sex" <-  sex
    nw %v% "blocks_attr" <- c(sample(c(1,2,3), bip, TRUE), sample(c(4,5,6,7), net_size - bip, TRUE))
    
    blocks_levels_2 <- matrix(FALSE, nrow = 3, 4)
    blocks_levels_2[cbind(c(3,2,2), c(1,2,3))] <- TRUE
    
    levels2 <- matrix(as.logical(round(runif(11*11))), nrow = 11, ncol = 11)
    levels2 <- levels2 | t(levels2)
    
    pmat <- 1 - matrix(c(1,0,0,0,1,0,1,0,0,0,0,1),nrow = 4, ncol = 3)
          
    nw_sim <- nw

    maxout <- matrix(0, nrow = net_size, ncol = 11)

    for(row_index in 1:11) {
      for(col_index in 1:11) {
        if(!levels2[row_index, col_index]) {
          maxout[sex == row_index, col_index] <- deg_bound
        }
      }
    }
    maxout <- maxout + round(5*(runif(length(maxout)) - 1/2))
    maxout[maxout < 0] <- 0
    
    for(i in 1:NITER) {
      nw_sim <- simulate(nw_sim ~ Form(~edges) + Persist(~edges),
                         coef = c(-6,2),
                         dynamic = TRUE,
                         time.slices = 3,                         
                         constraints = ~bd(attribs = attribs, maxout = maxout) + blocks(~blocks_attr, levels2 = blocks_levels_2) + strat(attr = "vattr", pmat = pmat),
                         output = "final")
      
      summ_stats_vattr <- summary(nw_sim ~ nodemix("vattr",levels2=TRUE))
      expect_true(all(summ_stats_vattr[c(1,5,7,12)] == 0))
      expect_true(all(summ_stats_vattr[-c(1,5,7,12)] > 0))
  
      summ_stats_blocks_attr <- summary(nw_sim ~ nodemix("blocks_attr",levels2=TRUE))
      expect_true(all(summ_stats_blocks_attr[c(3,5,8)] == 0))
      expect_true(all(summ_stats_blocks_attr[-c(3,5,8)] > 0))
      
      el <- as.edgelist(nw_sim)
      degs <- table(from = factor(c(el), levels = seq_len(net_size)), to = factor(sex[c(el[,c(2,1)])], levels = seq_len(11)))
      expect_true(all(degs <= maxout))
    }  
  }
})

test_that("discordBDStratTNT simulates reasonably with directed heterogeneous degree bounds", {
  for(deg_bound in c(1)) {
    net_size <- 2000L
  
    nw <- network.initialize(net_size, dir = TRUE)
  
    vattr <- sample(c("A","B","C"), net_size, TRUE)
    sex <- sample(c(1,2,3), net_size, TRUE)
    
    attribs <- matrix(FALSE, nrow = net_size, ncol = 3)
    attribs[cbind(seq_len(net_size), sex)] <- TRUE    
    
    nw %v% "vattr" <- vattr
    nw %v% "sex" <-  sex
    nw %v% "blocks_attr" <- sample(1:6, net_size, TRUE)
    
    blocks_levels_2 <- matrix(FALSE, 6, 6)
    blocks_levels_2[cbind(c(5,2,2,4), c(1,2,3,4))] <- TRUE
    
    levels2 <- matrix(c(1,0,0,0,0,1,1,0,0),3,3)
    levels2 <- levels2 > 0
    
    pmat <- 1 - matrix(c(1,0,0,0,0,0,0,1,0),3,3)
          
    nw_sim <- nw
    
    maxout <- matrix(0, nrow = net_size, ncol = 3)

    for(row_index in 1:3) {
      for(col_index in 1:3) {
        if(!levels2[row_index, col_index]) {
          maxout[sex == row_index, col_index] <- deg_bound
        }
      }
    }
    maxout <- maxout + round(5*(runif(length(maxout)) - 1/2))
    maxout[maxout < 0] <- 0
    
    maxin <- maxout + round(5*(runif(length(maxout)) - 1/2))
    maxin[maxin < 0] <- 0
    
    for(i in 1:NITER) {
      nw_sim <- simulate(nw_sim ~ Form(~edges) + Persist(~edges),
                         coef = c(-7,2),
                         dynamic = TRUE,
                         time.slices = 3,                         
                         constraints = ~bd(attribs = attribs, maxout = maxout, maxin = maxin) + blocks(~blocks_attr, levels2 = blocks_levels_2) + strat(attr = "vattr", pmat = pmat),
                         output = "final")
      
      summ_stats_vattr <- summary(nw_sim ~ nodemix("vattr",levels2=TRUE))
      expect_true(all(summ_stats_vattr[c(1,8)] == 0))
      expect_true(all(summ_stats_vattr[-c(1,8)] > 0))
  
      summ_stats_blocks_attr <- summary(nw_sim ~ nodemix("blocks_attr",levels2=TRUE))
      expect_true(all(summ_stats_blocks_attr[c(5,8,14,22)] == 0))
      expect_true(all(summ_stats_blocks_attr[-c(5,8,14,22)] > 0))
      
      el <- as.edgelist(nw_sim)
      out_degs <- table(from = factor(c(el[,1]), levels = seq_len(net_size)), to = factor(sex[c(el[,2])], levels = seq_len(3)))
      in_degs <- table(from = factor(c(el[,2]), levels = seq_len(net_size)), to = factor(sex[c(el[,1])], levels = seq_len(3)))
      expect_true(all(out_degs <= maxout))
      expect_true(all(in_degs <= maxin))      
    }  
  }
})

