#  File tests/testthat/test-discord-proposals.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################
#  File tests/testthat/test-discord-proposals.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

context("test-discord-proposals.R")

test_that("discordStratTNT behaves reasonably", {

  net_size <- 500L

  nw <- network.initialize(net_size, dir = FALSE)

  vattr <- sample(c("A","B","C"), net_size, TRUE)
  
  nw %v% "vattr" <- vattr
  
  pmat <- 1 - matrix(c(1,0,0,0,1,0,0,0,0),3,3)
    
  nw_sim <- nw
  
  for(i in 1:5) {
    nw_sim <- simulate(nw_sim ~ edges, 
                       coef = c(-3), 
                       time.slices = 5,
                       constraints = ~Strat(attr = "vattr", pmat = pmat),
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
  for(deg_bound in 1:5) {
    net_size <- 500L
  
    nw <- network.initialize(net_size, dir = FALSE)
  
    vattr <- sample(c("A","B","C"), net_size, TRUE)
    
    nw %v% "vattr" <- vattr
    
    fmat <- matrix(c(1,0,0,0,1,0,0,0,0),3,3)
        
    nw_sim <- nw
    
    for(i in 1:5) {
      nw_sim <- simulate(nw_sim ~ edges, 
                         coef = c(0),
                         time.slices = 5,                       
                         constraints = ~bd(maxout = deg_bound) + blocks(attr = "vattr", fmat = fmat),
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
  for(deg_bound in 1:5) {
    net_size <- 2000L
  
    nw <- network.initialize(net_size, dir = FALSE)
  
    vattr <- sample(c("A","B","C"), net_size, TRUE)
    sex <- sample(c("X","Y","Z"), net_size, TRUE)
    
    nw %v% "vattr" <- vattr
    nw %v% "sex" <-  sex
    
    fmat <- matrix(c(1,0,1,0,0,0,1,0,0),3,3)
    pmat <- 1 - matrix(c(1,0,0,0,1,0,0,0,0),3,3)
        
    nw_sim <- nw
    
    for(i in 1:5) {
      nw_sim <- simulate(nw_sim ~ edges, 
                         coef = c(0), 
                         time.slices = 5,
                         constraints = ~bd(maxout = deg_bound) + blocks(attr = "sex", fmat = fmat) + Strat(attr = "vattr", pmat = pmat),
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
  
  fmat <- matrix(0,length(sex_levels),length(sex_levels))
  colnames(fmat) <- rownames(fmat) <- sex_levels
  
  pmat <- matrix(0,length(uncommon_levels),length(uncommon_levels))
  colnames(pmat) <- rownames(pmat) <- uncommon_levels
  
  fmat["F","F"] <- 1
  fmat["M","M"] <- 1
  fmat[,"O"] <- 1
  fmat["O",] <- 1
  
  
  pmat["A","B"] <- 1
  pmat["B","A"] <- 1
  pmat["A","C"] <- 1
  pmat["C","A"] <- 1
  
  control <- control.simulate.network.tergm(MCMC.prop.weights = "discordBDStratTNT", 
                                            MCMC.prop.args = list(BD_attr = "sex",
                                                                  fmat = fmat,
                                                                  Strat_attr = "uncommon",
                                                                  pmat = pmat))
  
  out_one <- simulate(nw ~ edges, 
                      coef = c(2000), 
                      constraints = ~blocks(attr = "sex", fmat = fmat) + Strat(attr = "uncommon", pmat = pmat),
                      dynamic = TRUE,
                      output = "final")
  
  expect_true(summary(out_one ~ nodemix("sex", levels2=TRUE))["mix.sex.F.M"] == 100L)
  expect_true(summary(out_one ~ nodemix("uncommon", levels2=TRUE))["mix.uncommon.A.B"] == 100L)
  expect_true(network.edgecount(out_one) == 100L)
  
  out_two <- simulate(nw ~ edges, 
                      coef = c(2000), 
                      constraints = ~bd(maxout = 2) + blocks(attr = "sex", fmat = fmat) + Strat(attr = "uncommon", pmat = pmat),
                      dynamic = TRUE,
                      output = "final")
  
  expect_true(summary(out_two ~ nodemix("sex", levels2=TRUE))["mix.sex.F.M"] == network.edgecount(out_two))
  expect_true(summary(out_two ~ nodemix("uncommon", levels2=TRUE))["mix.uncommon.A.B"] == network.edgecount(out_two))
  expect_true(19 <= network.edgecount(out_two) && network.edgecount(out_two) <= 20L)
  
  ## now allow ties to O but continue to give them zero proposal weight
  fmat["O",] <- 0
  fmat[,"O"] <- 0

  out_three <- simulate(nw ~ edges, 
                        coef = c(2000), 
                        constraints = ~bd(maxout = 2) + blocks(attr = "sex", fmat = fmat) + Strat(attr = "uncommon", pmat = pmat),
                        dynamic = TRUE,
                        output = "final")
  
  expect_true(summary(out_three ~ nodemix("sex", levels2=TRUE))["mix.sex.F.M"] == network.edgecount(out_three))
  expect_true(summary(out_three ~ nodemix("uncommon", levels2=TRUE))["mix.uncommon.A.B"] == network.edgecount(out_three))
  expect_true(19 <= network.edgecount(out_three) && network.edgecount(out_three) <= 20L)

  ## now allow M-M ties and remove degree bound
  
  fmat["M","M"] <- 0
  
  out_four <- simulate(nw ~ edges, 
                       coef = c(2000), 
                       constraints = ~blocks(attr = "sex", fmat = fmat) + Strat(attr = "uncommon", pmat = pmat),
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
                       constraints = ~bd(maxout = 1) + blocks(attr = "sex", fmat = fmat) + Strat(attr = "uncommon", pmat = pmat),
                       dynamic = TRUE,
                       output = "final")

## A-B (M-F) and C-C (M-M) are only ties we should see
  expect_true(sum(summary(out_five ~ nodemix("uncommon", levels2=TRUE))[c("mix.uncommon.A.B", "mix.uncommon.C.C")]) == network.edgecount(out_five))  
  expect_true(network.edgecount(out_five) == 15L)  


  out_six <- simulate(out_five ~ edges, 
                      coef = c(0), 
                      constraints = ~bd(maxout = 1) + blocks(attr = "sex", fmat = fmat) + Strat(attr = "uncommon", pmat = pmat),
                      dynamic = TRUE,
                      output = "final")

  expect_true(sum(summary(out_six ~ nodemix("uncommon", levels2=TRUE))[c("mix.uncommon.A.B", "mix.uncommon.C.C")]) == network.edgecount(out_six))  
  
  out_seven <- simulate(out_six ~ edges, 
                        coef = c(0), 
                        constraints = ~bd(maxout = 1) + blocks(attr = "sex", fmat = fmat) + Strat(attr = "uncommon", pmat = pmat),
                        dynamic = TRUE,
                        output = "final")

  expect_true(sum(summary(out_seven ~ nodemix("uncommon", levels2=TRUE))[c("mix.uncommon.A.B", "mix.uncommon.C.C")]) == network.edgecount(out_seven))  

  out_eight <- simulate(out_seven ~ edges, 
                        coef = c(12), 
                        constraints = ~bd(maxout = 1) + blocks(attr = "sex", fmat = fmat) + Strat(attr = "uncommon", pmat = pmat),
                        dynamic = TRUE,
                        output = "final")

  expect_true(sum(summary(out_eight ~ nodemix("uncommon", levels2=TRUE))[c("mix.uncommon.A.B", "mix.uncommon.C.C")]) == network.edgecount(out_eight))  

  out_nine <- simulate(out_eight ~ edges, 
                       coef = c(0.5), 
                       constraints = ~bd(maxout = 1) + blocks(attr = "sex", fmat = fmat) + Strat(attr = "uncommon", pmat = pmat),
                       dynamic = TRUE,
                       output = "final")

  expect_true(sum(summary(out_nine ~ nodemix("uncommon", levels2=TRUE))[c("mix.uncommon.A.B", "mix.uncommon.C.C")]) == network.edgecount(out_nine))  

  out_ten <- simulate(out_nine ~ edges, 
                      coef = c(2000), 
                      constraints = ~bd(maxout = 1) + blocks(attr = "sex", fmat = fmat) + Strat(attr = "uncommon", pmat = pmat),
                      dynamic = TRUE,
                      output = "final")

  expect_true(sum(summary(out_ten ~ nodemix("uncommon", levels2=TRUE))[c("mix.uncommon.A.B", "mix.uncommon.C.C")]) == network.edgecount(out_ten))  
  expect_true(network.edgecount(out_ten) == 15L)  

  out_eleven <- simulate(out_ten ~ edges, 
                         coef = c(-2000), 
                         constraints = ~bd(maxout = 1) + blocks(attr = "sex", fmat = fmat) + Strat(attr = "uncommon", pmat = pmat),
                         dynamic = TRUE,
                         output = "final")

  expect_true(network.edgecount(out_eleven) == 0L)  

})

