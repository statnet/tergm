#  File tests/testthat/test-durational-terms.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################

test_that("durational terms behave correctly with summary and godfather", {
  # nwtype = 1 is undirected unipartite
  # nwtype = 2 is undirected bipartite
  # nwtype = 3 is directed unipartite
  for(nwtype in 1:3) {
  
    netsize <- 100
    bipsize <- 30
    
    init_tmin <- 0
    init_tmax <- 10
    init_time <- init_tmax
  
    skip_frac <- 0.3
  
    sampling_frac <- 0.8
  
    toggle_tmin <- 11
    toggle_tmax <- 30
  
    mean_age_levels <- 7
    mean_age_b1levels <- 6
    mean_age_b2levels <- 8
  
    max_levels <- 5
  
    if(nwtype == 1) {
      nw <- network.initialize(netsize,dir=FALSE)
    } else if(nwtype == 2) {
      nw <- network.initialize(netsize,dir=FALSE,bip=bipsize)
    } else {
      nw <- network.initialize(netsize,dir=TRUE)  
    }
    
    nw <- san(nw ~ edges, target.stats = c(3*netsize))
    
    el <- as.edgelist(nw)
  
    attrname <- "race"
  
    attr <- sample(letters[1:max_levels], netsize, TRUE)
    levels <- sort(unique(attr))
  
    nw %v% attrname <- attr
    
    if(nwtype == 1) {
      mean_age_attr <- rep(seq_len(mean_age_levels), length.out=netsize)
      n_nodefactor_stats <- mean_age_levels
      n_nodemix_stats <- mean_age_levels*(mean_age_levels + 1)/2
      indmat <- matrix(0, nrow=mean_age_levels, ncol=mean_age_levels)
      indmat[upper.tri(indmat,diag=TRUE)] <- seq_len(n_nodemix_stats)
      indmat <- indmat + t(indmat) - diag(diag(indmat))
      nltk <- ceiling(mean_age_levels/2)
      nodemix_ma_lev <- sample(seq_len(mean_age_levels), nltk, FALSE)
      nstats_nodemix_no2 <- length(nodemix_ma_lev)*(length(nodemix_ma_lev) + 1)/2
      nodemix_ma_lev2 <- sample(seq_len(nltk*(nltk + 1)/2), ceiling(nltk*(nltk + 1)/4), FALSE)
      indmat_subset <- indmat[nodemix_ma_lev, nodemix_ma_lev]
      nodemix_indices_1 <- unique(indmat_subset[upper.tri(indmat_subset, diag=TRUE)])
    } else if(nwtype == 2) {
      mean_age_attr <- c(rep(seq_len(mean_age_b1levels), length.out=bipsize), rep(seq_len(mean_age_b2levels), length.out=netsize-bipsize))
      n_nodefactor_stats <- max(mean_age_b1levels, mean_age_b2levels)
      n_nodemix_stats <- mean_age_b1levels*mean_age_b2levels
      indmat <- matrix(seq_len(n_nodemix_stats),nrow=mean_age_b1levels,ncol=mean_age_b2levels)
      nltkb1 <- ceiling(mean_age_b1levels/2)
      nltkb2 <- ceiling(mean_age_b2levels/2)
      nodemix_ma_b1lev <- sample(seq_len(mean_age_b1levels), nltkb1, FALSE)
      nodemix_ma_b2lev <- sample(seq_len(mean_age_b2levels), nltkb2, FALSE)
      nstats_nodemix_no2 <- length(nodemix_ma_b1lev)*length(nodemix_ma_b2lev)
      nodemix_ma_lev2 <- sample(seq_len(nltkb1*nltkb2), ceiling(nltkb1*nltkb2/2), FALSE)
      nodemix_indices_1 <- c(indmat[nodemix_ma_b1lev, nodemix_ma_b2lev])
    } else {
      mean_age_attr <- rep(seq_len(mean_age_levels), length.out=netsize)
      n_nodefactor_stats <- mean_age_levels
      n_nodemix_stats <- mean_age_levels*mean_age_levels
      indmat <- matrix(seq_len(n_nodemix_stats),nrow=mean_age_levels,ncol=mean_age_levels)
      nltk <- ceiling(mean_age_levels/2)
      nodemix_ma_lev <- sample(seq_len(mean_age_levels), nltk, FALSE)
      nstats_nodemix_no2 <- length(nodemix_ma_lev)*length(nodemix_ma_lev)
      nodemix_ma_lev2 <- sample(seq_len(nltk*nltk/2), ceiling(nltk*nltk/4), FALSE)
      nodemix_indices_1 <- c(indmat[nodemix_ma_lev, nodemix_ma_lev])
    }

    nodemix_indices_2 <- nodemix_indices_1[nodemix_ma_lev2]
    nodefactor_ma_lev <- sample(seq_len(n_nodefactor_stats), ceiling(n_nodefactor_stats/2), FALSE)
    
    nodemix_emptyvals <- runif(length(nodemix_indices_1))
    nodefactor_emptyvals <- runif(length(nodefactor_ma_lev))
    
    nw %v% "mean_age_attr" <- mean_age_attr
    
    init.lasttoggle.times <- sample(init_tmin:init_tmax, NROW(el), TRUE)
    
    lt <- cbind(el, init.lasttoggle.times)
      
    nw %n% "time" <- init_time  
    nw %n% "lasttoggle" <- lt
    
    if(nwtype == 1) {
      wts <- matrix(runif(netsize*netsize),netsize,netsize)
      wts <- wts + t(wts)
    } else if(nwtype == 2) {
      wts <- matrix(runif(bipsize*(netsize - bipsize)),nrow=bipsize,ncol=netsize-bipsize)  
    } else {
      wts <- matrix(runif(netsize*netsize),netsize,netsize)  
    }
    
    nw %n% "edgewts" <- wts
    
    if(nwtype == 1) {
      nw0 <- network.initialize(netsize,dir=FALSE)
    } else if(nwtype == 2) {
      nw0 <- network.initialize(netsize,dir=FALSE,bip=bipsize)    
    } else {
      nw0 <- network.initialize(netsize,dir=TRUE)  
    }

    nw0 %v% attrname <- attr
    nw0 %v% "mean_age_attr" <- mean_age_attr
    nw0 %n% "edgewts" <- replace(wts, TRUE, 0) # A 0 matrix with the same dimensions.
    
    ages_from = c(1, 1, 4, 2, 5, 10, 1, 5)
    ages_to = c(2, 3, 15, 5, 7, 12, Inf, Inf) 
    
    init_ages <- init_time - init.lasttoggle.times + 1
    
    init_age_counts <- integer(length(ages_from))
    for(i in seq_len(length(ages_from))) {
      init_age_counts[i] <- sum(init_ages >= ages_from[i] & init_ages < ages_to[i])
    }
  
    ## NB: degree.mean.age and degrange.mean.age do not currently support directed networks,
    ##     so code related to these terms only needs to work for undirected networks
    degree_vec <- 1:10
    degrange_from <- c(1:5, 5:1)
    degrange_to <- c(3,5,4,8,10,6,5,7,3,2)
  
    ## in attr case, it is the attr of the central node that matters; all neighbors count towards degree regardless
    ## of their attr type
    
    ## note degree from edgelist is easy, just sum number of occurrences of node
  
    drma_from_el_lt <- function(dvec_from, dvec_to, el_lt, time, zeroval = 0, attr=integer(netsize), log=FALSE) {
      levels <- sort(unique(attr))
      outvec <- numeric(length(dvec_from)*length(levels))
      
      degs <- integer(netsize)
      for(i in seq_len(netsize)) {
        degs[i] <- sum(el_lt[,1:2] == i)
      }
      
      for(j in seq_len(length(levels))) {
        for(i in seq_len(length(dvec_from))) {
          wod <- which(degs >= dvec_from[i] & degs < dvec_to[i])
          woa <- which(attr == levels[j])
          
          edgecounts <- ((el_lt[,1] %in% wod) & (el_lt[,1] %in% woa)) + ((el_lt[,2] %in% wod) & (el_lt[,2] %in% woa))
    
          if(any(edgecounts)) {
            transf_ages <- if(log) log(time - el_lt[,3] + 1) else time - el_lt[,3] + 1
            sumages <- sum(edgecounts*transf_ages)
            outvec[(j - 1)*length(dvec_from) + i] <- sumages/sum(edgecounts)
          } else {
            outvec[(j - 1)*length(dvec_from) + i] <- zeroval
          }
        }
      }
      outvec
    }
  
    # takes initial edgelist with lasttoggle times and the set of toggles and
    # computes edgelists with lasttoggles for all times between the first and last toggle,
    # assuming the toggles take place after the lasttoggles in lt_input
    map_toggles_to_edgelist_lts <- function(lt_input, toggles) {
      tmin <- min(toggles[,1])
      tmax <- max(toggles[,1])
      
      outlist <- list()
      
      lt_current <- lt_input
      for(timestep in tmin:tmax) {
        w <- which(toggles[,1] == timestep)
        
        if(length(w) > 0) {
          these_toggles <- toggles[w,,drop=FALSE]
          
          ## need to cancel "double toggles" so lasttoggle time doesn't change incorrectly
          keep <- rep(TRUE, length(w))
          for(i in seq_len(NROW(these_toggles) - 1)) {
            if(keep[i]) {
              for(j in (i + 1):NROW(these_toggles)) {
                if(these_toggles[i,2] == these_toggles[j,2] && these_toggles[i,3] == these_toggles[j,3] && keep[j]) {
                  keep[i] <- FALSE
                  keep[j] <- FALSE
                  break
                }
              }
            }
          }
          
          these_toggles <- these_toggles[keep,,drop=FALSE]
          
          for(i in seq_len(NROW(these_toggles))) {
            found <- FALSE
            for(j in seq_len(NROW(lt_current))) {
              if(lt_current[j,1] == these_toggles[i,2] && lt_current[j,2] == these_toggles[i,3]) {
                lt_current <- lt_current[-j,,drop=FALSE]
                found <- TRUE
                break
              }
            }
            if(!found) {
              lt_current <- rbind(lt_current, these_toggles[i, c(2,3,1)])
            }
          }
        }
        outlist[[timestep]] <- lt_current
      }
      outlist
    }
  
    elb <- el
    if(nwtype == 2) {
      elb[,2] <- elb[,2] - bipsize
    }
  
    expect_equal(unname(summary(dynamic=TRUE, nw0 ~ mean.age(emptyval=4321.4))), 4321.4)  
    expect_equal(unname(summary(nw ~ mean.age)), mean(init_time - init.lasttoggle.times + 1))
    
    expect_equal(unname(summary(dynamic=TRUE, nw0 ~ edge.ages)), 0)  
    expect_equal(unname(summary(nw ~ edge.ages)), sum(init_time - init.lasttoggle.times + 1))
    
    expect_equal(unname(summary(dynamic=TRUE, nw0 ~ edgecov.ages("edgewts"))), 0)  
    expect_equal(unname(summary(nw ~ edgecov.ages("edgewts"))), sum(wts[elb]*(init_time - init.lasttoggle.times + 1)))
    
    expect_equal(unname(summary(dynamic=TRUE, nw0 ~ edgecov.mean.age("edgewts", emptyval=pi^4/25))), pi^4/25)  
    expect_equal(unname(summary(nw ~ edgecov.mean.age("edgewts"))), sum(wts[elb]*(init_time - init.lasttoggle.times + 1))/sum(wts[elb]))
      
    expect_equal(unname(summary(dynamic=TRUE, nw0 ~ mean.age(emptyval=pi^3/83, log=TRUE))), pi^3/83)    
    expect_equal(unname(summary(nw ~ mean.age(log=TRUE))), mean(log(init_time - init.lasttoggle.times + 1)))
  
    expect_equal(unname(summary(dynamic=TRUE, nw0 ~ edgecov.mean.age("edgewts", emptyval=pi^4/25, log=TRUE))), pi^4/25)  
    expect_equal(unname(summary(nw ~ edgecov.mean.age("edgewts", log=TRUE))), sum(wts[elb]*log(init_time - init.lasttoggle.times + 1))/sum(wts[elb]))
      
    expect_equal(unname(summary(dynamic=TRUE, nw0 ~ edges.ageinterval(from=ages_from, to=ages_to))), integer(length(ages_from)))  
    expect_equal(unname(summary(nw ~ edges.ageinterval(from=ages_from, to=ages_to))), init_age_counts)  

    ## nodefactor/nodemix mean age tests
    nodefactor_ages <- numeric(n_nodefactor_stats)
    nodefactor_log_ages <- numeric(n_nodefactor_stats)
    nodefactor_counts <- integer(n_nodefactor_stats)
    
    attr_elt <- lt
    attr_elt[,1] <- mean_age_attr[attr_elt[,1]]
    attr_elt[,2] <- mean_age_attr[attr_elt[,2]]
    for(k in seq_len(NROW(attr_elt))) {
      nodefactor_ages[attr_elt[k,1]] <- nodefactor_ages[attr_elt[k,1]] + init_time - attr_elt[k,3] + 1
      nodefactor_log_ages[attr_elt[k,1]] <- nodefactor_log_ages[attr_elt[k,1]] + log(init_time - attr_elt[k,3] + 1)
      nodefactor_counts[attr_elt[k,1]] <- nodefactor_counts[attr_elt[k,1]] + 1
      nodefactor_ages[attr_elt[k,2]] <- nodefactor_ages[attr_elt[k,2]] + init_time - attr_elt[k,3] + 1
      nodefactor_counts[attr_elt[k,2]] <- nodefactor_counts[attr_elt[k,2]] + 1
      nodefactor_log_ages[attr_elt[k,2]] <- nodefactor_log_ages[attr_elt[k,2]] + log(init_time - attr_elt[k,3] + 1)
    }
    
    nodefactor_vals <- numeric(n_nodefactor_stats)
    nodefactor_log_vals <- numeric(n_nodefactor_stats)
    for(k in seq_len(n_nodefactor_stats)) {
      if(nodefactor_counts[k] == 0) {
        ## else 0 for now; emptyval is handled later
      } else {
        nodefactor_vals[k] <- nodefactor_ages[k]/nodefactor_counts[k]
        nodefactor_log_vals[k] <- nodefactor_log_ages[k]/nodefactor_counts[k]
      }
    }
    
    nodemix_ages <- numeric(n_nodemix_stats)
    nodemix_log_ages <- numeric(n_nodemix_stats)
    nodemix_counts <- integer(n_nodemix_stats)
    
    for(k in seq_len(NROW(attr_elt))) {
      this_nodemix_index <- indmat[attr_elt[k,1],attr_elt[k,2]]
      nodemix_ages[this_nodemix_index] <- nodemix_ages[this_nodemix_index] + init_time - attr_elt[k,3] + 1
      nodemix_log_ages[this_nodemix_index] <- nodemix_log_ages[this_nodemix_index] + log(init_time - attr_elt[k,3] + 1)
      nodemix_counts[this_nodemix_index] <- nodemix_counts[this_nodemix_index] + 1
    }
    
    nodemix_vals <- numeric(n_nodemix_stats)
    nodemix_log_vals <- numeric(n_nodemix_stats)
    for(k in seq_len(n_nodemix_stats)) {
      if(nodemix_counts[k] == 0) {
        ## else 0 for now since we're not passing an emptyval
      } else {
        nodemix_vals[k] <- nodemix_ages[k]/nodemix_counts[k]
        nodemix_log_vals[k] <- nodemix_log_ages[k]/nodemix_counts[k]
      }
    }

    ## set nonzero emptyvals here
    
    
    ## at least two each for nodefactor, nodemix mean age; with and without emptyval, with and without log,
    ## with and without (b1/b2)levels(2), also levels2 without levels for nodemix at least once
    
    ## will need to get extraction of correct values as below; ellt is named lt and constructed already above
    
    expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodefactor.mean.age("mean_age_attr"))), rep(0, n_nodefactor_stats))
    expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodefactor.mean.age("mean_age_attr", levels=nodefactor_ma_lev))), rep(0, length(nodefactor_ma_lev)))
    expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodefactor.mean.age("mean_age_attr", levels=nodefactor_ma_lev, emptyval=nodefactor_emptyvals))), nodefactor_emptyvals)
    expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodefactor.mean.age("mean_age_attr", levels=nodefactor_ma_lev, emptyval=nodefactor_emptyvals, log=TRUE))), nodefactor_emptyvals)

    expect_equal(unname(summary(nw ~ nodefactor.mean.age("mean_age_attr"))), nodefactor_vals)
    expect_equal(unname(summary(nw ~ nodefactor.mean.age("mean_age_attr", levels=nodefactor_ma_lev))), nodefactor_vals[nodefactor_ma_lev])

    ## substitute emptyvals for 0s
    nodefactor_vals[nodefactor_ma_lev][nodefactor_counts[nodefactor_ma_lev] == 0] <- nodefactor_emptyvals[nodefactor_counts[nodefactor_ma_lev] == 0]
    nodefactor_log_vals[nodefactor_ma_lev][nodefactor_counts[nodefactor_ma_lev] == 0] <- nodefactor_emptyvals[nodefactor_counts[nodefactor_ma_lev] == 0]
    
    expect_equal(unname(summary(nw ~ nodefactor.mean.age("mean_age_attr", levels=nodefactor_ma_lev, emptyval=nodefactor_emptyvals))), nodefactor_vals[nodefactor_ma_lev])
    expect_equal(unname(summary(nw ~ nodefactor.mean.age("mean_age_attr", levels=nodefactor_ma_lev, emptyval=nodefactor_emptyvals, log=TRUE))), nodefactor_log_vals[nodefactor_ma_lev])
 
    if(nwtype == 2) {
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodemix.mean.age("mean_age_attr"))), rep(0, n_nodemix_stats))
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodemix.mean.age("mean_age_attr", b1levels=nodemix_ma_b1lev, b2levels=nodemix_ma_b2lev))), rep(0, length(nodemix_indices_1)))
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodemix.mean.age("mean_age_attr", levels2=nodemix_ma_lev2))), rep(0, length(nodemix_ma_lev2)))
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodemix.mean.age("mean_age_attr", b1levels=nodemix_ma_b1lev, b2levels=nodemix_ma_b2lev, log=TRUE))), rep(0, length(nodemix_indices_1)))
   
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodemix.mean.age("mean_age_attr", b1levels=nodemix_ma_b1lev, b2levels=nodemix_ma_b2lev, log=TRUE, emptyval=nodemix_emptyvals))), nodemix_emptyvals)
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodemix.mean.age("mean_age_attr", b1levels=nodemix_ma_b1lev, b2levels=nodemix_ma_b2lev, levels2=nodemix_ma_lev2, log=TRUE, emptyval=nodemix_emptyvals[nodemix_ma_lev2]))), nodemix_emptyvals[nodemix_ma_lev2])
   
      expect_equal(unname(summary(nw ~ nodemix.mean.age("mean_age_attr"))), nodemix_vals)
      expect_equal(unname(summary(nw ~ nodemix.mean.age("mean_age_attr", b1levels=nodemix_ma_b1lev, b2levels=nodemix_ma_b2lev))), nodemix_vals[nodemix_indices_1])
      expect_equal(unname(summary(nw ~ nodemix.mean.age("mean_age_attr", levels2=nodemix_ma_lev2))), nodemix_vals[nodemix_ma_lev2])
      expect_equal(unname(summary(nw ~ nodemix.mean.age("mean_age_attr", b1levels=nodemix_ma_b1lev, b2levels=nodemix_ma_b2lev, log=TRUE))), nodemix_log_vals[nodemix_indices_1])
   
      ## substitute emptyvals for 0s
      nodemix_vals[nodemix_indices_1][nodemix_counts[nodemix_indices_1] == 0] <- nodemix_emptyvals[nodemix_counts[nodemix_indices_1] == 0]
      nodemix_log_vals[nodemix_indices_1][nodemix_counts[nodemix_indices_1] == 0] <- nodemix_emptyvals[nodemix_counts[nodemix_indices_1] == 0]
   
      expect_equal(unname(summary(nw ~ nodemix.mean.age("mean_age_attr", b1levels=nodemix_ma_b1lev, b2levels=nodemix_ma_b2lev, log=TRUE, emptyval=nodemix_emptyvals))), nodemix_log_vals[nodemix_indices_1])
      expect_equal(unname(summary(nw ~ nodemix.mean.age("mean_age_attr", b1levels=nodemix_ma_b1lev, b2levels=nodemix_ma_b2lev, levels2=nodemix_ma_lev2, log=TRUE, emptyval=nodemix_emptyvals[nodemix_ma_lev2]))), nodemix_log_vals[nodemix_indices_2])
    } else {    
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodemix.mean.age("mean_age_attr"))), rep(0, n_nodemix_stats))
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodemix.mean.age("mean_age_attr", levels=nodemix_ma_lev))), rep(0, length(nodemix_indices_1)))
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodemix.mean.age("mean_age_attr", levels2=nodemix_ma_lev2))), rep(0, length(nodemix_ma_lev2)))
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodemix.mean.age("mean_age_attr", levels=nodemix_ma_lev, log=TRUE))), rep(0, length(nodemix_indices_1)))
   
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodemix.mean.age("mean_age_attr", levels=nodemix_ma_lev, log=TRUE, emptyval=nodemix_emptyvals))), nodemix_emptyvals)
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ nodemix.mean.age("mean_age_attr", levels=nodemix_ma_lev, levels2=nodemix_ma_lev2, log=TRUE, emptyval=nodemix_emptyvals[nodemix_ma_lev2]))), nodemix_emptyvals[nodemix_ma_lev2])
   
      expect_equal(unname(summary(nw ~ nodemix.mean.age("mean_age_attr"))), nodemix_vals)
      expect_equal(unname(summary(nw ~ nodemix.mean.age("mean_age_attr", levels=nodemix_ma_lev))), nodemix_vals[nodemix_indices_1])
      expect_equal(unname(summary(nw ~ nodemix.mean.age("mean_age_attr", levels2=nodemix_ma_lev2))), nodemix_vals[nodemix_ma_lev2])
      expect_equal(unname(summary(nw ~ nodemix.mean.age("mean_age_attr", levels=nodemix_ma_lev, log=TRUE))), nodemix_log_vals[nodemix_indices_1])
   
      ## substitute emptyvals for 0s
      nodemix_vals[nodemix_indices_1][nodemix_counts[nodemix_indices_1] == 0] <- nodemix_emptyvals[nodemix_counts[nodemix_indices_1] == 0]
      nodemix_log_vals[nodemix_indices_1][nodemix_counts[nodemix_indices_1] == 0] <- nodemix_emptyvals[nodemix_counts[nodemix_indices_1] == 0]
   
      expect_equal(unname(summary(nw ~ nodemix.mean.age("mean_age_attr", levels=nodemix_ma_lev, log=TRUE, emptyval=nodemix_emptyvals))), nodemix_log_vals[nodemix_indices_1])
      expect_equal(unname(summary(nw ~ nodemix.mean.age("mean_age_attr", levels=nodemix_ma_lev, levels2=nodemix_ma_lev2, log=TRUE, emptyval=nodemix_emptyvals[nodemix_ma_lev2]))), nodemix_log_vals[nodemix_indices_2])
    }
  
    if(nwtype < 3) {
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ degree.mean.age(degree_vec, emptyval=0345.4))), rep(0345.4, length(degree_vec)))  
      expect_equal(unname(summary(nw ~ degree.mean.age(degree_vec))), drma_from_el_lt(degree_vec, degree_vec + 1, lt, init_time, 0))
    
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ degrange.mean.age(degrange_from, degrange_to, emptyval=1345.4))), rep(1345.4, length(degrange_from)))  
      expect_equal(unname(summary(nw ~ degrange.mean.age(degrange_from, degrange_to))), drma_from_el_lt(degrange_from, degrange_to, lt, init_time, 0))
    
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ degree.mean.age(degree_vec, attrname, emptyval=2345.4))), rep(2345.4, length(degree_vec)*length(levels)))  
      expect_equal(unname(summary(nw ~ degree.mean.age(degree_vec, attrname))), drma_from_el_lt(degree_vec, degree_vec + 1, lt, init_time, 0, attr))
    
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ degrange.mean.age(degrange_from, degrange_to, attrname, emptyval=3345.4))), rep(3345.4, length(degrange_from)*length(levels)))  
      expect_equal(unname(summary(nw ~ degrange.mean.age(degrange_from, degrange_to, attrname))), drma_from_el_lt(degrange_from, degrange_to, lt, init_time, 0, attr))
    
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ degree.mean.age(degree_vec, emptyval=4345.4, log=TRUE))), rep(4345.4, length(degree_vec)))  
      expect_equal(unname(summary(nw ~ degree.mean.age(degree_vec, log=TRUE))), drma_from_el_lt(degree_vec, degree_vec + 1, lt, init_time, 0, log=TRUE))
    
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ degrange.mean.age(degrange_from, degrange_to, emptyval=5345.4, log=TRUE))), rep(5345.4, length(degrange_from)))  
      expect_equal(unname(summary(nw ~ degrange.mean.age(degrange_from, degrange_to, log=TRUE))), drma_from_el_lt(degrange_from, degrange_to, lt, init_time, 0, log=TRUE))
    
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ degree.mean.age(degree_vec, attrname, emptyval=6345.4, log=TRUE))), rep(6345.4, length(degree_vec)*length(levels)))  
      expect_equal(unname(summary(nw ~ degree.mean.age(degree_vec, attrname, log=TRUE))), drma_from_el_lt(degree_vec, degree_vec + 1, lt, init_time, 0, attr, log=TRUE))
    
      expect_equal(unname(summary(dynamic=TRUE, nw0 ~ degrange.mean.age(degrange_from, degrange_to, attrname, emptyval=7345.4, log=TRUE))), rep(7345.4, length(degrange_from)*length(levels)))  
      expect_equal(unname(summary(nw ~ degrange.mean.age(degrange_from, degrange_to, attrname, log=TRUE))), drma_from_el_lt(degrange_from, degrange_to, lt, init_time, 0, attr, log=TRUE))
    }
  
    toggle_list <- list()
    toggle_mat <- NULL
    
    last_ellt <- lt
    
    for(timestep in toggle_tmin:toggle_tmax) {
      if(timestep > toggle_tmin && timestep < toggle_tmax && runif(1) < skip_frac) next
      
      ## generate ~ netsize toggles randomly, about half edges and half non-edges
      inds <- sample(seq_len(NROW(last_ellt)), as.integer((sampling_frac*runif(1)*netsize + 1)), TRUE)
      
      toggles1 <- cbind(timestep, last_ellt[inds, 1:2, drop=FALSE])
      
      ndyads <- as.integer(runif(1)*netsize + 1)
      toggles2 <- matrix(timestep, nrow=ndyads, ncol=3)
      if(nwtype == 1) {
        for(i in seq_len(ndyads)) {
          tail <- as.integer(runif(1)*netsize + 1)
          head <- as.integer(runif(1)*(netsize - 1) + 1)
          if(head == tail) head <- netsize
          if(tail > head) {
            tmp <- tail
            tail <- head
            head <- tmp
          }
          toggles2[i,2] <- tail
          toggles2[i,3] <- head
        }
      } else if(nwtype == 2) {
        for(i in seq_len(ndyads)) {
          tail <- as.integer(runif(1)*bipsize + 1)
          head <- as.integer(bipsize + runif(1)*(netsize - bipsize - 1) + 1)
          toggles2[i,2] <- tail
          toggles2[i,3] <- head
        }    
      } else {
        for(i in seq_len(ndyads)) {
          tail <- as.integer(runif(1)*netsize + 1)
          head <- as.integer(runif(1)*(netsize - 1) + 1)
          if(head == tail) head <- netsize
          toggles2[i,2] <- tail
          toggles2[i,3] <- head
        }    
      }
      t_toggles <- rbind(toggles1, toggles2)
      t_toggles <- t_toggles[sample(seq_len(ndyads)),,drop=FALSE]
      toggle_list[[timestep]] <- t_toggles
      toggle_mat <- rbind(toggle_mat, t_toggles)
      last_ellt <- map_toggles_to_edgelist_lts(last_ellt, t_toggles)[[timestep]]
    }
  
    # toggle off all edges at time toggle_tmax + 1
    toggle_mat <- rbind(toggle_mat, cbind(toggle_tmax + 1, last_ellt[,1:2,drop=FALSE]))
  
    toggles <- toggle_mat
  
    if(nwtype == 1) {
      rv <- tergm.godfather(nw ~ edges + 
                                 mean.age(emptyval=3425.432) + 
                                 edge.ages + 
                                 edgecov.ages("edgewts") + 
                                 edgecov.mean.age("edgewts", emptyval=exp(pi/5)/3.14) + 
                                 mean.age(emptyval=13.1, log=TRUE) + 
                                 edgecov.mean.age("edgewts", emptyval=0.874, log=TRUE) + 
                                 edges.ageinterval(from=ages_from, to=ages_to) +
                                 nodefactor.mean.age(attr="mean_age_attr") +
                                 nodemix.mean.age(attr="mean_age_attr") +
                                 nodefactor.mean.age(attr="mean_age_attr", log=TRUE) +
                                 nodemix.mean.age(attr="mean_age_attr", log=TRUE) +
                                 nodefactor.mean.age(attr="mean_age_attr", levels=nodefactor_ma_lev) +
                                 nodemix.mean.age(attr="mean_age_attr", levels=nodemix_ma_lev,emptyval=nodemix_emptyvals) +
                                 nodefactor.mean.age(attr="mean_age_attr", log=TRUE, levels=nodefactor_ma_lev,emptyval=nodefactor_emptyvals) +
                                 nodemix.mean.age(attr="mean_age_attr", log=TRUE, levels=nodemix_ma_lev, levels2=nodemix_ma_lev2) +
                                 degree.mean.age(degree_vec, emptyval=0345.4) + 
                                 degrange.mean.age(degrange_from, degrange_to, emptyval=1345.4) + 
                                 degree.mean.age(degree_vec, attrname, emptyval=2345.4) + 
                                 degrange.mean.age(degrange_from, degrange_to, attrname, emptyval=3345.4) + 
                                 degree.mean.age(degree_vec, emptyval=4345.4, log=TRUE) + 
                                 degrange.mean.age(degrange_from, degrange_to, emptyval=5345.4, log=TRUE) + 
                                 degree.mean.age(degree_vec, attrname, emptyval=6345.4, log=TRUE) + 
                                 degrange.mean.age(degrange_from, degrange_to, attrname, emptyval=7345.4, log=TRUE) +
                                 Form(~edges) + 
                                 Persist(~edges),
                            toggles=toggles,
                            start=toggle_tmin-1L,
                            end=toggle_tmax+1L)
    } else if(nwtype == 2) {
      rv <- tergm.godfather(nw ~ edges + 
                                 mean.age(emptyval=3425.432) + 
                                 edge.ages + 
                                 edgecov.ages("edgewts") + 
                                 edgecov.mean.age("edgewts", emptyval=exp(pi/5)/3.14) + 
                                 mean.age(emptyval=13.1, log=TRUE) + 
                                 edgecov.mean.age("edgewts", emptyval=0.874, log=TRUE) + 
                                 edges.ageinterval(from=ages_from, to=ages_to) +
                                 nodefactor.mean.age(attr="mean_age_attr") +
                                 nodemix.mean.age(attr="mean_age_attr") +
                                 nodefactor.mean.age(attr="mean_age_attr", log=TRUE) +
                                 nodemix.mean.age(attr="mean_age_attr", log=TRUE) +
                                 nodefactor.mean.age(attr="mean_age_attr", levels=nodefactor_ma_lev) +
                                 nodemix.mean.age(attr="mean_age_attr", b1levels=nodemix_ma_b1lev, b2levels=nodemix_ma_b2lev,emptyval=nodemix_emptyvals) +
                                 nodefactor.mean.age(attr="mean_age_attr", log=TRUE, levels=nodefactor_ma_lev,emptyval=nodefactor_emptyvals) +
                                 nodemix.mean.age(attr="mean_age_attr", log=TRUE, b1levels=nodemix_ma_b1lev, b2levels=nodemix_ma_b2lev, levels2=nodemix_ma_lev2) +
                                 degree.mean.age(degree_vec, emptyval=0345.4) + 
                                 degrange.mean.age(degrange_from, degrange_to, emptyval=1345.4) + 
                                 degree.mean.age(degree_vec, attrname, emptyval=2345.4) + 
                                 degrange.mean.age(degrange_from, degrange_to, attrname, emptyval=3345.4) + 
                                 degree.mean.age(degree_vec, emptyval=4345.4, log=TRUE) + 
                                 degrange.mean.age(degrange_from, degrange_to, emptyval=5345.4, log=TRUE) + 
                                 degree.mean.age(degree_vec, attrname, emptyval=6345.4, log=TRUE) + 
                                 degrange.mean.age(degrange_from, degrange_to, attrname, emptyval=7345.4, log=TRUE) +
                                 Form(~edges) + 
                                 Persist(~edges),
                            toggles=toggles,
                            start=toggle_tmin-1L,
                            end=toggle_tmax+1L)    
    } else {
      rv <- tergm.godfather(nw ~ edges + 
                                 mean.age(emptyval=3425.432) + 
                                 edge.ages + 
                                 edgecov.ages("edgewts") + 
                                 edgecov.mean.age("edgewts", emptyval=exp(pi/5)/3.14) + 
                                 mean.age(emptyval=13.1, log=TRUE) + 
                                 edgecov.mean.age("edgewts", emptyval=0.874, log=TRUE) + 
                                 edges.ageinterval(from=ages_from, to=ages_to) +
                                 nodefactor.mean.age(attr="mean_age_attr") +
                                 nodemix.mean.age(attr="mean_age_attr") +
                                 nodefactor.mean.age(attr="mean_age_attr", log=TRUE) +
                                 nodemix.mean.age(attr="mean_age_attr", log=TRUE) +
                                 nodefactor.mean.age(attr="mean_age_attr", levels=nodefactor_ma_lev) +
                                 nodemix.mean.age(attr="mean_age_attr", levels=nodemix_ma_lev,emptyval=nodemix_emptyvals) +
                                 nodefactor.mean.age(attr="mean_age_attr", log=TRUE, levels=nodefactor_ma_lev,emptyval=nodefactor_emptyvals) +
                                 nodemix.mean.age(attr="mean_age_attr", log=TRUE, levels=nodemix_ma_lev, levels2=nodemix_ma_lev2) +
                                 Form(~edges) + 
                                 Persist(~edges),
                            toggles=toggles,
                            start=toggle_tmin-1L,
                            end=toggle_tmax+1L)  
    }
    
    rv <- as.matrix(rv)  
    
    el_lts <- map_toggles_to_edgelist_lts(lt, toggles)
    
    el_lts_b <- el_lts
    if(nwtype == 2) {
      for(j in toggle_tmin:toggle_tmax) {
        el_lts_b[[j]][,2] <- el_lts_b[[j]][,2] - bipsize
      }    
    }
    
    for(j in toggle_tmin:toggle_tmax) {
      this_iter_has_edges <- NROW(el_lts[[j]]) > 0
      expect_equal(if(this_iter_has_edges) mean(j - el_lts[[j]][,3] + 1) else 3425.432, unname(rv[j - init_time,2]))
      expect_equal(sum(j - el_lts[[j]][,3] + 1), unname(rv[j - init_time,3]))
      expect_equal(sum(wts[el_lts_b[[j]][,1:2,drop=FALSE]]*(j - el_lts[[j]][,3] + 1)), unname(rv[j - init_time, 4]))
      expect_equal(if(this_iter_has_edges) sum(wts[el_lts_b[[j]][,1:2,drop=FALSE]]*(j - el_lts[[j]][,3] + 1))/sum(wts[el_lts_b[[j]][,1:2,drop=FALSE]]) else exp(pi/5)/3.14, unname(rv[j - init_time, 5]))
      expect_equal(if(this_iter_has_edges) mean(log(j - el_lts[[j]][,3] + 1)) else 13.1, unname(rv[j - init_time,6]))
      expect_equal(if(this_iter_has_edges) sum(wts[el_lts_b[[j]][,1:2,drop=FALSE]]*log(j - el_lts[[j]][,3] + 1))/sum(wts[el_lts_b[[j]][,1:2,drop=FALSE]]) else 0.874, unname(rv[j - init_time, 7]))
      
      for(k in seq_len(length(ages_from))) {
        expect_equal(sum(j - el_lts[[j]][,3] + 1 >= ages_from[k] & j - el_lts[[j]][,3] + 1 < ages_to[k]), unname(rv[j - init_time, 7 + k]))
      }
  
      ## nodefactor/nodemix mean age tests
      nodefactor_ages <- numeric(n_nodefactor_stats)
      nodefactor_log_ages <- numeric(n_nodefactor_stats)
      nodefactor_counts <- integer(n_nodefactor_stats)
      
      attr_elt <- el_lts[[j]]
      attr_elt[,1] <- mean_age_attr[attr_elt[,1]]
      attr_elt[,2] <- mean_age_attr[attr_elt[,2]]
      for(k in seq_len(NROW(attr_elt))) {
        nodefactor_ages[attr_elt[k,1]] <- nodefactor_ages[attr_elt[k,1]] + j - attr_elt[k,3] + 1
        nodefactor_log_ages[attr_elt[k,1]] <- nodefactor_log_ages[attr_elt[k,1]] + log(j - attr_elt[k,3] + 1)
        nodefactor_counts[attr_elt[k,1]] <- nodefactor_counts[attr_elt[k,1]] + 1
        nodefactor_ages[attr_elt[k,2]] <- nodefactor_ages[attr_elt[k,2]] + j - attr_elt[k,3] + 1
        nodefactor_counts[attr_elt[k,2]] <- nodefactor_counts[attr_elt[k,2]] + 1
        nodefactor_log_ages[attr_elt[k,2]] <- nodefactor_log_ages[attr_elt[k,2]] + log(j - attr_elt[k,3] + 1)
      }
  
      nodefactor_vals <- numeric(n_nodefactor_stats)
      nodefactor_log_vals <- numeric(n_nodefactor_stats)
      for(k in seq_len(n_nodefactor_stats)) {
        if(nodefactor_counts[k] == 0) {
          ## else 0 for now; emptyval is handled later
        } else {
          nodefactor_vals[k] <- nodefactor_ages[k]/nodefactor_counts[k]
          nodefactor_log_vals[k] <- nodefactor_log_ages[k]/nodefactor_counts[k]
        }
      }
  
      nodemix_ages <- numeric(n_nodemix_stats)
      nodemix_log_ages <- numeric(n_nodemix_stats)
      nodemix_counts <- integer(n_nodemix_stats)
      
      for(k in seq_len(NROW(attr_elt))) {
        this_nodemix_index <- indmat[attr_elt[k,1],attr_elt[k,2]]
        nodemix_ages[this_nodemix_index] <- nodemix_ages[this_nodemix_index] + j - attr_elt[k,3] + 1
        nodemix_log_ages[this_nodemix_index] <- nodemix_log_ages[this_nodemix_index] + log(j - attr_elt[k,3] + 1)
        nodemix_counts[this_nodemix_index] <- nodemix_counts[this_nodemix_index] + 1
      }
  
      nodemix_vals <- numeric(n_nodemix_stats)
      nodemix_log_vals <- numeric(n_nodemix_stats)
      for(k in seq_len(n_nodemix_stats)) {
        if(nodemix_counts[k] == 0) {
          ## else 0 for now since we're not passing an emptyval
        } else {
          nodemix_vals[k] <- nodemix_ages[k]/nodemix_counts[k]
          nodemix_log_vals[k] <- nodemix_log_ages[k]/nodemix_counts[k]
        }
      }

      expect_equal(nodefactor_vals, unname(rv[j - init_time, 7 + length(ages_from) + seq_len(n_nodefactor_stats)]))
      expect_equal(nodemix_vals, unname(rv[j - init_time, 7 + length(ages_from) + n_nodefactor_stats + seq_len(n_nodemix_stats)]))
      expect_equal(nodefactor_log_vals, unname(rv[j - init_time, 7 + length(ages_from) + n_nodefactor_stats + n_nodemix_stats + seq_len(n_nodefactor_stats)]))
      expect_equal(nodemix_log_vals, unname(rv[j - init_time, 7 + length(ages_from) + 2*n_nodefactor_stats + n_nodemix_stats + seq_len(n_nodemix_stats)]))

      nodemix_vals[nodemix_indices_1][nodemix_counts[nodemix_indices_1] == 0] <- nodemix_emptyvals[nodemix_counts[nodemix_indices_1] == 0]
      nodefactor_log_vals[nodefactor_ma_lev][nodefactor_counts[nodefactor_ma_lev] == 0] <- nodefactor_emptyvals[nodefactor_counts[nodefactor_ma_lev] == 0]

      expect_equal(nodefactor_vals[nodefactor_ma_lev], unname(rv[j - init_time, 7 + length(ages_from) + 2*n_nodefactor_stats + 2*n_nodemix_stats + seq_along(nodefactor_ma_lev)]))
      expect_equal(nodemix_vals[nodemix_indices_1], unname(rv[j - init_time, 7 + length(ages_from) + 2*n_nodefactor_stats + 2*n_nodemix_stats + length(nodefactor_ma_lev) + seq_along(nodemix_indices_1)]))
      expect_equal(nodefactor_log_vals[nodefactor_ma_lev], unname(rv[j - init_time, 7 + length(ages_from) + 2*n_nodefactor_stats + 2*n_nodemix_stats + length(nodefactor_ma_lev) + length(nodemix_indices_1) + seq_along(nodefactor_ma_lev)]))
      expect_equal(nodemix_log_vals[nodemix_indices_2], unname(rv[j - init_time, 7 + length(ages_from) + 2*n_nodefactor_stats + 2*n_nodemix_stats + 2*length(nodefactor_ma_lev) + length(nodemix_indices_1) + seq_along(nodemix_indices_2)]))

      if(nwtype < 3) {  
        last_index <- 7 + 2*n_nodefactor_stats + 2*length(nodefactor_ma_lev) + 2*n_nodemix_stats + nstats_nodemix_no2 + length(nodemix_ma_lev2) + length(ages_from)
    
        expect_equal(drma_from_el_lt(degree_vec, degree_vec + 1, el_lts[[j]], j, 0345.4), unname(rv[j - init_time, (last_index + 1):(last_index + length(degree_vec))]))
        last_index <- last_index + length(degree_vec)
        expect_equal(drma_from_el_lt(degrange_from, degrange_to, el_lts[[j]], j, 1345.4), unname(rv[j - init_time, (last_index + 1):(last_index + length(degrange_from))]))
        last_index <- last_index + length(degrange_from)
        expect_equal(drma_from_el_lt(degree_vec, degree_vec + 1, el_lts[[j]], j, 2345.4, attr), unname(rv[j - init_time, (last_index + 1):(last_index + length(degree_vec)*length(levels))]))
        last_index <- last_index + length(degree_vec)*length(levels)
        expect_equal(drma_from_el_lt(degrange_from, degrange_to, el_lts[[j]], j, 3345.4, attr), unname(rv[j - init_time, (last_index + 1):(last_index + length(degrange_from)*length(levels))]))
        last_index <- last_index + length(degrange_from)*length(levels)
    
        expect_equal(drma_from_el_lt(degree_vec, degree_vec + 1, el_lts[[j]], j, 4345.4, log=TRUE), unname(rv[j - init_time, (last_index + 1):(last_index + length(degree_vec))]))
        last_index <- last_index + length(degree_vec)
        expect_equal(drma_from_el_lt(degrange_from, degrange_to, el_lts[[j]], j, 5345.4, log=TRUE), unname(rv[j - init_time, (last_index + 1):(last_index + length(degrange_from))]))
        last_index <- last_index + length(degrange_from)
        expect_equal(drma_from_el_lt(degree_vec, degree_vec + 1, el_lts[[j]], j, 6345.4, attr, log=TRUE), unname(rv[j - init_time, (last_index + 1):(last_index + length(degree_vec)*length(levels))]))
        last_index <- last_index + length(degree_vec)*length(levels)
        expect_equal(drma_from_el_lt(degrange_from, degrange_to, el_lts[[j]], j, 7345.4, attr, log=TRUE), unname(rv[j - init_time, (last_index + 1):(last_index + length(degrange_from)*length(levels))]))
      }
    }
  
    expect_equal(3425.432, unname(rv[toggle_tmax + 1 - init_time,2]))  
    expect_equal(0, unname(rv[toggle_tmax + 1 - init_time,3]))
    expect_equal(0, unname(rv[toggle_tmax + 1 - init_time,4]))
    expect_equal(exp(pi/5)/3.14, unname(rv[toggle_tmax + 1 - init_time,5]))
    expect_equal(13.1, unname(rv[toggle_tmax + 1 - init_time,6]))
    expect_equal(0.874, unname(rv[toggle_tmax + 1 - init_time,7]))
    expect_equal(integer(length(ages_from)), unname(rv[toggle_tmax + 1 - init_time, 8:(8 + length(ages_from) - 1)]))
  
    expect_equal(rep(0,n_nodefactor_stats), unname(rv[toggle_tmax + 1 - init_time, 7 + length(ages_from) + seq_len(n_nodefactor_stats)]))
    expect_equal(rep(0,n_nodemix_stats), unname(rv[toggle_tmax + 1 - init_time, 7 + length(ages_from) + n_nodefactor_stats + seq_len(n_nodemix_stats)]))
    expect_equal(rep(0,n_nodefactor_stats), unname(rv[toggle_tmax + 1 - init_time, 7 + length(ages_from) + n_nodefactor_stats + n_nodemix_stats + seq_len(n_nodefactor_stats)]))
    expect_equal(rep(0,n_nodemix_stats), unname(rv[toggle_tmax + 1 - init_time, 7 + length(ages_from) + 2*n_nodefactor_stats + n_nodemix_stats + seq_len(n_nodemix_stats)]))
  
    expect_equal(rep(0, length(nodefactor_ma_lev)), unname(rv[toggle_tmax + 1 - init_time, 7 + length(ages_from) + 2*n_nodefactor_stats + 2*n_nodemix_stats + seq_along(nodefactor_ma_lev)]))
    expect_equal(nodemix_emptyvals, unname(rv[toggle_tmax + 1 - init_time, 7 + length(ages_from) + 2*n_nodefactor_stats + 2*n_nodemix_stats + length(nodefactor_ma_lev) + seq_along(nodemix_indices_1)]))
    expect_equal(nodefactor_emptyvals, unname(rv[toggle_tmax + 1 - init_time, 7 + length(ages_from) + 2*n_nodefactor_stats + 2*n_nodemix_stats + length(nodefactor_ma_lev) + length(nodemix_indices_1) + seq_along(nodefactor_ma_lev)]))
    expect_equal(rep(0, length(nodemix_indices_2)), unname(rv[toggle_tmax + 1 - init_time, 7 + length(ages_from) + 2*n_nodefactor_stats + 2*n_nodemix_stats + 2*length(nodefactor_ma_lev) + length(nodemix_indices_1) + seq_along(nodemix_indices_2)]))
  
    if(nwtype < 3) {  
      last_index <- 7 + 2*n_nodefactor_stats + 2*length(nodefactor_ma_lev) + 2*n_nodemix_stats + nstats_nodemix_no2 + length(nodemix_ma_lev2) + length(ages_from)

      expect_equal(rep(0345.4, length(degree_vec)), unname(rv[toggle_tmax + 1 - init_time, (last_index + 1):(last_index + length(degree_vec))]))
      last_index <- last_index + length(degree_vec)
      expect_equal(rep(1345.4, length(degrange_from)), unname(rv[toggle_tmax + 1 - init_time, (last_index + 1):(last_index + length(degrange_from))]))
      last_index <- last_index + length(degrange_from)
      expect_equal(rep(2345.4, length(degree_vec)*length(levels)), unname(rv[toggle_tmax + 1 - init_time, (last_index + 1):(last_index + length(degree_vec)*length(levels))]))
      last_index <- last_index + length(degree_vec)*length(levels)
      expect_equal(rep(3345.4, length(degrange_from)*length(levels)), unname(rv[toggle_tmax + 1 - init_time, (last_index + 1):(last_index + length(degrange_from)*length(levels))]))
      last_index <- last_index + length(degrange_from)*length(levels)
    
      expect_equal(rep(4345.4, length(degree_vec)), unname(rv[toggle_tmax + 1 - init_time, (last_index + 1):(last_index + length(degree_vec))]))
      last_index <- last_index + length(degree_vec)
      expect_equal(rep(5345.4, length(degrange_from)), unname(rv[toggle_tmax + 1 - init_time, (last_index + 1):(last_index + length(degrange_from))]))
      last_index <- last_index + length(degrange_from)
      expect_equal(rep(6345.4, length(degree_vec)*length(levels)), unname(rv[toggle_tmax + 1 - init_time, (last_index + 1):(last_index + length(degree_vec)*length(levels))]))
      last_index <- last_index + length(degree_vec)*length(levels)
      expect_equal(rep(7345.4, length(degrange_from)*length(levels)), unname(rv[toggle_tmax + 1 - init_time, (last_index + 1):(last_index + length(degrange_from)*length(levels))]))
    }
  }
})


test_that("edges.ageinterval behaves correctly inside a dissolution operator", {
  set.seed(0)

  logit <- function(p) log(p/(1-p))
  expit <- function(x) 1/(1+exp(-x))

  T <- 1000

  nw <- network(16, directed=FALSE)
  ndyads <- network.dyadcount(nw)

  nwd <- simulate(nw~Form(~edges)+Persist(~edges+edges.ageinterval(3,7)), dynamic=TRUE, output="networkDynamic", coef=c(-2,2,-1/2), time.slices=T, seed=0)
  spells <- as.data.frame(nwd)

  # Test dissolution hazards
  freq <- tabulate(spells$duration)
  surv <- rev(cumsum(rev(freq)))

  haz <- freq/surv

  basehaz <- weighted.mean(haz[-(3:6)], surv[-(3:6)])
  addhaz <- weighted.mean(haz[3:3], surv[3:3])

  expect_equal(basehaz, 1-expit(2), tolerance=0.05)
  expect_equal(addhaz, 1-expit(2-1/2), tolerance=0.15)

  pform <- sapply(1:T, function(t){
    eid0 <- spells$edge.id[spells$onset<=t-1 & t-1<spells$terminus]
    eid1 <- spells$edge.id[spells$onset<=t & t<spells$terminus]
    sum(! eid1%in%eid0)/(ndyads-length(eid0))
  })

  expect_equal(mean(pform), expit(-2), tolerance=.05)

  nwd.diss <- simulate(nw~Form(~edges)+Diss(~edges+edges.ageinterval(3,7)), dynamic=TRUE, output="networkDynamic", coef=c(-2,-2,+1/2), time.slices=T, seed=0)
  expect_equal(nwd, nwd.diss, ignore_attr=TRUE)
})
