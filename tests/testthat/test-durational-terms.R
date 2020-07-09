#  File tests/testthat/test-durational-terms.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2003-2020 Statnet Commons
#######################################################################

context("test-durational-terms.R")

test_that("durational terms behave correctly with summary and godfather", {
  
  netsize <- 100

  init_tmin <- 0
  init_tmax <- 10
  init_time <- init_tmax

  skip_frac <- 0.3

  sampling_frac <- 0.8

  toggle_tmin <- 11
  toggle_tmax <- 30

  max_levels <- 5

  nw <- network.initialize(netsize,dir=FALSE)
  nw <- san(nw ~ edges, target.stats = c(3*netsize))
  
  el <- as.edgelist(nw)

  attrname <- "race"

  attr <- sample(letters[1:max_levels], netsize, TRUE)
  levels <- sort(unique(attr))

  nw %v% attrname <- attr
  
  init.lasttoggle.times <- sample(init_tmin:init_tmax, NROW(el), TRUE)
  
  lt <- cbind(el, init.lasttoggle.times)
    
  nw %n% "time" <- init_time  
  nw %n% "lasttoggle" <- lt
  
  wts <- matrix(runif(netsize*netsize),netsize,netsize)
  wts <- wts + t(wts)
  
  nw %n% "edgewts" <- wts
  
  nw0 <- network.initialize(netsize,dir=FALSE)
  nw0 %v% attrname <- attr
  
  ages_from = c(1, 1, 4, 2, 5, 10, 1, 5)
  ages_to = c(2, 3, 15, 5, 7, 12, Inf, Inf) 
  
  init_ages <- init_time - init.lasttoggle.times + 1
  
  init_age_counts <- integer(length(ages_from))
  for(i in seq_len(length(ages_from))) {
    init_age_counts[i] <- sum(init_ages >= ages_from[i] & init_ages < ages_to[i])
  }

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

  expect_equal(unname(summary(nw0 ~ mean.age(emptyval=4321.4))), 4321.4)  
  expect_equal(unname(summary(nw ~ mean.age)), mean(init_time - init.lasttoggle.times + 1))
  
  expect_equal(unname(summary(nw0 ~ edge.ages)), 0)  
  expect_equal(unname(summary(nw ~ edge.ages)), sum(init_time - init.lasttoggle.times + 1))
  
  expect_equal(unname(summary(nw0 ~ edgecov.ages("edgewts"))), 0)  
  expect_equal(unname(summary(nw ~ edgecov.ages("edgewts"))), sum(wts[el]*(init_time - init.lasttoggle.times + 1)))
  
  expect_equal(unname(summary(nw0 ~ edgecov.mean.age("edgewts", emptyval=pi^4/25))), pi^4/25)  
  expect_equal(unname(summary(nw ~ edgecov.mean.age("edgewts"))), sum(wts[el]*(init_time - init.lasttoggle.times + 1))/sum(wts[el]))
    
  expect_equal(unname(summary(nw0 ~ mean.age(emptyval=pi^3/83, log=TRUE))), pi^3/83)    
  expect_equal(unname(summary(nw ~ mean.age(log=TRUE))), mean(log(init_time - init.lasttoggle.times + 1)))

  expect_equal(unname(summary(nw0 ~ edgecov.mean.age("edgewts", emptyval=pi^4/25, log=TRUE))), pi^4/25)  
  expect_equal(unname(summary(nw ~ edgecov.mean.age("edgewts", log=TRUE))), sum(wts[el]*log(init_time - init.lasttoggle.times + 1))/sum(wts[el]))
    
  expect_equal(unname(summary(nw0 ~ edges.ageinterval(from=ages_from, to=ages_to))), integer(length(ages_from)))  
  expect_equal(unname(summary(nw ~ edges.ageinterval(from=ages_from, to=ages_to))), init_age_counts)  

  expect_equal(unname(summary(nw0 ~ degree.mean.age(degree_vec, emptyval=0345.4))), rep(0345.4, length(degree_vec)))  
  expect_equal(unname(summary(nw ~ degree.mean.age(degree_vec))), drma_from_el_lt(degree_vec, degree_vec + 1, lt, init_time, 0))

  expect_equal(unname(summary(nw0 ~ degrange.mean.age(degrange_from, degrange_to, emptyval=1345.4))), rep(1345.4, length(degrange_from)))  
  expect_equal(unname(summary(nw ~ degrange.mean.age(degrange_from, degrange_to))), drma_from_el_lt(degrange_from, degrange_to, lt, init_time, 0))

  expect_equal(unname(summary(nw0 ~ degree.mean.age(degree_vec, attrname, emptyval=2345.4))), rep(2345.4, length(degree_vec)*length(levels)))  
  expect_equal(unname(summary(nw ~ degree.mean.age(degree_vec, attrname))), drma_from_el_lt(degree_vec, degree_vec + 1, lt, init_time, 0, attr))

  expect_equal(unname(summary(nw0 ~ degrange.mean.age(degrange_from, degrange_to, attrname, emptyval=3345.4))), rep(3345.4, length(degrange_from)*length(levels)))  
  expect_equal(unname(summary(nw ~ degrange.mean.age(degrange_from, degrange_to, attrname))), drma_from_el_lt(degrange_from, degrange_to, lt, init_time, 0, attr))

  expect_equal(unname(summary(nw0 ~ degree.mean.age(degree_vec, emptyval=4345.4, log=TRUE))), rep(4345.4, length(degree_vec)))  
  expect_equal(unname(summary(nw ~ degree.mean.age(degree_vec, log=TRUE))), drma_from_el_lt(degree_vec, degree_vec + 1, lt, init_time, 0, log=TRUE))

  expect_equal(unname(summary(nw0 ~ degrange.mean.age(degrange_from, degrange_to, emptyval=5345.4, log=TRUE))), rep(5345.4, length(degrange_from)))  
  expect_equal(unname(summary(nw ~ degrange.mean.age(degrange_from, degrange_to, log=TRUE))), drma_from_el_lt(degrange_from, degrange_to, lt, init_time, 0, log=TRUE))

  expect_equal(unname(summary(nw0 ~ degree.mean.age(degree_vec, attrname, emptyval=6345.4, log=TRUE))), rep(6345.4, length(degree_vec)*length(levels)))  
  expect_equal(unname(summary(nw ~ degree.mean.age(degree_vec, attrname, log=TRUE))), drma_from_el_lt(degree_vec, degree_vec + 1, lt, init_time, 0, attr, log=TRUE))

  expect_equal(unname(summary(nw0 ~ degrange.mean.age(degrange_from, degrange_to, attrname, emptyval=7345.4, log=TRUE))), rep(7345.4, length(degrange_from)*length(levels)))  
  expect_equal(unname(summary(nw ~ degrange.mean.age(degrange_from, degrange_to, attrname, log=TRUE))), drma_from_el_lt(degrange_from, degrange_to, lt, init_time, 0, attr, log=TRUE))


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
    t_toggles <- rbind(toggles1, toggles2)
    t_toggles <- t_toggles[sample(seq_len(ndyads)),,drop=FALSE]
    toggle_list[[timestep]] <- t_toggles
    toggle_mat <- rbind(toggle_mat, t_toggles)
    last_ellt <- map_toggles_to_edgelist_lts(last_ellt, t_toggles)[[timestep]]
  }

  # toggle off all edges at time toggle_tmax + 1
  toggle_mat <- rbind(toggle_mat, cbind(toggle_tmax + 1, last_ellt[,1:2]))

  toggles <- toggle_mat

  rv <- tergm.godfather(nw ~ edges + 
                             mean.age(emptyval=3425.432) + 
                             edge.ages + 
                             edgecov.ages("edgewts") + 
                             edgecov.mean.age("edgewts", emptyval=exp(pi/5)/3.14) + 
                             mean.age(emptyval=13.1, log=TRUE) + 
                             edgecov.mean.age("edgewts", emptyval=0.874, log=TRUE) + 
                             edges.ageinterval(from=ages_from, to=ages_to) +
                             degree.mean.age(degree_vec, emptyval=0345.4) + 
                             degrange.mean.age(degrange_from, degrange_to, emptyval=1345.4) + 
                             degree.mean.age(degree_vec, attrname, emptyval=2345.4) + 
                             degrange.mean.age(degrange_from, degrange_to, attrname, emptyval=3345.4) + 
                             degree.mean.age(degree_vec, emptyval=4345.4, log=TRUE) + 
                             degrange.mean.age(degrange_from, degrange_to, emptyval=5345.4, log=TRUE) + 
                             degree.mean.age(degree_vec, attrname, emptyval=6345.4, log=TRUE) + 
                             degrange.mean.age(degrange_from, degrange_to, attrname, emptyval=7345.4, log=TRUE) +
                             Form(~edges) + 
                             Diss(~edges),
                        toggles=toggles)

  rv <- as.matrix(rv)
  
  
  el_lts <- map_toggles_to_edgelist_lts(lt, toggles)
  
  # note that the averages will not function correctly if we hit an empty network,
  # but parameters have been chosen such that we never should, except when it is
  # forced on the final time step, which we test differently after this loop
  for(j in toggle_tmin:toggle_tmax) {
    expect_equal(mean(j - el_lts[[j]][,3] + 1), unname(rv[j - init_time,2]))
    expect_equal(sum(j - el_lts[[j]][,3] + 1), unname(rv[j - init_time,3]))
    expect_equal(sum(wts[el_lts[[j]][,1:2]]*(j - el_lts[[j]][,3] + 1)), unname(rv[j - init_time, 4]))
    expect_equal(sum(wts[el_lts[[j]][,1:2]]*(j - el_lts[[j]][,3] + 1))/sum(wts[el_lts[[j]][,1:2]]), unname(rv[j - init_time, 5]))
    expect_equal(mean(log(j - el_lts[[j]][,3] + 1)), unname(rv[j - init_time,6]))
    expect_equal(sum(wts[el_lts[[j]][,1:2]]*log(j - el_lts[[j]][,3] + 1))/sum(wts[el_lts[[j]][,1:2]]), unname(rv[j - init_time, 7]))
    
    for(k in seq_len(length(ages_from))) {
      expect_equal(sum(j - el_lts[[j]][,3] + 1 >= ages_from[k] & j - el_lts[[j]][,3] + 1 < ages_to[k]), unname(rv[j - init_time, 7 + k]))
    }
  
    last_index <- 7 + length(ages_from)

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

  expect_equal(3425.432, unname(rv[toggle_tmax + 1 - init_time,2]))  
  expect_equal(0, unname(rv[toggle_tmax + 1 - init_time,3]))
  expect_equal(0, unname(rv[toggle_tmax + 1 - init_time,4]))
  expect_equal(exp(pi/5)/3.14, unname(rv[toggle_tmax + 1 - init_time,5]))
  expect_equal(13.1, unname(rv[toggle_tmax + 1 - init_time,6]))
  expect_equal(0.874, unname(rv[toggle_tmax + 1 - init_time,7]))
  expect_equal(integer(length(ages_from)), unname(rv[toggle_tmax + 1 - init_time, 8:(8 + length(ages_from) - 1)]))


  last_index <- 7 + length(ages_from)
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
})

