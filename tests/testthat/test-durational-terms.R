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

test_that("mean.age, edge.ages, edgecov.ages, and edgecov.mean.age behave correctly with summary and godfather", {
  el <- matrix(c(1, 2,
                 1, 3,
                 2, 3,
                 5, 6,
                 7, 8,
                 7, 10), ncol = 2, byrow = TRUE)
  
  init.lasttoggle.times <- c(0,1,2,3,3,5)
  
  lt <- cbind(el, init.lasttoggle.times)
  
  nw <- network(el, dir = FALSE)
  
  nw %n% "time" <- 5
  
  nw %n% "lasttoggle" <- lt
  
  wts <- matrix(runif(100),10,10)
  wts <- wts + t(wts)
  
  nw %n% "edgewts" <- wts
  
  nw0 <- network.initialize(10,dir=FALSE)
  
  expect_equal(unname(summary(nw0 ~ mean.age(emptyval=4321.4))), 4321.4)  
  expect_equal(unname(summary(nw ~ mean.age)), mean(5 - init.lasttoggle.times + 1))
  
  expect_equal(unname(summary(nw0 ~ edge.ages)), 0)  
  expect_equal(unname(summary(nw ~ edge.ages)), sum(5 - init.lasttoggle.times + 1))
  
  expect_equal(unname(summary(nw0 ~ edgecov.ages("edgewts"))), 0)  
  expect_equal(unname(summary(nw ~ edgecov.ages("edgewts"))), sum(wts[el]*(5 - init.lasttoggle.times + 1)))
  
  expect_equal(unname(summary(nw0 ~ edgecov.mean.age("edgewts", emptyval=pi^4/25))), pi^4/25)  
  expect_equal(unname(summary(nw ~ edgecov.mean.age("edgewts"))), sum(wts[el]*(5 - init.lasttoggle.times + 1))/sum(wts[el]))
    
  expect_equal(unname(summary(nw0 ~ mean.age(emptyval=pi^3/83, log=TRUE))), pi^3/83)    
  expect_equal(unname(summary(nw ~ mean.age(log=TRUE))), mean(log(5 - init.lasttoggle.times + 1)))

  expect_equal(unname(summary(nw0 ~ edgecov.mean.age("edgewts", emptyval=pi^4/25, log=TRUE))), pi^4/25)  
  expect_equal(unname(summary(nw ~ edgecov.mean.age("edgewts", log=TRUE))), sum(wts[el]*log(5 - init.lasttoggle.times + 1))/sum(wts[el]))
    
  toggles <- matrix(c(6, 2, 3,
                      7, 1, 4,
                      7, 1, 5,
                      8, 7, 8,
                      8, 7, 10,
                      8, 9, 10,
                      9, 1, 4,
                      9, 5, 6,
                      10, 5, 6,
                      10, 5, 6,
                      11, 1, 2,
                      11, 1, 2,
                      12, 5, 8,
                      12, 5, 8,
                      12, 5, 8,
                      13, 1, 3,
                      13, 1, 3,
                      13, 1, 3,
                      15, 1, 2,
                      15, 1, 5,
                      15, 5, 8,
                      15, 9, 10), ncol = 3, byrow = TRUE)
  
  rv <- tergm.godfather(nw ~ edges + mean.age(emptyval=3425.432) + edge.ages + edgecov.ages("edgewts") + edgecov.mean.age("edgewts", emptyval=exp(pi/5)/3.14) + mean.age(emptyval=13.1, log=TRUE) + edgecov.mean.age("edgewts", emptyval=0.874, log=TRUE), toggles=toggles)
  rv <- as.matrix(rv)
  
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
  
  el_lts <- map_toggles_to_edgelist_lts(lt, toggles)
  
  for(j in 6:14) {
    expect_equal(mean(j - el_lts[[j]][,3] + 1), unname(rv[j - 5,2]))
    expect_equal(sum(j - el_lts[[j]][,3] + 1), unname(rv[j - 5,3]))
    expect_equal(sum(wts[el_lts[[j]][,1:2]]*(j - el_lts[[j]][,3] + 1)), unname(rv[j - 5, 4]))
    expect_equal(sum(wts[el_lts[[j]][,1:2]]*(j - el_lts[[j]][,3] + 1))/sum(wts[el_lts[[j]][,1:2]]), unname(rv[j - 5, 5]))
    expect_equal(mean(log(j - el_lts[[j]][,3] + 1)), unname(rv[j - 5,6]))
    expect_equal(sum(wts[el_lts[[j]][,1:2]]*log(j - el_lts[[j]][,3] + 1))/sum(wts[el_lts[[j]][,1:2]]), unname(rv[j - 5, 7]))
  }

  expect_equal(3425.432, unname(rv[10,2]))  
  expect_equal(0, unname(rv[10,3]))
  expect_equal(0, unname(rv[10,4]))
  expect_equal(exp(pi/5)/3.14, unname(rv[10,5]))
  expect_equal(13.1, unname(rv[10,6]))
  expect_equal(0.874, unname(rv[10,7]))
})
