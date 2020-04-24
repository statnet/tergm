#  File tests/tergm_offset_tests.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2019 Statnet Commons
#######################################################################
library(statnet.common)
#opttest({
  library(tergm)
  data(florentine)
  net <- flobusiness
  summary(net ~ edges+degree(3))
  set.seed(3)
  
  # default initialization
  mod1 <- tergm(flobusiness ~ FormE(~edges) + offset(FormE(~degree(3))) + offset(DissE(~edges)), 
                 offset.coef=c(0.8,log(9)), 
                 targets="formation",
                 estimate="EGMME"
  )
  
  # init.method set to zeros works
  mod2 <- tergm(flobusiness ~ FormE(~edges + offset(degree(3))) + offset(DissE(~edges)), 
                 offset.coef=c(0.8,log(9)), 
                 targets="formation",
                 estimate="EGMME",control=control.tergm(init.method='zeros'))
  
  # this works, auto defaulting SAN coefs if they are different from stergm init.form
  mod3 <- tergm(flobusiness ~ FormE(~edges) + FormE(~offset(degree(3))) + offset(DissE(~edges)), 
                 offset.coef=c(0.8,log(9)), 
                 targets="formation", estimate="EGMME",   
                 control=control.tergm(init=c(-3,0.8,log(9)))
  )
  
  # we can explicitly specify the target and target stats
  mod4 <- tergm(flobusiness ~ FormE(~edges) + offset(FormE(~degree(3))) + offset(DissE(~edges)), 
                 offset.coef=c(0.8,log(9)), 
                 targets=~edges, target.stats = 15,
                 estimate="EGMME",   
                 control=control.tergm(init=c(-3,0.8,log(9)))
  )
  
  sapply(list(mod1, mod2, mod3, mod4), function(x) x$fit$coef)
  
  for (mod in list(mod1, mod2, mod3, mod4)) {
    print(apply(simulate(mod, monitor=~edges+degree(3), output="stats", time.slices=200), 2, mean))
  }
#}, testname='target_offset')
