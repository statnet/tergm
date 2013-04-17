#  File R/summary.statistics.networkDynamic.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2013 Statnet Commons
#######################################################################
summary.statistics.networkDynamic <- function(object, at,..., basis=NULL){
  if(!is.null(basis)) object <- ergm.update.formula(object, basis~.)
  t(rbind(sapply(at,
                 function(t){
                   nw <- network.extract.with.lasttoggle(ergm.getnetwork(object), t)
                   f <- ergm.update.formula(object, nw~., from.new="nw")
                   summary(f,...)
                 }
                 )
          )
    )
}
