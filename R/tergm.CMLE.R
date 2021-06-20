#  File R/tergm.CMLE.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################
tergm.CMLE <- function(formula, times, ..., control) {
  nw <- eval_lhs.formula(formula)

  if(!is(nw, "tergm_NetSeries")){
    if(inherits(nw, "network.list") || (is.list(nw) && !is.network(nw) && is.network(nw[[1]]))){
      NetSeries <- NetSeries(nw, NA.impute=control$CMLE.NA.impute)
    }else if(inherits(nw,"networkDynamic")){
      NetSeries <- NetSeries(nw, times, NA.impute=control$CMLE.NA.impute)
    }else{
      stop("Unsupported specification for the network series. See help for ",sQuote("NetSeries")," for arguments.")
    }

    formula <- nonsimp_update.formula(formula, NetSeries~., from.new="NetSeries")
  }

  fit <- ergm(formula, ..., control=control$CMLE.ergm)
  class(fit) <- c("tergm_CMLE", "tergm", class(fit))
  fit
}
