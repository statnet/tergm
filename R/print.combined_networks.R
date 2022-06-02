#  File R/print.combined_networks.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2022 Statnet Commons
################################################################################
#' @describeIn combine_networks A wrapper around
#'   [network::print.network()] to print constituent network
#'   information and omit some internal variables.
#'
#' @param x,object a combined network.
#' @param ... additional arguments to methods.
#' @export
print.combined_networks<-function(x, ...) {
  cat(" Combined ", length((x%n%".subnetattr")[[1]]$n), " networks on ", sQuote(names(x$gal$.subnetattr)), ":\n", sep="")
  nattrs <- (x%n%".subnetattr")[[1]]
  nids <- format(seq_along(nattrs$n))
  for(i in seq_along(nattrs$n)){
    cat("  ", nids[[i]],": n = ", nattrs$n[[i]], ", directed = ", nattrs$directed[[i]], ", bipartite = ", nattrs$bipartite[[i]], ", loops = ", nattrs$loops[[i]], "\n", sep="")
  }
  cat("\n")

  delete.network.attribute(x, c(".subnetcache", ".subnetattr"))
  NextMethod()
}

#' @describeIn combine_networks A wrapper around
#'   [network::summary.network()] to print constituent network
#'   information and omit some internal variables.
#' @export
summary.combined_networks<-function (object, ...) {
  object <- NextMethod(object)
  structure(object, class = c("summary.combined_networks", class(object)))
}

#' @describeIn combine_networks A wrapper around
#'   [network::print.summary.network()] to print constituent network
#'   information and omit some internal variables.
#' @export
print.summary.combined_networks<-function(x, ...) {
  cat("Combined ", length((x%n%".subnetattr")[[1]]$n), " networks on ", sQuote(names(x$gal$.subnetattr)), ":\n", sep="")
  nattrs <- (x%n%".subnetattr")[[1]]
  nids <- format(seq_along(nattrs$n))
  for(i in seq_along(nattrs$n)){
    cat("  ", nids[[i]],": n = ", nattrs$n[[i]], ", directed = ", nattrs$directed[[i]], ", bipartite = ", nattrs$bipartite[[i]], ", loops = ", nattrs$loops[[i]], "\n", sep="")
  }
  cat("\n")

  delete.network.attribute(x, c(".subnetcache", ".subnetattr"))
  NextMethod()
}

