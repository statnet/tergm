#  File R/InitErgmTerm.netseries.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################

#' A network series specification for conditional modeling.
#'
#' A function for specifying the LHS of a temporal network series ERGM.
#'
#' @param ... series specification, in one of three formats:
#' 
#'   1. A list of identically- dimensioned and directed networks.
#'
#'   1. Several networks as arguments.
#'
#'   1. A [`networkDynamic`] object and a numeric vector of time indices.
#'
#' @param order how many previous networks to store as an accessible
#'   covariate of the model.
#'
#' @param NA.impute How missing dyads in transitioned-from networks
#'   are be imputed when using conditional estimation. See argument
#'   \code{imputers} of [impute.network.list()] for details.
#'
#' @return A network object with temporal metadata.
#'
#' @note It is not recommended to modify the network returned by
#'   `NetSeries` except by adding and removing edges, and even that
#'   must be done with some care, to avoid putting it into an
#'   inconsistent state.
#'
#'   It is almost always better to modify the original networks and
#'   regenerate the series.
#'
#' @seealso [`ergmTerm`] for specific terms.
#' 
#' @examples
#'
#' data(samplk)
#'
#' # Method 1: list of networks
#' monks <- NetSeries(list(samplk1,samplk2,samplk3))
#' ergm(monks ~ Form(~edges)+Diss(~edges))
#' ergm(monks ~ Form(~edges)+Persist(~edges))
#'
#' # Method 2: networks as arguments
#' monks <- NetSeries(samplk1,samplk2,samplk3)
#' ergm(monks ~ Form(~edges)+Diss(~edges))
#' ergm(monks ~ Form(~edges)+Persist(~edges))
#'
#' # Method 3: networkDynamic and time points:
#' ## TODO
#'
#' @export
NetSeries <- function(..., order=1, NA.impute=NULL){
  if(order>1) stop("Higher-order network models are not supported at this time.")
  
  args <- list(...)
  #' @importFrom methods is
  if(is(args[[1]], "networkDynamic")){
    nw <- args[[1]]
    times <- if(length(args)>=2) args[[2]]
             else{warning("the times argument was not provided to specify sampling time points for networkDynamic object. Modeling transition from time 0 to 1."); c(0,1)}
    # Grab only the needed vertices.
    subnw <- network.extract(nw, onset=min(times), terminus=max(times)+min(abs(diff(times)))/2)
    # Grab the vector of vertex activity indicators for each time
    # point, bind them into an n*T matrix. If any rows have
    # variability, we have a changing composition.
    if(any(apply(sapply(times, function(t) networkDynamic::is.active(subnw, at=t, v=1:network.size(subnw))), 1, var)>0)) warning("Active vertex set varies from time point to time point. Estimation may not work.")
    
    nwl <- lapply(times, function(t) networkDynamic::network.collapse(subnw, at=t, retain.all.vertices=TRUE))
  }else if(all(sapply(args, is, "network"))){
    nwl <- args
    times <- seq_along(nwl)
  }else if(is.list(args[[1]]) && all(sapply(args[[1]], is, "network"))){
    nwl <- args[[1]]
    times <- seq_along(nwl)
  }else stop("Unrecognized format for network series specification. See help for information.")

  nwl0 <- impute.network.list(nwl[-length(nwl)], NA.impute, nwl.append=nwl[length(nwl)])

  nwl0.NA <- sapply(nwl0, network.naedgecount)>0 # Update which networks have missing dyads.
  
  if(any(nwl0.NA)) stop("Transitioned-from network(s) at time(s) ", paste.and(times[-length(times)][nwl0.NA]), " has missing dyads that cannot be imputed using the selected imputation options. Fix or add imputation options via CMLE.NA.impute control parameter.")

  nwl <- lapply(seq_along(nwl)[-1], function(t){
    nwl[[t]] %n% ".PrevNets" <- list(nwl0[[t-1]])
    nwl[[t]] %n% ".TimeID" <- t
    nwl[[t]] %n% ".Time" <- times[t]
    nwl[[t]] %n% ".TimeDelta" <- times[t] - times[t-1]
    nwl[[t]]
  })

  # nwl now has all networks in the series but the first, which each network's previous network attached as a network attribute. In other words, it's a list of *transitions*.

  if(!all_identical(sapply(nwl,network.size)) || !all_identical(sapply(nwl,is.directed)) || !all_identical(sapply(nwl, `%n%`, "bipartite"))){
    stop("Networks in the network series passed must all be of the same size, directedness, and bipartite status.")
  }
  
  # Now, just combine them using the Networks() constructor.
  nw <- Networks(nwl)
  nw %n% ".combiner" <- c("NetSeries", nw %n% ".combiner")
  # Add previous networks combined.
  nw %n% ".PrevNet" <- Networks(nwl0)
  nw %ergmlhs% "constraints" <- c(ergm_flatten_conterm_list(nw %ergmlhs% "constraints") %||% term_list(list()),
                                  term_list(quote(discord(".PrevNet")), env = baseenv())) |> unique()

  structure(nw, class = c("tergm_NetSeries", class(nw)))
}

#' @rdname NetSeries
#' @description `unNetSeries()` extracts the networks in the series into a list.
#'
#' @param object a multinetwork network returned by `NetSeries()`
#'
#' @export
unNetSeries <- function(object) {
  assert_combined_network(object, "NetSeries", FALSE)
  nwl <- uncombine_network(object, populate = TRUE)
  c(nwl[[1L]] %n% ".PrevNets", nwl)
}
