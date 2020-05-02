#  File R/tergm.getMCMCsample.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2019 Statnet Commons
#######################################################################

#' Collects a sample of networks and returns the statistics of each sample
#' 
#' \code{tergm_MCMC_sample} is a low-level internal function not intended to
#' be called directly by end users. It collects a sample of networks and
#' returns the statistics of each sample, along with
#' a toggle matrix of the changes needed from the original network to each in
#' the sample.
#' 
#' This function is normally called inside \code{\link{simulate.tergm}} to
#' prepare inputs for the C sampling code and return its results
#' 
#' @aliases tergm_MCMC_sample tergm_MCMC_slave
#' @param nw a \code{\link{network}} object
#' @param model the model, as returned by \code{\link{ergm_model}}
#' @param model.mon the optional monitoring model, as returned by \code{\link{ergm_model}}
#' @param proposal a list of parameters needed for
#' proposals
#' @param eta vector of natural parameters.
#' @param control list of control paramters, probably from
#' \code{\link{control.tergm}}
#' @param verbose logical; whether this and other functions should be verbose
#' @return returns the MCMC sample as a list containing: \itemize{
#' \item statsmatrix: the matrix of sampled statistics for 'model'
#' RELATIVE TO INITIAL NETWORK \item newnetwork : the final network from the
#' sampling process \item changed : a toggle matrix, where the first column is
#' the timestamp of the toggle and the 2nd and 3rd columns are the head & tail
#' of the toggle; this is only returned if `control$changes` is not NULL
#' \item maxchanges : the "MCMC Dyn workspace"; see 'maxchanges' in the input
#' param list }
#' @seealso \code{\link{simulate.tergm}}
#' @keywords internal
#' @export
tergm_MCMC_sample <- function(nw, model, model.mon = NULL,
                               proposal, control,
                               theta,
                               verbose=FALSE,...,
                               eta = ergm.eta(theta, model$etamap)
                               ){
  # this is where we combine models and pad out eta 
  # with 0s as necessary to accomodate the monitoring model
  model.comb <- c(model, model.mon)
  proposal$aux.slots <- model.comb$slots.extra.aux$proposal
  
  eta.comb <- c(eta, rep(0, NVL(model.mon$etamap$etalength, 0)))

  # always collect if monitoring model is passed
  control$collect <- NVL(control$collect, TRUE) || !is.null(model.mon)

  #
  #   Check for truncation of the returned edge list
  #  
  
  state <- ergm_state(nw, model=model.comb, proposal=proposal, stats=rep(0,nparam(model.comb, canonical=TRUE)))
  
  z <- tergm_MCMC_slave(state, eta.comb, control, verbose)

  # check that this is the correct replacement for as.network(pending_update_network(...))
  state <- z$state
  
  diffedgelist<-if(control$changes) {
    if(z$diffnwtime[1]>0){
      tmp <- cbind(z$diffnwtime[2:(z$diffnwtime[1]+1)],z$diffnwtails[2:(z$diffnwtails[1]+1)],z$diffnwheads[2:(z$diffnwheads[1]+1)],z$diffnwdirs[2:(z$diffnwdirs[1]+1)])
      colnames(tmp) <- c("time","tail","head","to")
      tmp
    }else{
      tmp <- matrix(0, ncol=4, nrow=0)
      colnames(tmp) <- c("time","tail","head","to")
      tmp
    }
  }else{
    NULL
  }
  mode(diffedgelist) <- "integer" # Might save some memory.

  statsmatrix <- z$statsmatrix
  
  if(!is.null(statsmatrix)) colnames(statsmatrix) <- model.comb$coef.names

  # this is where we separate monitored stats from generative stats if model.mon is passed
  if(is.null(model.mon)) {
    statsmatrix.gen <- statsmatrix
    statsmatrix.mon <- NULL
  } else {
    statsmatrix.gen <- statsmatrix[,1:(NCOL(statsmatrix) - model.mon$etamap$etalength),drop=FALSE]
    statsmatrix.mon <- statsmatrix[,(NCOL(statsmatrix) - model.mon$etamap$etalength + 1):NCOL(statsmatrix),drop=FALSE]
  }

  list(statsmatrix.gen=statsmatrix.gen,
       statsmatrix.mon=statsmatrix.mon,
       newnetwork=state,
       changed=diffedgelist,
       maxchanges=control$MCMC.maxchanges)
}

#' @rdname tergm_MCMC_sample
#' @description \code{tergm_MCMC_slave} is an even
#'   lower-level function that actually calls the C code.
#' @useDynLib tergm
#' @export
tergm_MCMC_slave <- function(state, eta, control, verbose){
  collect <- if(!is.null(control$collect)) control$collect else TRUE

  ## should be maxedges, init.maxchanges, and max.maxchanges
  ## add copying logic to MCMCDyn1step_advance, as well as addtnl arg
  maxedges <- control$MCMC.maxedges
  maxchanges <- control$MCMC.maxchanges
  
  z <- .Call("MCMCDyn_wrapper",
             state,
             as.double(deInf(eta)),
             # MCMC settings.
             as.integer(control$time.samplesize),
             as.integer(control$MCMC.burnin.min),
             as.integer(control$MCMC.burnin.max),
             as.double(control$MCMC.burnin.pval),
             as.double(control$MCMC.burnin.add),
             as.integer(control$time.burnin),
             as.integer(control$time.interval),
             # output settings.
             as.integer(collect),
             as.integer(maxedges),
             as.integer(maxchanges),
             as.integer(control$changes),
             as.integer(verbose),
             PACKAGE="tergm")

  if(z$status != 0) stop("MCMCDyn failed with error code ", z$status)
  
  z$state <- update(z$state)
  
  statsmatrix <-
    if(collect) matrix(z$s, nrow=control$time.samplesize+1,
                            ncol=nparam(state,canonical=TRUE),
                            byrow = TRUE)[-1,,drop=FALSE]
    else
      NULL
  
  c(z,
    list(statsmatrix = statsmatrix))
}
