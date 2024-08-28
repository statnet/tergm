#  File R/tergm.getMCMCsample.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2023 Statnet Commons
################################################################################

#' Collects a sample of networks and returns the statistics of each sample
#'
#' \code{tergm_MCMC_sample} is a low-level internal function not intended to
#' be called directly by end users. It collects a sample of networks and
#' returns the statistics of each sample, along with a toggle matrix of the
#' changes needed from the original network to each in the sample.
#'
#' This function is normally called inside \code{\link{simulate.tergm}} functions
#' to prepare inputs for the C sampling code and return its results
#'
#' @aliases tergm_MCMC_sample tergm_MCMC_slave
#' @param nw a \code{\link{network}} object
#' @param model the model, as returned by \code{\link{ergm_model}}
#' @param model.mon the optional monitoring model, as returned by \code{\link{ergm_model}}
#' @param proposal the proposal, as returned by \code{\link{ergm_proposal}}
#' @param theta the vector of curved parameters
#' @param eta the vector of natural parameters
#' @param control the list of control parameters
#' @template verbose
#' @return returns the MCMC sample as a list containing:
#' \itemize{
#'   \item statsmatrix.gen: the matrix of sampled statistics for \code{model},
#'     relative to the initial network
#'   \item statsmatrix.mon: the matrix of sampled statistics for \code{model.mon},
#'     relative to the initial network
#'   \item newnetwork: \code{ergm_state} with the final network from the
#'     sampling process
#'   \item changed: a matrix of changes, where the first column is
#'     the timestamp of the change, the second and third columns are the tail and head
#'     (respectively) of the changed dyad, and the fourth column is the edge state to which
#'     the dyad was changed; this is only returned if \code{control$changes} is \code{TRUE}
#'   \item maxchanges: the \code{maxchanges} value from the control list
#' }
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

  if(z$status)
    stop(switch(z$status,
                paste0("Number of edges in a simulated network exceeds the maximum set by the ", sQuote("MCMC.maxedges"), " control parameter."), # 1: MCMCDyn_TOO_MANY_EDGES
                "A Metropolis-Hastings proposal has failed.", # 2: MCMCDyn_MH_FAILED
                paste0("Logging of changes in the network has been requested, and the storage capacity specified by ", sQuote("MCMC.maxchanges"), " has been exceeded.") # 3: MCMCDyn_TOO_MANY_CHANGES
                )
         )

  state <- z$state

  diffedgelist<-if(control$changes) {
    if(z$diffnwtime[1]>0){
      tmp <- cbind(
        z$diffnwtime[2:(z$diffnwtime[1]+1)],
        z$diffnwtails[2:(z$diffnwtails[1]+1)],
        z$diffnwheads[2:(z$diffnwheads[1]+1)],
        z$diffnwdirs[2:(z$diffnwdirs[1]+1)]
      )
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

  if(!is.null(statsmatrix)) colnames(statsmatrix) <- param_names(model.comb, canonical = TRUE)

  # this is where we separate monitored stats from generative stats if model.mon is passed
  if(is.null(model.mon)) {
    statsmatrix.gen <- statsmatrix
    statsmatrix.mon <- NULL
  } else {
    statsmatrix.gen <- statsmatrix[,seq_len(nparam(model, canonical = TRUE)),drop=FALSE]
    statsmatrix.mon <- statsmatrix[,-seq_len(nparam(model, canonical = TRUE)),drop=FALSE]
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
  on.exit(ergm_Cstate_clear())

  collect <- if(!is.null(control$collect)) control$collect else TRUE

  maxedges <- NVL(control$MCMC.maxedges, Inf)
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
             as.integer(deInf(maxedges, "maxint")),
             as.integer(maxchanges),
             as.integer(control$changes),
             as.integer(verbose),
             PACKAGE="tergm")

  if(z$status) return(z) # If there is an error.

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
