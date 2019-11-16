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
  eta.comb <- c(eta, rep(0, NVL(model.mon$etamap$etalength, 0)))

  # always collect if monitoring model is passed
  control$collect <- NVL(control$collect, TRUE) || !is.null(model.mon)

  #
  #   Check for truncation of the returned edge list
  #  
  Clist <- ergm.Cprepare(nw, model.comb)
  
  z <- tergm_MCMC_slave(Clist, proposal, eta.comb, control, verbose)

  newnetwork<-as.network(pending_update_network(nw, z))
  if(is.durational(model.comb)){
    newnetwork %n% "time" <- z$time
    newnetwork %n% "lasttoggle" <- z$lasttoggle
  }
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

  # this is where we separate monitored stats from fd stats if model.mon is passed
  if(is.null(model.mon)) {
    statsmatrix.fd <- statsmatrix
    statsmatrix.mon <- NULL
  } else {
    statsmatrix.fd <- statsmatrix[,1:(NCOL(statsmatrix) - model.mon$etamap$etalength),drop=FALSE]
    statsmatrix.mon <- statsmatrix[,(NCOL(statsmatrix) - model.mon$etamap$etalength + 1):NCOL(statsmatrix),drop=FALSE]
  }

  list(statsmatrix.fd=statsmatrix.fd,
       statsmatrix.mon=statsmatrix.mon,
       newnetwork=newnetwork,
       changed=diffedgelist,
       maxchanges=control$MCMC.maxchanges)
}

#' @rdname tergm_MCMC_sample
#' @description \code{tergm_MCMC_slave} is an even
#'   lower-level function that actually calls the C code.
#' @param Clist the "Clist" for the network and model under consideration, 
#'   as returned by \code{\link{ergm.Cprepare}}
#' @useDynLib tergm
#' @export
tergm_MCMC_slave <- function(Clist, proposal, eta, control, verbose){
  collect <- if(!is.null(control$collect)) control$collect else TRUE

  maxedges <- control$MCMC.init.maxedges
  maxchanges <- control$MCMC.init.maxchanges
  repeat{
    #FIXME: Separate MCMC control parameters and properly attach them.
    
    # lasttoggle must hold the inputs while also having room for the outputs
    lasttoggle <- c(NROW(Clist$lasttoggle), Clist$lasttoggle)    
    lasttoggle <- c(lasttoggle, rep(0, 3*maxedges + 1 - length(lasttoggle)))
    
    z <- .C("MCMCDyn_wrapper",
            # Observed network.
            as.integer(Clist$tails), as.integer(Clist$heads),
            time = if(is.null(Clist$time)) as.integer(0) else as.integer(Clist$time),
            lasttoggle = as.integer(lasttoggle),  
            as.integer(Clist$nedges),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            # terms and proposals.
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring),
            as.character(proposal$name), as.character(proposal$pkgname),
            as.double(Clist$inputs), as.double(ergm:::.deinf(eta)),
            # Degree bounds.
            as.integer(proposal$arguments$constraints$bd$attribs), 
            as.integer(proposal$arguments$constraints$bd$maxout), as.integer(proposal$arguments$constraints$bd$maxin),
            as.integer(proposal$arguments$constraints$bd$minout), as.integer(proposal$arguments$constraints$bd$minin),
            as.integer(proposal$arguments$constraints$bd$condAllDegExact), as.integer(length(proposal$arguments$constraints$bd$attribs)),
            # MCMC settings.
            as.integer(control$time.samplesize), as.integer(control$MCMC.burnin.min), as.integer(control$MCMC.burnin.max), as.double(control$MCMC.burnin.pval), as.double(control$MCMC.burnin.add),
            as.integer(control$time.burnin), as.integer(control$time.interval),
            # Space for output.
            collect = as.integer(collect), s = if(collect) double(Clist$nstats*(control$time.samplesize+1)) else double(0),
            as.integer(maxedges),
            newnwtails = integer(maxchanges), newnwheads = integer(maxchanges), 
            as.integer(maxchanges),
            as.integer(control$changes),
            diffnwtime = if(control$changes) integer(maxchanges) else integer(0),
            diffnwtails = if(control$changes) integer(maxchanges) else integer(0),
            diffnwheads = if(control$changes) integer(maxchanges) else integer(0),
            diffnwdirs = if(control$changes) integer(maxchanges) else integer(0),
            as.integer(verbose),
            status = integer(1), # 0 = OK, MCMCDyn_TOO_MANY_EDGES = 1, MCMCDyn_MH_FAILED = 2, MCMCDyn_TOO_MANY_CHANGES = 3
            PACKAGE="tergm")

    if(z$status==0) break;
    if(z$status==1){
      maxedges <- 5*maxedges
      if(verbose>0) message("Too many edges encountered in the simulation. Increasing capacity to ", maxedges)
    }
    if(z$status==3){
      maxchanges <- 5*maxchanges
      if(verbose>0) message("Too many changes elapsed in the simulation. Increasing capacity to ", maxchanges)
    }
  }
  
  statsmatrix <-
    if(collect) matrix(z$s, nrow=control$time.samplesize+1,
                            ncol=Clist$nstats,
                            byrow = TRUE)[-1,,drop=FALSE]
    else
      NULL
  
  # Blank all elements of z that we don't want to bother returning.
  zn <- names(z)
  for(i in rev(seq_along(zn))){ # Do in reverse, to preserve indexing.
    if(! zn[i] %in% c("time", "lasttoggle", "newnwtails", "newnwheads", "diffnwtime", "diffnwtails", "diffnwheads", "diffnwdirs", "status"))
      z[[i]] <- NULL
  }

  # subselect the portion of z$lasttoggle that corresponds to actual edges and not just buffer
  z$lasttoggle <- matrix(z$lasttoggle[2:(3*z$lasttoggle[1] + 1)],nrow=z$lasttoggle[1])
  
  c(z,
    list(statsmatrix = statsmatrix))
}
