#  File R/stergm.getMCMCsample.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2017 Statnet Commons
#######################################################################
############################################################################
# The <stergm.getMCMCsample> function collects a sample of networks and
# returns the formation and dissolution statistics of each sample, along with
# a toggle matrix of the changes needed from the original network to each
# in the sample
#
# --PARAMETERS--
#   nw             : a network object
#   model.form     : a formation model, as returned by <ergm_model>
#   model.diss     : a dissolution model, as returned by <ergm_model>
#   proposal.form: a list of parameters needed for proposals of the
#                    formations
#   proposal.diss: a list of parameters needed for proposals of the
#                    dissolutions
#   init.form    : the initial model coefficients for 'model.form'
#   init.diss    : the initial model coefficients for 'model.diss'
#   control     : the list of parameters controlling the MCMC algorithm;
#                    recognized components include:
#      maxchanges    :  this ends up defining the "MCMCDyn workspace", but
#                       in a peculiar way: the 'maxchanges' you give, let's
#                       call this number 'mc' is fed into this equation:
#                           mc < 5^x * mc - 55
#                       the 'x' that solves the equation is then used and
#                       and 'maxchanges' and the MCMCDyn workspace is defined as
#                           5^(x-1) * (the original mc)
#                       if the computed Clist.form$nedges is > than the original
#                       'maxchanges', 'maxchanges' is ignored and this odd little
#                       process occurs to the 'nedeges'
#      parallel      :  the number of threads in which to run sampling
#      samplesize    :  the desired MCMC sample size
#      changes       :  whether to return the toggle matrix; non-NULL values
#                       will include 'changed' in the return list, NULL
#                       will not
#      MH.burnin     :  this is accepted as 'MH_interval' which determines
#                       the number of proposals in each MCMC step
#      time.burnin   :  the number of MCMC steps to disregard for the burnin
#                       period 
#      time.interval :  the number of MCMC steps to disregard in between
#                       sampled networks
#   verbose        : whether this and other functions should be verbose
#
# Note:  In reality, there should be many fewer arguments to this function,
# since most info should be passed via Clist (this is, after all, what Clist
# is for:  Holding all arguments required for the .C call).  In particular,
# the elements of proposal.form, control, verbose, should certainly
# be part of Clist.  But this is a project for another day!
#
# --RETURNED--
#   the MCMC sample as a list containing:
#     statsmatrix.form: the matrix of sampled statistics for 'model.form' RELATIVE TO INITIAL NETWORK
#     statsmatrix.diss: the matrix of sampled statistics for 'model.form' RELATIVE TO INITIAL NETWORK
#     statsmatrix.mon: the matrix of sampled statistics for 'model.mon' RELATIVE TO INITIAL NETWORK
#     newnetwork      : the final network from the sampling process
#     changed         : a toggle matrix, where the first column is
#                       the timestamp of the toggle and the 2nd and 3rd
#                       columns are the head & tail of the toggle; this
#                       is only returned if 'control'$changes is not NULL
#     maxchanges      : the "MCMC Dyn workspace"; see 'maxchanges' in the
#                       input param list
#
############################################################################





#' Collects a sample of networks and returns the formation and dissolution
#' statistics of each sample
#' 
#' \code{stergm.getMCMCsample} is a low-level internal function not intended to
#' be called directly by end users. It collects a sample of networks and
#' returns the formation and dissolution statistics of each sample, along with
#' a toggle matrix of the changes needed from the original network to each in
#' the sample.
#' 
#' This function is normally called inside \code{\link{simulate.stergm}} to
#' prepare inputs for the C sampling code and return its results
#' 
#' @aliases stergm.getMCMCsample stergm.getMCMCsample.slave
#' @param nw a \code{\link{network}} object
#' @param model.form,model.diss,model.mon formation, dissolution, and
#' monitoring model, as returned by \code{\link{ergm_model}}
#' @param proposal.form,proposal.diss a list of parameters needed for
#' proposals of the formations and dissolutions
#' @param eta.form,eta.diss vectors of natural parameters.
#' @param control list of control paramters, probably from
#' \code{\link{control.stergm}}
#' @param verbose logical; whether this and other functions should be verbose
#' @return returns the MCMC sample as a list containing: \itemize{
#' \item statsmatrix.form: the matrix of sampled statistics for 'model.form'
#' RELATIVE TO INITIAL NETWORK \item statsmatrix.diss: the matrix of sampled
#' statistics for 'model.form' RELATIVE TO INITIAL NETWORK
#' \item statsmatrix.mon: the matrix of sampled statistics for 'model.mon'
#' RELATIVE TO INITIAL NETWORK \item newnetwork : the final network from the
#' sampling process \item changed : a toggle matrix, where the first column is
#' the timestamp of the toggle and the 2nd and 3rd columns are the head & tail
#' of the toggle; this is only returned if `control$changes` is not NULL
#' \item maxchanges : the "MCMC Dyn workspace"; see 'maxchanges' in the input
#' param list }
#' @seealso \code{\link{simulate.stergm}}
#' @keywords internal
#' @export stergm.getMCMCsample
stergm.getMCMCsample <- function(nw, model.form, model.diss, model.mon,
                                  proposal.form, proposal.diss, eta.form, eta.diss, control, 
                                  verbose){


  #
  #   Check for truncation of the returned edge list
  #
  Clist.form <- ergm.Cprepare(nw, model.form)
  Clist.diss <- ergm.Cprepare(nw, model.diss)
  Clist.mon <- if(!is.null(model.mon)) ergm.Cprepare(nw, model.mon) else NULL
  
  z <- stergm.getMCMCsample.slave(Clist.form, Clist.diss, Clist.mon, proposal.form, proposal.diss, eta.form, eta.diss, control, verbose)

  newnetwork<-newnw.extract(nw,z)
  if(is.durational(model.form) || is.durational(model.diss) || is.durational(model.mon)){
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

  statsmatrix.form <- z$statsmatrix.form
  statsmatrix.diss <- z$statsmatrix.diss
  statsmatrix.mon <- z$statsmatrix.mon
  
  if(!is.null(statsmatrix.form)) colnames(statsmatrix.form) <- model.form$coef.names
  if(!is.null(statsmatrix.diss)) colnames(statsmatrix.diss) <- model.diss$coef.names
  if(!is.null(model.mon)) colnames(statsmatrix.mon) <- model.mon$coef.names

  list(statsmatrix.form=statsmatrix.form, statsmatrix.diss=statsmatrix.diss, statsmatrix.mon=statsmatrix.mon,
       newnetwork=newnetwork,
       changed=diffedgelist,
       maxchanges=control$MCMC.maxchanges)
}

#' @rdname stergm.getMCMCsample
#' @description \code{stergm.getMCMCsample.slave} is an even
#'   lower-level function that actually calls the C code.
#' @param Clist.form,Clist.diss,Clist.mon formation, dissolution, and
#'   monitoring "Clist", as returned by \code{\link{ergm.Cprepare}}
#' @useDynLib tergm
#' @export stergm.getMCMCsample.slave
stergm.getMCMCsample.slave <- function(Clist.form, Clist.diss, Clist.mon, proposal.form, proposal.diss, eta.form, eta.diss, control, verbose){
  collect.form<-if(!is.null(control$collect.form)) control$collect.form else TRUE
  collect.diss<-if(!is.null(control$collect.diss)) control$collect.diss else TRUE
  maxedges <- control$MCMC.init.maxedges
  maxchanges <- control$MCMC.init.maxchanges
  repeat{
    #FIXME: Separate MCMC control parameters and properly attach them.
    
    z <- .C("MCMCDyn_wrapper",
            # Observed network.
            as.integer(Clist.form$tails), as.integer(Clist.form$heads),
            time = if(is.null(Clist.form$time)) as.integer(0) else as.integer(Clist.form$time),
            lasttoggle = as.integer(NVL(Clist.form$lasttoggle,Clist.diss$lasttoggle,Clist.mon$lasttoggle,0)),  
            as.integer(Clist.form$nedges),
            as.integer(Clist.form$n),
            as.integer(Clist.form$dir), as.integer(Clist.form$bipartite),
            # Formation terms and proposals.
            as.integer(Clist.form$nterms), 
            as.character(Clist.form$fnamestring),
            as.character(Clist.form$snamestring),
            as.character(proposal.form$name), as.character(proposal.form$pkgname),
            as.double(Clist.form$inputs), as.double(ergm:::.deinf(eta.form)),
            # Dissolution terms and proposals.
            as.integer(Clist.diss$nterms), 
            as.character(Clist.diss$fnamestring),
            as.character(Clist.diss$snamestring),
            as.character(proposal.diss$name), as.character(proposal.diss$pkgname),
            as.double(Clist.diss$inputs), as.double(ergm:::.deinf(eta.diss)),
            # Monitored terms.
            if(!is.null(Clist.mon)) as.integer(Clist.mon$nterms) else as.integer(0), 
            if(!is.null(Clist.mon)) as.character(Clist.mon$fnamestring) else character(0),
            if(!is.null(Clist.mon)) as.character(Clist.mon$snamestring) else character(0),
            if(!is.null(Clist.mon)) as.double(Clist.mon$inputs) else double(0),
            # Degree bounds.
            as.integer(proposal.form$arguments$constraints$bd$attribs), 
            as.integer(proposal.form$arguments$constraints$bd$maxout), as.integer(proposal.form$arguments$constraints$bd$maxin),
            as.integer(proposal.form$arguments$constraints$bd$minout), as.integer(proposal.form$arguments$constraints$bd$minin),
            as.integer(proposal.form$arguments$constraints$bd$condAllDegExact), as.integer(length(proposal.form$arguments$constraints$bd$attribs)),
            # MCMC settings.
            as.integer(control$time.samplesize), as.integer(control$MCMC.burnin.min), as.integer(control$MCMC.burnin.max), as.double(control$MCMC.burnin.pval), as.double(control$MCMC.burnin.add),
            as.integer(control$time.burnin), as.integer(control$time.interval),
            # Space for output.
            collect.form = as.integer(collect.form), s.form = if(collect.form) double(Clist.form$nstats*(control$time.samplesize+1)) else double(0),
            collect.diss = as.integer(collect.diss), s.diss = if(collect.diss) double(Clist.diss$nstats*(control$time.samplesize+1)) else double(0),
            s.mon = if(!is.null(Clist.mon)) double(Clist.mon$nstats*(control$time.samplesize+1)) else double(0),
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
  
  statsmatrix.form <-
    if(collect.form) matrix(z$s.form, nrow=control$time.samplesize+1,
                            ncol=Clist.form$nstats,
                            byrow = TRUE)[-1,,drop=FALSE]
    else
      NULL
  
  statsmatrix.diss <-
    if(collect.diss) matrix(z$s.diss, nrow=control$time.samplesize+1,
                            ncol=Clist.diss$nstats,
                            byrow = TRUE)[-1,,drop=FALSE]
  else
    NULL

  statsmatrix.mon <-
    if(!is.null(Clist.mon))
      matrix(z$s.mon, nrow=control$time.samplesize+1,
             ncol=Clist.mon$nstats,
             byrow = TRUE)[-1,,drop=FALSE]
    else
      NULL

  # Blank all elements of z that we don't want to bother returning.
  zn <- names(z)
  for(i in rev(seq_along(zn))){ # Do in reverse, to preserve indexing.
    if(! zn[i] %in% c("time", "lasttoggle", "newnwtails", "newnwheads", "diffnwtime", "diffnwtails", "diffnwheads", "diffnwdirs", "status"))
      z[[i]] <- NULL
  }  
  c(z,
    list(statsmatrix.form = statsmatrix.form,
         statsmatrix.diss = statsmatrix.diss,
         statsmatrix.mon = statsmatrix.mon))
}
