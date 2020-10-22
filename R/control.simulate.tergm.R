#  File R/control.simulate.tergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################

#' Auxiliary for Controlling Temporal ERGM Simulation
#' 
#' Auxiliary function as user interface for fine-tuning TERGM simulation.
#' 
#' This function is only used within a call to the \code{\link{simulate}}
#' function.  See the \code{usage} section in \code{\link{simulate.tergm}} for
#' details.
#' 
#' @param
#'   MCMC.burnin.min,MCMC.burnin.max,MCMC.burnin.pval,MCMC.burnin.add
#'   Number of Metropolis-Hastings steps 
#'   per time step used in simulation. By default, this
#'   is determined adaptively by keeping track of increments in the
#'   Hamming distance between the transitioned-from network and the
#'   network being sampled. Once \code{MCMC.burnin.min} steps have elapsed, the
#'   increments are tested against 0, and when their average number
#'   becomes statistically indistinguishable from 0 (with the p-value
#'   being greater than \code{MCMC.burnin.pval}), or
#'   \code{MCMC.burnin.max} steps are proposed, whichever comes first,
#'   the simulation is stopped after an additional
#'   \code{MCMC.burnin.add} times the number of elapsed steps had been
#'   taken.  (Stopping immediately would bias the sampling.)
#' 
#'   To use a fixed number of steps, set both \code{MCMC.burnin.min}
#'   and \code{MCMC.burnin.max} to the desired number of steps.
#' @param MCMC.prop.weights Specifies the
#'   proposal distribution used in the MCMC Metropolis-Hastings
#'   algorithm. Possible
#'   choices are \code{"TNT"} or \code{"random"}; the
#'   \code{"default"}.  The \code{TNT} (tie / no tie) option puts
#'   roughly equal weight on selecting a dyad with or without a tie as
#'   a candidate for toggling, whereas the \code{random} option puts
#'   equal weight on all possible dyads, though the interpretation of
#'   \code{random} may change according to the constraints in place.
#'   When no constraints are in place, the default is TNT, which
#'   appears to improve Markov chain mixing particularly for networks
#'   with a low edge density, as is typical of many realistic social
#'   networks.
#' @param MCMC.prop.args An alternative,
#'   direct way of specifying additional arguments to the proposal.
#' @param MCMC.maxchanges Maximum number of toggles changes for
#'   which to allocate space.
#' @param MCMC.packagenames Names of packages in which to look for
#'   change statistic functions in addition to those
#'   autodetected. This argument should not be needed outside of very
#'   strange setups.
#' @param term.options A list of additional arguments to be passed to term initializers. It can also be set globally via `option(ergm.term=list(...))`.
#' @param MCMC.maxedges Maximum number of edges expected in
#'   network.
#' @param MCMC.burnin,MCMC.burnin.mul No longer used. See
#'   \code{MCMC.burnin.min}, \code{MCMC.burnin.max},
#'   \code{MCMC.burnin.pval}, \code{MCMC.burnin.pval}, and
#'   \code{MCMC.burnin.add}.
#' @return A list with arguments as components.
#' @seealso \code{\link{simulate.tergm}},
#'   \code{\link{simulate.formula}}.  \code{\link{control.tergm}}
#'   performs a similar function for \code{\link{tergm}}.
#' @keywords models
#' @export control.simulate.tergm
control.simulate.tergm<-function(MCMC.burnin.min=NULL,
                                  MCMC.burnin.max=NULL,
                                  MCMC.burnin.pval=NULL,
                                  MCMC.burnin.add=NULL,
                                  MCMC.prop.weights=NULL,MCMC.prop.args=NULL,
                                  MCMC.maxedges=NULL,
                                  MCMC.packagenames=NULL,

                                  term.options=NULL,
                                  
                                  MCMC.maxchanges=NULL){
    control<-list()
    for(arg in names(formals(sys.function())))
      control[arg]<-list(get(arg))

    set.control.class("control.simulate.tergm")
  }


#' @rdname control.simulate.tergm
#' @export control.simulate.network.tergm
control.simulate.network.tergm<-function(MCMC.burnin.min=1000,
                                   MCMC.burnin.max=100000,
                                   MCMC.burnin.pval=0.5,
                                   MCMC.burnin.add=1,
                                   MCMC.burnin=NULL, MCMC.burnin.mul=NULL,
                                   MCMC.prop.weights="default",MCMC.prop.args=NULL,
                                   MCMC.maxedges=Inf,
                                   MCMC.packagenames=c(),

                                   term.options=NULL,
                                   
                                   MCMC.maxchanges=1000000){
    control<-list()
    for(arg in names(formals(sys.function())))
      control[arg]<-list(get(arg))

    set.control.class("control.simulate.network.tergm")
  }
