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
#'   \code{MCMC.burnin.add} times the number of elapsed steps have been
#'   taken.  (Stopping immediately would bias the sampling.)
#' 
#'   To use a fixed number of steps, set \code{MCMC.burnin.min}
#'   and \code{MCMC.burnin.max} to the same value.
#'
#' @param MCMC.prop.weights Specifies the proposal weighting scheme to 
#'   be used in the MCMC Metropolis-Hastings algorithm.  Possible
#'   choices may be determined by calling \code{\link{ergm_proposal_table}}.
#'
#' @param MCMC.prop.args An alternative, direct way of specifying 
#'   additional arguments to the proposal.
#'
#' @param MCMC.prop Hints and/or constraints for selecting and initializing the proposal.
#'
#' @param MCMC.maxchanges Maximum number of changes for
#'   which to allocate space.
#'
#' @param MCMC.packagenames Names of packages in which to look for
#'   change statistic functions in addition to those
#'   autodetected. This argument should not be needed outside of very
#'   strange setups.
#'
#' @param term.options A list of additional arguments to be passed to term initializers. 
#'   It can also be set globally via \code{options(ergm.term=list(...))}.
#'
#' @param MCMC.maxedges Maximum number of edges expected in
#'   network.
#'
#' @return A list with arguments as components.
#'
#' @seealso \code{\link{simulate.tergm}},
#'   \code{\link{simulate.formula}}.  \code{\link{control.tergm}}
#'   performs a similar function for \code{\link{tergm}}.
#'
#' @keywords models
#'
#' @export control.simulate.tergm
control.simulate.tergm <- function(MCMC.burnin.min = NULL,
                                   MCMC.burnin.max = NULL,
                                   MCMC.burnin.pval = NULL,
                                   MCMC.burnin.add = NULL,
                                   
                                   MCMC.prop = NULL,
                                   MCMC.prop.weights = NULL,
                                   MCMC.prop.args = NULL,
                                   
                                   MCMC.maxedges = NULL,
                                   MCMC.maxchanges = NULL,
                                   
                                   term.options = NULL,
                                   
                                   MCMC.packagenames = NULL) {

  control <- list()
  for(arg in names(formals(sys.function())))
    control[arg] <- list(get(arg))

  set.control.class("control.simulate.tergm")
}


#' @rdname control.simulate.tergm
#' @export control.simulate.formula.tergm
control.simulate.formula.tergm <- function(MCMC.burnin.min = 1000,
                                           MCMC.burnin.max = 100000,
                                           MCMC.burnin.pval = 0.5,
                                           MCMC.burnin.add = 1,
                                           
                                           MCMC.prop = ~discord + sparse,
                                           MCMC.prop.weights = "default",
                                           MCMC.prop.args = NULL,
                                           
                                           MCMC.maxedges = Inf,
                                           MCMC.maxchanges = 1000000,
                                           
                                           term.options = NULL,
                                           
                                           MCMC.packagenames = c()) {
                                           
  control <- list()
  for(arg in names(formals(sys.function())))
    control[arg] <- list(get(arg))

  set.control.class("control.simulate.formula.tergm")
}
