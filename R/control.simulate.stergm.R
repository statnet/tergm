#  File R/control.simulate.stergm.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2022 Statnet Commons
################################################################################
########################################################################
# The <control.simulate.X> functions each create a list of paramaters
# for customizing simulation rountines
#
# --PARAMETERS--
#   prop.weights.form/
#   ptop.weighs.diss: specifies the method used to allocate probabilities
#                     of being proposed to dyads for the formation/dis-
#                     solution processes; options are "TNT",
#                     "random", and "default"; default="default" which
#                      picks a reasonable default considering any constraints
#   final           : whether only the final network of the simulation
#                     process should be returned (T or F); default= FALSE
#                     in which case, models, coefficients, stats matrices,
#                     and the toggle matrix are returned
#   maxchanges      : the maximum number of changes for which to allocate
#                     space; default=1000000
#
# --IGNORED--
#   prop.args.form/
#   prop.args.diss: an alternative, direct way of specifying additional
#                   arguments to proposal; as far as I can tell, the
#                   only use for 'prop.args' is to supply the name
#                   of a nodal attribute for use in the
#                   <InitErgmProposal.nobetweengroupties> function, but this
#                   function will never be called in the path from
#                   <simulate.stergm> which is the only code using this
#                   control list.
#
# --RETURNED--
#   a list of the above parameters
#
#########################################################################

#' @rdname control.simulate.stergm
#' @export control.simulate.network
control.simulate.network <- function(MCMC.burnin.min = 1000,
                                     MCMC.burnin.max = 100000,
                                     MCMC.burnin.pval = 0.5,
                                     MCMC.burnin.add = 1,
                                     
                                     MCMC.prop.form = ~discord + sparse,
                                     MCMC.prop.diss = ~discord + sparse,
                                     
                                     MCMC.prop.weights.form = "default",
                                     MCMC.prop.weights.diss = "default",
                                                                        
                                     MCMC.prop.args.form = NULL,
                                     MCMC.prop.args.diss = NULL,                                  
                                   
                                     MCMC.maxedges = Inf,
                                     MCMC.maxchanges = 1000000,

                                     term.options = NULL,
                                     
                                     MCMC.packagenames = c()) {
                                     
  control <- list()
  for(arg in names(formals(sys.function())))
    control[arg] <- list(get(arg))

  set.control.class("control.simulate.network")
}


#' Auxiliary for Controlling Separable Temporal ERGM Simulation
#' 
#' Auxiliary function as user interface for fine-tuning STERGM simulation.
#' 
#' This function is only used within a call to the \code{\link{simulate}}
#' function.  See the \code{usage} section in \code{\link{simulate.stergm}} for
#' details.
#'
#' These functions are included for backwards compatibility, and users are  
#' encouraged to use \code{control.simulate.tergm} or
#' \code{control.simulate.formula.tergm} with the \code{\link{simulate.tergm}}
#' family of functions instead.  When a
#' \code{control.simulate.stergm} or \code{control.simulate.network} object 
#' is passed to one of the \code{\link{simulate.stergm}} functions, 
#' the corresponding \code{\link{simulate.tergm}} function is invoked,
#' and uses the formation proposal control arguments, ignoring the 
#' dissolution proposal control arguments.
#' 
#' Note:  The old \code{dissolution} formula in \code{stergm} represents
#' tie persistence.  As a result it maps to the new \code{Persist()} operator
#' in \code{tergm}, NOT the \code{Diss()} operator
#' 
#' @param
#'   MCMC.burnin.min,MCMC.burnin.max,MCMC.burnin.pval,MCMC.burnin.add
#'   Number of Metropolis-Hastings steps per time step used in simulation. By default, this
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
#' @param MCMC.prop.weights.form Specifies the proposal weighting scheme to 
#'   be used in the MCMC Metropolis-Hastings algorithm.  Possible
#'   choices may be determined by calling \code{\link{ergm_proposal_table}}.
#' 
#' @param MCMC.prop.args.form An alternative,
#'   direct way of specifying additional arguments to proposals.
#' 
#' @param MCMC.prop.form Hints and/or constraints for selecting and initializing the proposal.
#' 
#' @param MCMC.prop.weights.diss,MCMC.prop.args.diss,MCMC.prop.diss Ignored. These are included
#'        for backwards compatibility of calls to \code{control}
#'        functions only; they have no effect on \code{simulate} behavior.
#'
#' @param MCMC.maxchanges Maximum number of changes for
#'   which to allocate space.
#' 
#' @template control_MCMC_packagenames
#'
#' @template term_options
#' 
#' @template control_MCMC_maxedges
#' 
#' @return A list with arguments as components.
#' 
#' @seealso \code{\link{simulate.stergm}},
#'   \code{\link{simulate.formula}}.  \code{\link{control.stergm}}
#'   performs a similar function for \code{\link{stergm}}.
#'
#' @keywords models
#'
#' @export control.simulate.stergm
control.simulate.stergm <- function(MCMC.burnin.min = NULL,
                                    MCMC.burnin.max = NULL,
                                    MCMC.burnin.pval = NULL,
                                    MCMC.burnin.add = NULL,
                                  
                                    MCMC.prop.form = NULL,
                                    MCMC.prop.diss = NULL,                                  
                                    
                                    MCMC.prop.weights.form = NULL,
                                    MCMC.prop.weights.diss = NULL,

                                    MCMC.prop.args.form = NULL,                                    
                                    MCMC.prop.args.diss = NULL,                                  
                                    
                                    MCMC.maxedges = NULL,
                                    MCMC.maxchanges = NULL,
                                    
                                    term.options = NULL,
                                    
                                    MCMC.packagenames = NULL) {

  control <- list()
  for(arg in names(formals(sys.function())))
    control[arg] <- list(get(arg))

  set.control.class("control.simulate.stergm")
}
