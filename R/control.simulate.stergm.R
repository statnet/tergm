#  File R/control.simulate.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2018 Statnet Commons
#######################################################################
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
control.simulate.network<-function(MCMC.burnin.min=1000,
                                   MCMC.burnin.max=100000,
                                   MCMC.burnin.pval=0.5,
                                   MCMC.burnin.add=1,
                                   MCMC.burnin=NULL, MCMC.burnin.mul=NULL,
                                   MCMC.prop.weights.form="default",MCMC.prop.args.form=NULL,
                                   MCMC.prop.weights.diss="default",MCMC.prop.args.diss=NULL,                                  
                                   MCMC.init.maxedges=20000,
                                   MCMC.packagenames=c(),
                                   
                                   MCMC.init.maxchanges=1000000){
    if(!is.null(MCMC.burnin) || !is.null(MCMC.burnin.mul)) stop("Control parameters MCMC.burnin and MCMC.burnin.mul are no longer used. See help for EGMME.MCMC.burnin.min, EGMME.MCMC.burnin.max, EGMME.MCMC.burnin.pval, EGMME.MCMC.burnin.pval, and CMLE.MCMC.burnin and CMLE.MCMC.interval for their replacements.")
    control<-list()
    for(arg in names(formals(sys.function())))
      control[arg]<-list(get(arg))

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
#' @param
#'   MCMC.burnin.min,MCMC.burnin.max,MCMC.burnin.pval,MCMC.burnin.add
#'   Number of Metropolis-Hastings steps per phase (formation and
#'   dissolution) per time step used in simulation. By default, this
#'   is determined adaptively by keeping track of increments in the
#'   Hamming distance between the transitioned-from network and the
#'   network being sampled (formation network or dissolution
#'   network). Once \code{MCMC.burnin.min} steps have elapsed, the
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
#' @param MCMC.prop.weights.form,MCMC.prop.weights.diss Specifies the
#'   proposal distribution used in the MCMC Metropolis-Hastings
#'   algorithm for formation and dissolution, respectively. Possible
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
#' @param MCMC.prop.args.form,MCMC.prop.args.diss An alternative,
#'   direct way of specifying additional arguments to proposals.
#' @param MCMC.init.maxchanges Maximum number of toggles changes for
#'   which to allocate space.
#' @param MCMC.packagenames Names of packages in which to look for
#'   change statistic functions in addition to those
#'   autodetected. This argument should not be needed outside of very
#'   strange setups.
#' @param MCMC.init.maxedges Maximum number of edges expected in
#'   network.
#' @param MCMC.burnin,MCMC.burnin.mul No longer used. See
#'   \code{MCMC.burnin.min}, \code{MCMC.burnin.max},
#'   \code{MCMC.burnin.pval}, \code{MCMC.burnin.pval}, and
#'   \code{MCMC.burnin.add}.
#' @return A list with arguments as components.
#' @seealso \code{\link{simulate.stergm}},
#'   \code{\link{simulate.formula}}.  \code{\link{control.stergm}}
#'   performs a similar function for \code{\link{stergm}}.
#' @keywords models
#' @export control.simulate.stergm
control.simulate.stergm<-function(MCMC.burnin.min=NULL,
                                  MCMC.burnin.max=NULL,
                                  MCMC.burnin.pval=NULL,
                                  MCMC.burnin.add=NULL,
                                  MCMC.burnin=NULL, MCMC.burnin.mul=NULL,
                                  MCMC.prop.weights.form=NULL,MCMC.prop.args.form=NULL,
                                  MCMC.prop.weights.diss=NULL,MCMC.prop.args.diss=NULL,                                  
                                  MCMC.init.maxedges=NULL,
                                  MCMC.packagenames=NULL,

                                  MCMC.init.maxchanges=NULL){
    if(!is.null(MCMC.burnin) || !is.null(MCMC.burnin.mul)) stop("Control parameters MCMC.burnin and MCMC.burnin.mul are no longer used. See help for EGMME.MCMC.burnin.min, EGMME.MCMC.burnin.max, EGMME.MCMC.burnin.pval, EGMME.MCMC.burnin.pval, and CMLE.MCMC.burnin and CMLE.MCMC.interval for their replacements.")
    control<-list()
    for(arg in names(formals(sys.function())))
      control[arg]<-list(get(arg))

    set.control.class("control.simulate.stergm")
  }


