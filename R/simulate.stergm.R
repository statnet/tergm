#  File R/simulate.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2019 Statnet Commons
#######################################################################
#========================================================================
# This file contains the following 3 functions for simulating stergms
#           <simulate.stergm>
#           <simulate.network>
#           <simulate.networkDynamic>
#========================================================================


########################################################################
# Each of the <simulate.X> functions collects a given number of networks
# drawn from the given distribution on the set of all networks; these
# may be returned as only the vector/matrix of sufficient statistics or
# as the networks and their statistics
#
# --PARAMETERS--
#   object       : either a stergm or a formation formula of the form
#                  'nw ~ term(s)'
#   dissolution  : for a formula 'object', this is the corresponding
#                  dissolution formula
#   nsim         : the number of networks to draw; default=1
#   seed         : an integer at which to set the random generator;
#                  default=NULL
#   theta.form   : the initial theta formation coefficients;
#                  default='object'$coef.form for stergm objects and
#                  default= a vector of 0's for formula objects
#   theta.diss   : the initial theta dissolution coefficients;
#                  default='object'$coef.diss for stergm objects and
#                  default= a vector of 0's for formula objects
#   time.burnin  : the number of MCMC steps to disregard before any MCMC
#                  sampling is done; default=0
#   time.interval: the number of MCMC steps between sampled networks;
#                  default=1
#   constraints  : a one-sided formula specifying the constraints on the
#                  support of the distribution of networks being simulated;
#                  default='object'$constraints for stergms, "~." for formulas
#   control      : a list of control parameters for algorithm tuning;
#                  default=<control.simulate.stergm>
#   changes      : whether 'changed', the change matrix of timestamps and
#                  changes, should be included in the return list (T or F);
#                  'changes' will be switched to FALSE if either of
#                  'time.burnin' or 'time.interval' do not have their default
#                  values; default=TRUE
#   duration.dependent : whether the model/formula are durational dependent
#                  default=NULL
#   verbose      : whether to print out information on the status of
#                  the simulations; default=FALSE
#
# --RETURNED--
#   only 
#      nw:  the final network from the simulation routine
#   if 'control$final'=TRUE (the default is FALSE)
#   otherwise
#     outlist: a network.list object as a list containing:
#        formation  : the formation formula
#        dissolution: the dissolution formula
#        coef.form  : the passed in or defaulted 'coef.form'
#        coef.diss  : the passed in or defaulted 'coef.diss'
#        networks   : the list of simulated networks
#        constraints: the constraints formula
#        stats.form : the matrix of sampled statistics for 'model.form'
#        stats.diss : the matrix of sampled statistics for 'model.form'
#        changed    : a toggle matrix, where the first column is
#                     the timestamp of the toggle and the 2nd and 3rd
#                     columns are the head & tail of the toggle; this
#                     is only returned if the input param 'changes'
#                     ends up being TRUE (see above) 
#                    'changed' will also have 2 attributes:
#            start  : 1
#            end    : the number of simulations
#        maxchanges : the size of "MCMC Dyn workspace"
#
###############################################################################



#' STERGM wrappers for TERGM simulation
#' 
#' The \code{simulate.network} and \code{simulate.networkDynamic} wrappers
#' are provided for backwards compatibility.  It is recommended that new
#' code make use of the \code{simulate_formula.network} and
#' \code{simulate_formula.networkDynamic} functions instead.  See
#' \code{\link{simulate.tergm}} for details on these new functions.
#'
#' Note that return values may be structured differently than in past versions.
#' 
#' @aliases simulate.stergm simulate.network simulate.networkDynamic
#'
#' @param object an object of type \code{\link{network}} or \code{\link{networkDynamic}}
#' @param formation,dissolution One-sided \code{\link{ergm}}-style formulas for
#' the formation and dissolution models, respectively.
#' @param nsim Number of replications (separate chains of networks) of the
#' process to run and return. The \code{\link{networkDynamic}} method only
#' supports \code{nsim=1}.
#' @param seed Random number integer seed.  See \code{\link[base]{set.seed}}.
#' @param coef.form Parameters for the formation model.
#' @param coef.diss Parameters for the dissolution model.
#' @param constraints A one-sided formula specifying one or more constraints on
#' the support of the distribution of the networks being modeled, using syntax
#' similar to the \code{formula} argument. Multiple constraints may be given,
#' separated by \dQuote{+} operators.  Together with the model terms in the
#' formula and the reference measure, the constraints define the distribution
#' of networks being modeled.
#' 
#' It is also possible to specify a proposal function directly by passing a
#' string with the function's name. In that case, arguments to the proposal
#' should be specified through the \code{prop.args} argument to
#' \code{\link{control.ergm}}.
#' 
#' The default is \code{~.}, for an unconstrained model.
#' 
#' See the [ERGM constraints][ergm-constraints] documentation for the
#' constraints implemented in the **[ergm][ergm-package]** package.
#' Other packages may add their own constraints.
#' 
#' Note that not all possible combinations of constraints are supported.
#' @param monitor A one-sided formula specifying one or more terms whose
#' value is to be monitored.
#' @param time.slices Number of time slices (or statistics) to return from each
#' replication of the dynamic process. See below for return types. Defaults to
#' 1, which, if \code{time.burnin==0} and \code{time.interval==1} (the
#' defaults), advances the process one time step.
#' @param time.start An optional argument specifying the time point at which
#' the simulation is to start. See Details for further information.
#' @param time.burnin Number of time steps to discard before starting to
#' collect network statistics. Actual network will only be returned if
#' \code{time.burnin==0}.
#' @param time.interval Number of time steps between successive recordings of
#' network statistics. Actual network will only be returned if
#' \code{time.interval==1}.
#' @param time.offset Argument specifying the offset between the point when the
#' state of the network is sampled (\code{time.start}) and the the beginning of
#' the spell that should be recorded for the newly simulated network state.
#' @param control A list of control parameters for algorithm tuning.
#' Constructed using \code{\link{control.simulate.network}}.  These are mapped
#' to \code{\link{control.simulate.network.tergm}} controls by assigning
#' \code{MCMC.prop.args.form} and \code{MCMC.prop.weights.form} to 
#' \code{MCMC.prop.args} and \code{MCMC.prop.weights} respectively, and assigning
#' \code{MCMC.init.maxedges} and \code{MCMC.init.maxchanges} to \code{MCMC.maxedges}
#' and \code{MCMC.maxchanges} respectively.
#' @param output A character vector specifying output type: one of
#' "networkDynamic" (the default), "stats", and "changes", with partial
#' matching allowed. See Value section for details.
#' @param stats.form,stats.diss Logical: Whether to return
#' formation/dissolution model statistics. This is not the recommended method:
#' use \code{monitor} argument instead.  Note that if either \code{stats.form}
#' or \code{stats.diss} is \code{TRUE}, all generative model statistics will be
#' returned as a single matrix.
#' @param duration.dependent Logical: Whether the model terms in formula or
#' model are duration dependent. E.g., if a duration-dependent term is used in
#' estimation/simulation model, the probability of forming or dissolving a tie
#' may dependent on the age the dyad status.
#' @param verbose Logical: If TRUE, extra information is printed as the Markov
#' chain progresses.
#' @param \dots Further arguments passed to or used by methods.
#' @return Depends on the \code{output} argument:
#'  \item{"stats"}{If \code{stats.form==FALSE} and
#'    \code{stats.diss==FALSE}, an \code{\link{mcmc}} matrix with
#'    monitored statistics, and if either of them is \code{TRUE}, a
#'    list containing elements \code{stats} for statistics specified in the
#'    \code{monitor} argument, and \code{stats.fd}
#'    for the combined formation and dissolution statistics.
#'    If \code{stats.form==FALSE} and
#'    \code{stats.diss==FALSE} and no monitored statistics are specified,
#'    an empty list is returned, with a warning.
#'    When \code{nsim>1}, an \code{\link{mcmc.list}} (or list of them) of
#'    the statistics is returned instead.}
#'  \item{"networkDynamic"}{A \code{\link[networkDynamic]{networkDynamic}}
#'    object representing the simulated process, with ties present in the
#'    initial network having onset \code{-Inf} and ties present at the end
#'    of the simulation having terminus \code{+Inf}. The method for
#'    \code{\link[networkDynamic]{networkDynamic}} returns the initial
#'    \code{\link[networkDynamic]{networkDynamic}} with simulated changes
#'    applied to it. The \code{\link{net.obs.period}} network attribute is
#'    updated (or added if not existing) to reflect the time period that was
#'    simulated. If the network does not have any \code{\link[networkDynamic]{persistent.ids}} 
#'    defined for vertices, a vertex.pid will be attached in a vertex attribute 
#'    named \code{'tergm_pid'} to facilitate 'bookkeeping' between the networkDynamic 
#'    argument and the simulated network time step.
#'    Additionally, attributes (\code{\link{attr}}, not network
#'    attributes) are attached as follows:
#'    \describe{
#'      \item{\code{formula}, \code{monitor}:}{Model
#'	and monitoring formulas	used in the simulation,	respectively.}
#'      \item{\code{stats}, \code{stats.fd}:}{Network statistics as above.}
#'      \item{\code{coef}:}{Coefficients used in the simulation.}
#'      \item{\code{changes}:}{A four-column matrix summarizing the changes in the
#'	\code{"changes"} output. (This may be removed in the future.)}
#'      \item{\code{formation}, \code{dissolution}:}{Formation and dissolution
#'	formula arguments.}
#'      \item{\code{coef.form}, \code{coef.diss}:}{Coefficient arguments.}
#'    } 
#'    When \code{nsim>1}, a \code{\link{network.list}} of these
#'    \code{\link{networkDynamic}}s is returned.
#'  }
#'  
#'  \item{"changes"}{An integer matrix with four columns (\code{time},
#'      \code{tail}, \code{head}, and \code{to}), giving the time-stamped
#'      changes relative to the current network. \code{to} is \code{1} if
#'      a tie was formed and \code{0} if a tie was dissolved. The
#'      convention for \code{time} is that it gives the time point during
#'      which the change is effective. For example, a row
#'      \code{c(5,2,3,1)} indicates that between time \eqn{4} and \eqn{5},
#'      a tie from node \eqn{2} to node \eqn{3} was formed, so that it was
#'      absent at time point \eqn{4} and present at time point \eqn{5};
#'      while a row \code{c(5,2,3,0)} indicates that in that time, that
#'      tie was dissolved, so that it is was present at time point \eqn{4}
#'      and absent at time point \eqn{5}.
#'      Additionally, same attributes (\code{\link{attr}}, not network
#'      attributes) as with \code{output=="networkDynamic"} are attached.
#'      When \code{nsim>1}, a list of these change matrices is returned.}
#'    \item{"final"}{A \code{\link[network]{network}}
#'    object representing the last network in the series generated. This
#'    is not implemented in the method for
#'    \code{\link[networkDynamic]{networkDynamic}}.
#'    \code{\link{lasttoggle}} attributes are also included.
#'    Additionally, attributes (\code{\link{attr}}, not network
#'    attributes) are attached as follows:
#'    \describe{
#'      \item{formula, monitor:}{Model
#'	and monitoring formulas	used in the simulation,	respectively.}
#'      \item{stats, stats.fd:}{Network statistics as above.}
#'      \item{coef:}{Coefficients used in the simulation.}
#'      \item{changes}{A four-column matrix summarizing the changes in the
#'	\code{"changes"} output. (This may be removed in the future.)}
#'      \item{\code{formation}, \code{dissolution}:}{Formation and dissolution
#'	formula arguments.}
#'      \item{\code{coef.form}, \code{coef.diss}:}{Coefficient arguments.}
#'    }
#'    When \code{nsim>1}, a \code{\link{network.list}} of these
#'    \code{\link{network}}s is returned.
#'  }



#' @rdname simulate.stergm
#' @export
simulate.network <- function(object, nsim=1, seed=NULL,
                             formation, dissolution,
                             coef.form,coef.diss,
                             constraints = ~.,
                             monitor = NULL,
                             time.slices = 1, time.start=NULL, time.burnin=0, time.interval=1, time.offset=1,
                             control=control.simulate.network(),
                             output=c("networkDynamic", "stats", "changes", "final"),
                             stats.form = FALSE,
                             stats.diss = FALSE,
                             duration.dependent=NULL,
                             verbose=FALSE,...) {
  check.control.class("simulate.network", "STERGM simulate.network")

  control$MCMC.prop.args <- control$MCMC.prop.args.form
  control$MCMC.prop.weights <- control$MCMC.prop.weights.form
  
  control$MCMC.maxedges <- control$MCMC.init.maxedges
  control$MCMC.maxchanges <- control$MCMC.init.maxchanges
  
  control <- set.control.class("control.simulate.network.tergm")

  rv <- simulate(object ~ FormE(formation) + DissE(dissolution), nsim=nsim, seed=seed, coef = c(coef.form, coef.diss), constraints=constraints,
           monitor=monitor, time.slices=time.slices, time.start=time.start, time.burnin=time.burnin, time.interval=time.interval, time.offset=time.offset, control=control, output=match.arg(output), stats = stats.form || stats.diss, duration.dependent=duration.dependent, verbose=verbose, dynamic=TRUE, ...)
           
  if(is(rv, "network")) {
    attributes(rv) <- c(attributes(rv), list(formation = formation, dissolution = dissolution, coef.form = coef.form, coef.diss = coef.diss))
  }
  
  rv
}
#' @rdname simulate.stergm
#' @export
simulate.networkDynamic <- function(object, nsim=1, seed=NULL,
                                    formation = attr(object, "formation"), dissolution = attr(object, "dissolution"),
                                    coef.form = attr(object, "coef.form"), coef.diss = attr(object, "coef.diss"),
                                    constraints = NVL(attr(object, "constraints"),~.),
                                    monitor = attr(object, "monitor"),
                                    time.slices = 1, time.start=NULL, time.burnin=0, time.interval=1, time.offset=1,
                                    control=control.simulate.network(),
                                    output=c("networkDynamic", "stats", "changes"),
                                    stats.form = FALSE,
                                    stats.diss = FALSE,
                                    duration.dependent = NULL,
                                    verbose=FALSE, ...){
  check.control.class("simulate.network", "STERGM simulate.network")

  control$MCMC.prop.args <- control$MCMC.prop.args.form
  control$MCMC.prop.weights <- control$MCMC.prop.weights.form
  
  control$MCMC.maxedges <- control$MCMC.init.maxedges
  control$MCMC.maxchanges <- control$MCMC.init.maxchanges
  
  control <- set.control.class("control.simulate.network.tergm")



  rv <- simulate(object ~ FormE(formation) + DissE(dissolution), nsim=nsim, seed=seed, coef = c(coef.form, coef.diss), constraints=constraints,
           monitor=monitor, time.slices=time.slices, time.start=time.start, time.burnin=time.burnin, time.interval=time.interval, time.offset=time.offset, control=control, output=match.arg(output), stats = stats.form || stats.diss, duration.dependent=duration.dependent, verbose=verbose, dynamic=TRUE, ...)

  if(is(rv, "networkDynamic")) {
    attributes(rv) <- c(attributes(rv), list(formation = formation, dissolution = dissolution, coef.form = coef.form, coef.diss = coef.diss))
  }
  
  rv
}

.set.default.net.obs.period <- function(nw, time.start=NULL){
  # get net.obs.period from nw, if it exists
  if (!is.null(nw%n%'net.obs.period')) return(nw)
  
  nwtime <- nw %n% "time"
  delete.network.attribute(nw,"time")
  
  nwtime <- if(is.null(nwtime)) NVL(time.start,0) else{
    if(!is.null(time.start)){
      if(time.start!=nwtime) warning("Argument time.start of ",time.start," specified for a network that already has a time stamp of ",nwtime, ". Overriding the time stamp.")
      time.start
    }else nwtime
  }

  set.network.attribute(nw, 'net.obs.period', list(observations=list(c(nwtime,nwtime+1)),mode="discrete",time.increment=1,time.unit="step"))
}

.get.last.obs.time <- function(nw, time.user=NULL){
  # get net.obs.period from nw
  net.obs.period<-nw%n%'net.obs.period'
  spells <- do.call(rbind,net.obs.period$observations)
  last.spell <- spells[which.max(apply(spells,1,mean)),]
  nwend <-
    if(last.spell[1]==last.spell[2] || net.obs.period$mode=="continuous") last.spell[2]
  # If in discrete mode and the last spell is not a point spell, then
  # find the latest time point that is an integer number of
  # time.icrements away from the onset, while still being strictly
  # less than terminus. For example, with time.interval=1, c(0,2) -> 1, c(0, 1.5) -> 1.
    else last.spell[1]+ceiling((last.spell[2]-last.spell[1])/net.obs.period$time.increment-1)*net.obs.period$time.increment

  
  if(!is.null(time.user)){
    if(time.user<nwend & nwend!=Inf) stop("Attempting to resume from a time point prior to the end of the previous simulation is not supported at this time.", call.=FALSE)
    if(time.user>nwend) warning("Argument time.start of ", time.user," specified for a network that already has a time stamp of ",nwend,". Overriding the time stamp.", call.=FALSE)
    time.user
  }else nwend
}

.get.first.obs.time <- function(nw){
  # get net.obs.period from nw
  net.obs.period<-nw%n%'net.obs.period'
  spells <- do.call(rbind,net.obs.period$observations)
  first.spell <- spells[which.min(apply(spells,1,mean)),]
  nwstart <- first.spell[1]
}

# add another observation spell to the end; note that the first *simulated* network is at start+1
.add.net.obs.period.spell <- function(nw, time.start, time.steps){
  nop <- nw%n%'net.obs.period'
  # make sure we don't add a duplicate of existing last spell
  # why didn't I (skye) define $observations as a spell list so we could use existing tools?  ... sigh
  newSpl<-c(time.start+1,time.start+time.steps+1)
  if (length(nop$observations)>0 && all(nop$observations[[length(nop$observations)]]==newSpl)){
    # don't do anything
  } else {
    # tack the new spell on the end of the observations list
    nop$observations<-c(nop$observations,list(newSpl))
  }
  set.network.attribute(nw, 'net.obs.period', nop)
}
