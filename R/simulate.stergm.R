#  File R/simulate.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
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
#' collect network statistics.
#' @param time.interval Number of time steps between successive recordings of
#' network statistics.
#' @param time.offset Argument specifying the offset between the point when the
#' state of the network is sampled (\code{time.start}) and the the beginning of
#' the spell that should be recorded for the newly simulated network state.
#' @param control A list of control parameters for algorithm tuning,
#' constructed using \code{\link{control.simulate.network}}.  These are mapped
#' to \code{\link{control.simulate.network.tergm}} controls by assigning:
#' \itemize{
#'   \item \code{MCMC.prop.form} to \code{MCMC.prop},
#'   \item \code{MCMC.prop.args.form} to \code{MCMC.prop.args},
#'   \item \code{MCMC.prop.weights.form} to \code{MCMC.prop.weights},
#'   \item \code{MCMC.init.maxedges} to \code{MCMC.maxedges}, and
#'   \item \code{MCMC.init.maxchanges} to \code{MCMC.maxchanges}.
#' }
#' @param output A character vector specifying output type: one of
#' \code{"networkDynamic"} (the default), \code{"stats"}, \code{"changes"},
#' \code{"final"}, and \code{"ergm_state"}, with partial matching allowed.
#' @param stats.form,stats.diss Logical: Whether to return
#' formation/dissolution model statistics. This is not the recommended method:
#' use the \code{monitor} argument instead.  Note that if either \code{stats.form}
#' or \code{stats.diss} is \code{TRUE}, all generative model statistics will be
#' returned.
#' @param verbose Logical: If \code{TRUE}, extra information is printed as the Markov
#' chain progresses.
#' @param \dots Further arguments passed to or used by methods.
#' @return Depends on the \code{output} argument.  See \code{\link{simulate.tergm}}
#'         for details.  Note that some formation/dissolution separated
#'         information is also attached to the return value for calls made through
#'         \code{simulate.network} and \code{simulate.networkDynamic} in
#'         an attempt to increase backwards compatibility.
#' @examples
#' logit<-function(p)log(p/(1-p))
#' coef.form.f<-function(coef.diss,density) -log(((1+exp(coef.diss))/(density/(1-density)))-1)
#' 
#' # Construct a network with 20 nodes and 20 edges
#' n<-20
#' target.stats<-edges<-20
#' g0<-network.initialize(n,dir=TRUE)
#' g1<-san(g0~edges,target.stats=target.stats,verbose=TRUE)
#' 
#' S<-10
#' 
#' # To get an average duration of 10...
#' duration<-10
#' coef.diss<-logit(1-1/duration)
#' 
#' # To get an average of 20 edges...
#' dyads<-network.dyadcount(g1)
#' density<-edges/dyads
#' coef.form<-coef.form.f(coef.diss,density)
#' 
#' # ... coefficients.
#' print(coef.form)
#' print(coef.diss)
#' 
#' # Simulate a networkDynamic
#' dynsim<-simulate(g1,formation=~edges,dissolution=~edges,
#'                  coef.form=coef.form,coef.diss=coef.diss,
#'                  time.slices=S,verbose=TRUE)
#' 
#' # "Resume" the simulation.
#' dynsim2<-simulate(dynsim,formation=~edges,dissolution=~edges,time.slices=S,verbose=TRUE)


#' @rdname simulate.stergm
#' @export
simulate.network <- function(object, nsim=1, seed=NULL,
                             formation, dissolution,
                             coef.form, coef.diss,
                             constraints = ~.,
                             monitor = NULL,
                             time.slices = 1, time.start=NULL, time.burnin=0, time.interval=1, time.offset=1,
                             control=control.simulate.network(),
                             output=c("networkDynamic", "stats", "changes", "final", "ergm_state"),
                             stats.form = FALSE,
                             stats.diss = FALSE,
                             verbose=FALSE,...) {
  check.control.class("simulate.network", "STERGM simulate.network")

  control$MCMC.prop <- control$MCMC.prop.form
  control$MCMC.prop.args <- control$MCMC.prop.args.form
  control$MCMC.prop.weights <- control$MCMC.prop.weights.form
  
  control$MCMC.maxedges <- control$MCMC.init.maxedges
  control$MCMC.maxchanges <- control$MCMC.init.maxchanges
  
  control <- set.control.class("control.simulate.network.tergm")

  output <- match.arg(output)

  rv <- simulate(trim_env(~Form(formation) + Diss(dissolution), keep = c("formation", "dissolution")), basis = object, nsim=nsim, seed=seed, coef = c(coef.form, coef.diss), constraints=constraints,
                 monitor=monitor, time.slices=time.slices, time.start=time.start, time.burnin=time.burnin, time.interval=time.interval, time.offset=time.offset, control=control, output=output, stats = stats.form || stats.diss, verbose=verbose, dynamic=TRUE, ...)
         
  if(output != "stats") {
    attributes(rv) <- c(attributes(rv), list(coef.form = coef.form, coef.diss = coef.diss))
    stats.gen <- attr(rv, "stats.gen")
    if(NCOL(stats.gen) > 0) {
      attr(rv, "stats.form") <- stats.gen[,grepl("Form", colnames(stats.gen)),drop=FALSE]
      attr(rv, "stats.diss") <- stats.gen[,grepl("Diss", colnames(stats.gen)),drop=FALSE]    
    }
  } else {
    stats.gen <- rv$stats.gen
    if(NCOL(stats.gen) > 0) {
      rv$stats.form <- stats.gen[,grepl("Form", colnames(stats.gen)),drop=FALSE]
      rv$stats.diss <- stats.gen[,grepl("Diss", colnames(stats.gen)),drop=FALSE]
    }
  }
  
  rv
}
#' @rdname simulate.stergm
#' @export
simulate.networkDynamic <- function(object, nsim=1, seed=NULL,
                                    formation, dissolution,
                                    coef.form = attr(object, "coef.form"), coef.diss = attr(object, "coef.diss"),
                                    constraints = ~.,
                                    monitor = NULL,
                                    time.slices = 1, time.start=NULL, time.burnin=0, time.interval=1, time.offset=1,
                                    control=control.simulate.network(),
                                    output=c("networkDynamic", "stats", "changes", "final", "ergm_state"),
                                    stats.form = FALSE,
                                    stats.diss = FALSE,
                                    verbose=FALSE, ...){
  check.control.class("simulate.network", "STERGM simulate.network")

  control$MCMC.prop <- control$MCMC.prop.form
  control$MCMC.prop.args <- control$MCMC.prop.args.form
  control$MCMC.prop.weights <- control$MCMC.prop.weights.form
  
  control$MCMC.maxedges <- control$MCMC.init.maxedges
  control$MCMC.maxchanges <- control$MCMC.init.maxchanges
  
  control <- set.control.class("control.simulate.network.tergm")

  output <- match.arg(output)

  rv <- simulate(trim_env(~Form(formation) + Diss(dissolution), keep = c("formation", "dissolution")), basis = object, nsim=nsim, seed=seed, coef = c(coef.form, coef.diss), constraints=constraints,
                 monitor=monitor, time.slices=time.slices, time.start=time.start, time.burnin=time.burnin, time.interval=time.interval, time.offset=time.offset, control=control, output=output, stats = stats.form || stats.diss, verbose=verbose, dynamic=TRUE, ...)

  if(output != "stats") {
    attributes(rv) <- c(attributes(rv), list(coef.form = coef.form, coef.diss = coef.diss))
    stats.gen <- attr(rv, "stats.gen")
    if(NCOL(stats.gen) > 0) {
      attr(rv, "stats.form") <- stats.gen[,grepl("Form", colnames(stats.gen)),drop=FALSE]
      attr(rv, "stats.diss") <- stats.gen[,grepl("Diss", colnames(stats.gen)),drop=FALSE]    
    }
  } else {
    stats.gen <- rv$stats.gen
    if(NCOL(stats.gen) > 0) {
      rv$stats.form <- stats.gen[,grepl("Form", colnames(stats.gen)),drop=FALSE]
      rv$stats.diss <- stats.gen[,grepl("Diss", colnames(stats.gen)),drop=FALSE]
    }
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
