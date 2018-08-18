#  File R/simulate.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2018 Statnet Commons
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



#' Draw from the distribution of an Separable Temporal Exponential Family
#' Random Graph Model
#' 
#' \code{\link[stats]{simulate}} is used to draw from separable temporal
#' exponential family random network models in their natural parameterizations.
#' See \code{\link{stergm}} for more information on these models.
#' 
#' The dynamic process is run forward and the results are returned. For the
#' method for \code{\link{networkDynamic}}, the simulation is resumed from the
#' last generated time point of \code{object}, by default with the same model
#' and parameters.
#' 
#' The starting network for the \code{\link{stergm}} object method
#' (\code{simulate.stergm}) is determined by the \code{nw.start} argument.
#' 
#'   \itemize{
#'    \item{If \code{time.start} is specified, it is used as the initial
#'      time index of the simulation.}
#'    \item{If \code{time.start} is not specified (is \code{NULL}), then
#'      if the \code{object} carries a time stamp from which to start
#'      or resume the simulation, either in the form
#'      of a \code{"time"} network attribute (for the
#'      \code{\link{network}} method --- see the
#'      \code{\link{lasttoggle}} "API") or
#'      in the form of an \code{\link{net.obs.period}} network attribute (for the
#'      \code{\link{networkDynamic}} method), this attribute will be used. (If
#'      specified, \code{time.start} will override it with a warning.)
#'    }
#'    \item{Othewise, the simulation starts at 0.}
#'  }
#' 
#' @aliases simulate.stergm simulate.network simulate.networkDynamic
#' @param object an object of type \code{\link{stergm}} giving a model fit or
#' of type \code{\link{network}} giving the initial network.
#' 
#' \code{simulate.network} understands the \code{\link{lasttoggle}} "API".
#' @param formation,dissolution One-sided \code{\link{ergm}}-style formulas for
#' the formation and dissolution models, respectively.
#' @param nsim Number of replications (separate chains of networks) of the
#' process to run and return. The \code{\link{networkDynamic}} method only
#' supports \code{nsim=1}.
#' @param seed Random number integer seed.  See \code{\link[base]{set.seed}}.
#' @param coef.form Parameters for the model from which the post-formation
#' network is drawn.
#' @param coef.diss As \code{coef.form}, but for the post-dissolution network.
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
#' For STERGMs in particular, the constraints apply to the post-formation and
#' the post-dissolution network, rather than the final network. This means, for
#' example, that if the degree of all vertices is constrained to be less than
#' or equal to three, and a vertex begins a time step with three edges, then,
#' even if one of its edges is dissolved during its time step, it won't be able
#' to form another edge until the next time step. This behavior may change in
#' the future.
#' 
#' Note that not all possible combinations of constraints are supported.
#' @param monitor Either a one-sided formula specifying one or more terms whose
#' value is to be monitored, or a string containing \code{"formation"} or
#' \code{"dissolution"}, to monitor their respective terms, or \code{"all"} to
#' monitor distinct terms from both.
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
#' Constructed using \code{\link{control.simulate.stergm}} or
#' \code{\link{control.simulate.network}}.
#' @param statsonly Deprecated in favor of \code{output}.
#' @param output A character vector specifying output type: one of
#' "networkDynamic" (the default), "stats", and "changes", with partial
#' matching allowed. See Value section for details.
#' @param nw.start A specification for the starting network to be used by
#' \code{simulate.stergm}, optional for EGMME fits, but required for CMLE and
#' CMPLE fits: \describe{ \item{a numeric index }{use \code{i}th time-point's
#' network, where the first network in the series used to fit the model is
#' defined to be at the first time point;}
#' \item{`i`}{use \code{i}th
#' time-point's network, where the first network in the series used to fit the
#' model is defined to be at the first time point;}
#' \item{`"first"` or `"last"`}{the first or last time point used in
#' fitting the model; or}
#' \item{`network`}{specify the network directly.}
#' }
#' \code{\link{networkDynamic}}s
#' cannot be used as starting networks for \code{simulate.stergm} at this time.
#' (They can be used as starting networks for \code{simulate.networkDynamic},
#' of course.)
#' @param stats.form,stats.diss Logical: Whether to return
#' formation/dissolution model statistics. This is not the recommended method:
#' use \code{monitor} argument instead.
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
#'    \code{monitor} argument, and \code{stats.form} and \code{stats.diss}
#'    for the respective formation and dissolution statistics.
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
#'      \item{\code{formation}, \code{dissolution}, \code{monitor}:}{Formation, dissolution,
#'	and monitoring formulas	used in the simulation,	respectively.}
#'      \item{\code{stats}, \code{stats.form}, \code{stats.diss}:}{Network statistics as above.}
#'      \item{\code{coef.form}, \code{coef.diss}:}{Coefficients used in the simulation.}
#'      \item{\code{changes}:}{A four-column matrix summarizing the changes in the
#'	\code{"changes"} output. (This may be removed in the future.)}
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
#'      \item{formation, dissolution, monitor:}{Formation, dissolution,
#'	and monitoring formulas	used in the simulation,	respectively.}
#'      \item{stats, stats.form, stats.diss:}{Network statistics as above.}
#'      \item{coef.form, coef.diss:}{Coefficients used in the simulation.}
#'      \item{changes}{A four-column matrix summarizing the changes in the
#'	\code{"changes"} output. (This may be removed in the future.)}
#'    }
#'    When \code{nsim>1}, a \code{\link{network.list}} of these
#'    \code{\link{network}}s is returned.
#'  }
#' @examples
#' 
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
#' dynsim2<-simulate(dynsim,time.slices=S,verbose=TRUE)
#' @importFrom stats simulate
#' @export simulate.stergm
simulate.stergm<-function(object, nsim=1, seed=NULL,
                          coef.form=object$formation.fit$coef,coef.diss=object$dissolution.fit$coef,
                          constraints = object$constraints,
                          monitor = object$targets,
                          time.slices = 1, time.start=NULL, time.burnin=0, time.interval=1,
                          control=control.simulate.stergm(),
                          statsonly=NULL,
                          output=c("networkDynamic", "stats", "changes", "final"),
                          nw.start = NULL,
                          stats.form = FALSE,
                          stats.diss = FALSE,
                          duration.dependent = NULL,
                          verbose=FALSE, ...){
  .dep_method("simulate","stergm")
  check.control.class(c("simulate.stergm","simulate.network"), "simulate.stergm")
  
  control.transfer <- list(MCMC.prop.weights="MCMC.prop.weights",
                           MCMC.prop.args="MCMC.prop.args",
                           MCMC.packagenames="MCMC.packagenames",
                           MCMC.init.maxedges="MCMC.init.maxedges",
                           MCMC.init.maxchanges="MCMC.init.maxchanges",
                           EGMME.MCMC.burnin.min="MCMC.burnin.min",
                           EGMME.MCMC.burnin.max="MCMC.burnin.max",
                           EGMME.MCMC.burnin.pval="MCMC.burnin.pval",
                           EGMME.MCMC.burnin.add="MCMC.burnin.add")
  for(arg in names(control.transfer))
    if(is.null(control[[control.transfer[[arg]]]]))
      control[control.transfer[[arg]]] <- list(object$control[[arg]])

  control <- set.control.class("control.simulate.network")

  if(is.null(nw.start)){
    if(is.network(object$network)) nw.start <- object$network
    else stop('Simulating from STERGM CMLE fit requires the starting network to be specified in the nw.start argument: "first", "last", a numeric index of the network in the series (with "first"==1), or a network (NOT networkDynamic at this time).')
  }else if(is.numeric(nw.start)){
    nw.start <- object$network[[nw.start]]
    if(!is.network(nw.start)) stop("Invalid starting network specification.")
  }else if(is.character(nw.start)){
    nw.start <- switch(nw.start,
                       first = object$network[[1]],
                       last = object$network[[length(object$network)]],
                       stop("Invalid starting network specification."))
  }else if(is.networkDynamic(nw.start)){
    stop("Using a networkDynamic to start a simulation from a STERGM is not supported at this time.")
  }
    if(is.null(duration.dependent)){
    duration.dependent <-   is.lasttoggle(nw.start,object$formation,object$dissolution,object$monitor)
  }
  
  if(!duration.dependent)
    nw.start %n% "lasttoggle" <- NULL
  
  simulate.network(nw.start,formation=object$formation,dissolution=object$dissolution,nsim=nsim,coef.form=coef.form, coef.diss=coef.diss, constraints=constraints, monitor=monitor, time.start=time.start, time.slices=time.slices, time.burnin=time.burnin, time.interval=time.interval,control=control, statsonly=statsonly, output=output, stats.form = stats.form, stats.diss = stats.diss, duration.dependent=duration.dependent, verbose=verbose,...)
}



#' @rdname simulate.stergm
#' @export simulate.network
simulate.network <- function(object, nsim=1, seed=NULL,
                             formation, dissolution,
                             coef.form,coef.diss,
                             constraints = ~.,
                             monitor = NULL,
                             time.slices = 1, time.start=NULL, time.burnin=0, time.interval=1, time.offset=1,
                             control=control.simulate.network(),
                             statsonly=NULL,
                             output=c("networkDynamic", "stats", "changes", "final"),
                             stats.form = FALSE,
                             stats.diss = FALSE,
                             duration.dependent=NULL,
                             verbose=FALSE,...) {
  .dep_method("simulate","network")
  if(length(list(...))) stop("Unknown arguments: ",names(list(...)))
  check.control.class("simulate.network", "STERGM simulate.network")

  if(!is.null(statsonly)){
    warning("Argument `statsonly' for STERGM simulate() is deprecated and may be removed in a future version. Use `output' instead.")
    output <- if(statsonly) "stats" else "networkDynamic"
  }

  output <- match.arg(output)
  
  if(!is.null(seed)) set.seed(as.integer(seed))
  # output is a "main call" parameter, since it affects what to
  # compute rather than just how to compute it, but it's convenient to
  # have it as a part of the control data structure.
  if((time.burnin!=0 || time.interval!=1) && output!="stats" && output!="final"){
    stop("Generating a networkDynamic or change list output is incompatible with a time.burnin!=1 or a time.interval!=1. Only network statistics can be returned with these settings.")
  }
  
  control$changes <- output != "stats"

  if(output=="stats" && is.null(monitor) && !stats.form && !stats.diss){
    warning("Requesting a statistics matrix as output, but there are no statistics to return (monitor is NULL, and stats.form and stats.diss are both FALSE). Output will be empty.")
  }

  
  nw <- as.network(object)
  if(!is.network(nw)){
    stop("A network object must be given")
  }

  formation<-nonsimp_update.formula(formation,nw~., from.new="nw")
  dissolution<-nonsimp_update.formula(dissolution,nw~., from.new="nw")

  unset.offset.call <- function(object){
    if(inherits(object,"call") && object[[1]]=="offset")
      object[[2]]
    else
      object
  }
  
  if(is.character(monitor)){
    monitor <- switch(monitor,
                      formation = formation,
                      dissolution = dissolution,
                      all = append_rhs.formula(~nw, unique(lapply(c(term.list.formula(formation[[3]]),term.list.formula(dissolution[[3]])), unset.offset.call)))
                      )
  }
  
  if(!is.null(monitor)) monitor<-nonsimp_update.formula(monitor,nw~., from.new="nw")
  
  if(is.null(duration.dependent)){ # in the case of simulating from nw directly 
    duration.dependent <- is.lasttoggle(nw,formation,dissolution,monitor)
    if(duration.dependent)
      nw %n% "lasttoggle" <- NVL(nw %n% "lasttoggle",rep(round(-.Machine$integer.max/2), network.dyadcount(nw)))  else nw %n% "lasttoggle" <- NULL
      
  }



  model.form <- ergm_model(formation, nw, role="formation")
  if(!missing(coef.form) && nparam(model.form)!=length(coef.form)) stop("coef.form has ", length(coef.form), " elements, while the model requires ",nparam(model.form)," parameters.")

  model.diss <- ergm_model(dissolution, nw, role="dissolution")
  if(!missing(coef.diss) && nparam(model.diss)!=length(coef.diss)) stop("coef.diss has ", length(coef.diss), " elements, while the model requires ",nparam(model.diss)," parameters.")

  model.mon <- if(!is.null(monitor)) ergm_model(monitor, nw, role="target") else NULL
  
  if(missing(coef.form)) {
    coef.form <- rep(0,nparam(model.form, canonical=TRUE))
    warning("No parameter values given, using Bernouli formation.\nThis means that every time step, half the non-tie dyads will gain a tie!")
  }

  if(missing(coef.diss)) {
    coef.diss <- rep(0,nparam(model.diss, canonical=TRUE))
    warning("No parameter values given, using Bernoulli dissolution.\nThis means that every time step, half the ties get dissolved!\n")
  }
    
  proposal.form <- ergm_proposal(constraints,control$MCMC.prop.args.form,nw,
                                weights=control$MCMC.prop.weights.form,class="f")
  proposal.diss <- ergm_proposal(constraints,control$MCMC.prop.args.diss,nw,
                                weights=control$MCMC.prop.weights.diss,class="d")

  eta.form <- ergm.eta(coef.form, model.form$etamap)
  eta.diss <- ergm.eta(coef.diss, model.diss$etamap)
  
  control$time.burnin <- time.burnin
  control$time.interval <- time.interval
  control$time.samplesize <- time.slices
  control$collect.form <- stats.form
  control$collect.diss <- stats.diss
  
  out <- replicate(nsim, {
    nw <- .set.default.net.obs.period(nw, time.start)
    nw %n% "time" <- start <- .get.last.obs.time(nw, time.start)
    z <- stergm_MCMC_sample(nw, model.form, model.diss, model.mon,
                              proposal.form, proposal.diss,
                              eta.form=eta.form, eta.diss=eta.diss, control=control, verbose=verbose)
    
    stats.form <- if(control$collect.form) mcmc(sweep(z$statsmatrix.form,2,summary(formation),"+"),start=time.burnin+1,thin=time.interval)
    stats.diss <- if(control$collect.diss) mcmc(sweep(z$statsmatrix.diss,2,summary(dissolution),"+"),start=time.burnin+1,thin=time.interval)
    stats <- if(!is.null(model.mon)) mcmc(sweep(z$statsmatrix.mon,2,summary(monitor),"+"),start=time.burnin+1,thin=time.interval)

    out <-
      switch(output,
             networkDynamic = {
               changes <- z$changed
               # update the times on the list of changes returned to match the update time requested by time.offset
               # TODO: this is a horrible hack, should really make the change inside simulate.stergm
               # assume that simulate.stergm has added +1 to all the time values, so subtract 1 for an offset of 0
               changes[,1]<-changes[,1]-1+time.offset
               nwd <- to.networkDynamic.lasttoggle(nw)
               nwd <- networkDynamic.apply.changes(nwd,changes)
               attributes(nwd) <- c(attributes(nwd), # Don't clobber existing attributes!
                                    list(formation = formation,
                                         dissolution = dissolution,
                                         stats.form = stats.form,
                                         stats.diss = stats.diss,
                                         monitor = monitor,
                                         stats = stats,
                                         coef.form=coef.form,
                                         coef.diss=coef.diss,
                                         constraints=constraints,
                                         changes = changes))
               nwd <- .add.net.obs.period.spell(nwd, start-1+time.offset, time.slices)
               nwd
             },
             changes = {
               changes <- z$changed
               # todo: update the times on the list of changes returned to match the update time requested by time.offset
               # this is a horrible hack, should really make the change inside simulate.stergm
               # assume that simulate.stergm has added +1 to all the time values, so subtract 1 for an offset of 0
               changes[,1]<-changes[,1]-1+time.offset
               attributes(changes) <- c(attributes(changes), # Don't clobber existing attributes!
                                        list(formation = formation,
                                             dissolution = dissolution,
                                             stats.form = stats.form,
                                             stats.diss = stats.diss,
                                             monitor = monitor,
                                             constraints = constraints,
                                             stats = stats,
                                             coef.form=coef.form,
                                             coef.diss=coef.diss,
                                             start = nw%n%"time" + 0,
                                             end = nw%n%"time" + time.slices))
               changes
             },
             stats = {
               list(stats.form = stats.form,stats.diss = stats.diss, stats = stats)
             },
             final = {
               changes <- z$changed
               # update the times on the list of changes returned to match the update time requested by time.offset
               # TODO: this is a horrible hack, should really make the change inside simulate.stergm
               # assume that simulate.stergm has added +1 to all the time values, so subtract 1 for an offset of 0
               changes[,1]<-changes[,1]-1+time.offset
               newnw <- z$newnetwork
               attributes(newnw) <- c(attributes(newnw), # Don't clobber existing attributes!
                                      list(formation = formation,
                                           dissolution = dissolution,
                                           stats.form = stats.form,
                                           stats.diss = stats.diss,
                                           monitor = monitor,
                                           stats = stats,
                                           coef.form=coef.form,
                                           coef.diss=coef.diss,
                                           constraints=constraints,
                                           start = nw%n%"time" + 0,
                                           end = nw%n%"time" + time.slices,
                                           changes = changes))
               newnw
             })
  },
                   simplify = FALSE)

  if(nsim==1){
    out<-out[[1]]
    if(output == "stats"){
      for(name in names(out)) # Strip the unreturned stats matrices.
        if(is.null(out[[name]]))
          out[[name]] <- NULL
      if(length(out)==1)
        out <- out[[1]] # If there is only one, just return it.
    }
    out
  }else{
    switch(output,
           networkDynamic = {
             # If we've returned a list of networkDynamics, then it's a network list.
             # FIXME: Should we reserve network.list for serial network data?
             class(out) <- "network.list"
             out
           },
           stats = {
             # If we returned several mcmc() objects, merge them into mcmc.lists.
             outl <- list()
             for(name in names(out[[1]])){
               if(!is.null(out[[1]][[name]]))
                 outl[[name]] <- do.call(mcmc.list, lapply(out, "[[", name))
             }
             if(length(outl)==1)
               outl <- outl[[1]] # If there is only one, just return it.
             outl
           },
           changes = {
             out
           },
           final = {
             # If we've returned a list of networks, then it's a network list.
             # FIXME: Should we reserve network.list for serial network data?
             class(out) <- "network.list"
             out
           })
  }
}
#' @rdname simulate.stergm
#' @export simulate.networkDynamic
simulate.networkDynamic <- function(object, nsim=1, seed=NULL,
                                    formation = attr(object, "formation"), dissolution = attr(object, "dissolution"),
                                    coef.form = attr(object, "coef.form"), coef.diss = attr(object, "coef.diss"),
                                    constraints = NVL(attr(object, "constraints"),~.),
                                    monitor = attr(object, "monitor"),
                                    time.slices = 1, time.start=NULL, time.burnin=0, time.interval=1, time.offset=1,
                                    control=control.simulate.network(),
                                    statsonly=NULL,
                                    output=c("networkDynamic", "stats", "changes"),
                                    stats.form = FALSE,
                                    stats.diss = FALSE,
                                    duration.dependent = NULL,
                                    verbose=FALSE, ...){
  .dep_method("simulate","networkDynamic")
  
  if(nsim>1) stop("Simulating more than one chain of networks is not supported at this time. If you want to simulate over multiple time steps, use the time.slices argument.")

  if(!is.null(statsonly)){
    warning("Argument `statsonly' for STERGM simulate() is deprecated and may be removed in a future version. Use `output' instead.")
    output <- if(statsonly) "stats" else "networkDynamic"
  }

  # Resolve the starting time by setting the initial (implicit) net.obs.period.
  object <- .set.default.net.obs.period(object, time.start)
  start <- .get.last.obs.time(object, time.start)
  
  # if the network does not have a vertex pid, create one
  if(is.null(object%n%'vertex.pid')){
    set.vertex.attribute(object,'tergm_pid',1:network.size(object))
    object%n%'vertex.pid'<-'tergm_pid'
  }
  
  if(verbose) cat("extracting state of networkDynamic at time ",start,"\n")
  
  # extract nwd to nw
  
  if(is.null(duration.dependent))
    duration.dependent <- is.lasttoggle(object,formation,dissolution,monitor)
  
  nw <- network.extract.with.lasttoggle(object, start, duration.dependent)

  output <- match.arg(output)
  # get back a 'changes' matrix for the next sim step with columns 'time','tail','head','to'
  sim <- simulate.network(nw, nsim=1, seed=NULL,
                          formation=formation, dissolution=dissolution,
                          coef.form=coef.form,coef.diss=coef.diss,
                          constraints = constraints,
                          monitor = monitor,
                          time.slices=time.slices, time.start=time.start, time.burnin=time.burnin, time.interval=time.interval, time.offset=time.offset,
                          control=control,          
                          output=switch(output, networkDynamic = "changes", output),
                          stats.form = stats.form,
                          stats.diss = stats.diss,
                          duration.dependent=duration.dependent,
                          verbose=verbose, ...)
  
  ## Map the vertex IDs back to the original dynamic network from the static sim:
  #  using the pid methods to be safe for when pids are non-numeric
  # pids are either pre-existing, or set by network.extract.with.lasttoggle
  #sim[,"tail"] <- vActiveIDs[sim[,"tail"]]
  #sim[,"head"] <- vActiveIDs[sim[,"head"]]
  sim[,"tail"] <- get.vertex.id(object,get.vertex.pid(nw,sim[,"tail"]))
  sim[,"head"] <- get.vertex.id(object,get.vertex.pid(nw,sim[,"head"]))
  
  ## If all the user wants is statistics or a list of toggles, we are done.
  if(output!="networkDynamic") return(sim)

  if(verbose) cat("Updating networkDynamic ")
  
  object  <- networkDynamic.apply.changes(object, sim)
  # set up net.obs.period list to describe time period simulated
  object <- .add.net.obs.period.spell(object, start+time.offset-1, time.slices)
  
  if(verbose){
    obs<-(object%n%'net.obs.period')$observations
    cat("with simulated time: (",obs[[length(obs)]],").\n")
  }
  
  attributes(object) <- c(attributes(object), # Don't clobber existing attributes!
                          list(formation = nonsimp_update.formula(formation,nw~., from.new="nw"),
                               dissolution = nonsimp_update.formula(dissolution,nw~., from.new="nw"),
                               stats.form = rbind(if(isTRUE(attr(sim,"formation")==attr(object,"formation"))) attr(object,"stats.form"),attr(sim,"stats.form")),
                               stats.diss = rbind(if(isTRUE(attr(sim,"dissolution")==attr(object,"dissolution"))) attr(object,"stats.diss"),attr(sim,"stats.diss")),
                               monitor = attr(sim,"monitor"),
                               stats = rbind(if(isTRUE(attr(sim,"monitor")==attr(object,"monitor"))) attr(object,"stats"),attr(sim,"stats")),
                               coef.form=coef.form,
                               coef.diss=coef.diss,
                               constraints=constraints,
                               changes = rbind(attr(object,"changes"),matrix(c(sim), nrow=nrow(sim),ncol=ncol(sim),dimnames=list(rownames(sim),colnames(sim))))
                               ))
  object
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
