#  File R/simulate.tergm.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2022 Statnet Commons
################################################################################

#' Draw from the distribution of a Temporal Exponential Family
#' Random Graph Model
#' 
#' \code{\link[stats]{simulate}} is used to draw from temporal
#' exponential family random network models in their natural parameterizations.
#' See \code{\link{tergm}} for more information on these models.
#' 
#' The dynamic process is run forward and the results are returned. For the
#' method for \code{\link{networkDynamic}}, the simulation is resumed from the
#' last generated time point of \code{basis} (or the left hand side of \code{object}
#' if \code{basis} is missing), by default with the same model
#' and parameters.
#' 
#' The starting network for the \code{\link{tergm}} object method
#' (\code{simulate.tergm}) is determined by the \code{nw.start} argument.
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
#' @aliases simulate.tergm simulate_formula.network simulate_formula.networkDynamic
#' @param object for \code{simulate.tergm}, an object of type \code{\link{tergm}} giving a model fit;
#'   for \code{simulate_formula.network} and \code{simulate_formula.networkDynamic}, a formula specifying
#'   the model
#' 
#' \code{simulate_formula.network} understands the \code{\link{lasttoggle}} "API".
#' @param nsim Number of replications (separate chains of networks) of the
#' process to run and return. The \code{\link{networkDynamic}} method only
#' supports \code{nsim=1}.
#' @template seed
#' @param coef Parameters for the model.
#' @param constraints A one-sided formula specifying one or more constraints on
#' the support of the distribution of the networks being modeled. Multiple constraints may be given,
#' separated by \dQuote{+} operators.  Together with the model terms in the
#' formula and the reference measure, the constraints define the distribution
#' of networks being modeled.
#' 
#' The default is \code{~.}, for an unconstrained model.
#' 
#' See the [ERGM constraints][ergmConstraint] documentation for the
#' constraints implemented in the **[ergm][ergm-package]** package.
#' Other packages may add their own constraints.
#' 
#' Note that not all possible combinations of constraints are supported.
#' @param monitor A one-sided formula specifying one or more terms whose
#' value is to be monitored.  If \code{monitor} is specified as a character
#' (one of \code{"formation"}, \code{"dissolution"}, and \code{"all"}) then
#' the function \code{\link{.extract.fd.formulae}} is used to determine the
#' corresponding formula; the user should be aware of its behavior and limitations.
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
#' @param control A list of control parameters for algorithm tuning.
#' Constructed using \code{\link{control.simulate.tergm}} or
#' \code{\link{control.simulate.formula.tergm}}.  For backwards compatibility,
#' control lists from \code{\link{control.simulate.stergm}} and
#' \code{\link{control.simulate.network}} are allowed in calls to
#' \code{simulate.tergm}; they are mapped to \code{control.simulate.tergm}
#' by assigning:
#' \itemize{
#'   \item \code{MCMC.prop.form} to \code{MCMC.prop},
#'   \item \code{MCMC.prop.args.form} to \code{MCMC.prop.args},
#'   \item \code{MCMC.prop.weights.form} to \code{MCMC.prop.weights}.
#' }
#' @param output A character vector specifying output type: one of
#' `"networkDynamic"` (the default), `"stats"`, `"changes"`, `"final"`, and
#' `"ergm_state"`, with partial matching allowed. See Value section for details.
#' @param nw.start A specification for the starting network to be used by
#' \code{simulate.tergm}, optional for EGMME fits, but required for CMLE and
#' CMPLE fits: \describe{ \item{a numeric index `i`}{use \code{i}th time-point's
#' network, where the first network in the series used to fit the model is
#' defined to be at the first time point;}
#' \item{`"first"` or `"last"`}{the first or last time point used in
#' fitting the model; or}
#' \item{`network`}{specify the network directly.}
#' }
#' \code{\link{networkDynamic}}s
#' cannot be used as starting networks for \code{simulate.tergm} at this time.
#' (They can be used as starting networks for \code{simulate_formula.networkDynamic},
#' of course.)
#' @param stats Logical: Whether to return
#' model statistics. This is not the recommended method:
#' use \code{monitor} argument instead.
#' @template verbose
#' @param \dots Further arguments passed to or used by methods.
#' @param basis For the \code{network} and \code{networkDynamic} methods,
#'   the network to start the simulation from.  (If \code{basis} is missing,
#'   the default is the left hand side of the \code{object} argument.)
#' @param dynamic Logical; if \code{TRUE}, dynamic simulation is performed in
#'   \code{tergm}; if \code{FALSE} (the default), ordinary \code{ergm}
#'   simulation is performed instead.  Note that when \code{dynamic=FALSE},
#'   default argument values for \code{ergm}'s \code{simulate} methods
#'   are used.
#' @return Depends on the \code{output} argument:
#'  \item{"stats"}{If \code{stats == FALSE}, an \code{\link{mcmc}} matrix with
#'    monitored statistics, and if \code{stats == TRUE}, a
#'    list containing elements \code{monitor} for statistics specified in the
#'    \code{monitor} argument, and \code{model} for the model statistics.
#'    If \code{stats == FALSE} and no monitored statistics are specified,
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
#'    and monitoring formulas used in the simulation, respectively.}
#'      \item{\code{stats}:}{Network statistics as above.}
#'      \item{\code{coef}:}{Coefficients used in the simulation.}
#'      \item{\code{changes}:}{A four-column matrix summarizing the changes in the
#'    \code{"changes"} output. (This may be removed in the future.)}
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
#'      Additionally, the same attributes (\code{\link{attr}}, not network
#'      attributes) as with \code{output=="networkDynamic"} are attached.
#'      When \code{nsim>1}, a list of these change matrices is returned.}
#'  \item{"final"}{A \code{\link[network]{network}}
#'    object representing the last network in the series generated.
#'    \code{\link{lasttoggle}} and \code{time} attributes are also included.
#'    Additionally, the same attributes (\code{\link{attr}}, not network
#'    attributes) as with \code{output=="networkDynamic"} are attached.
#'    When \code{nsim>1}, a \code{\link{network.list}} of these
#'    \code{\link{network}}s is returned.
#'  }
#'  \item{"ergm_state"}{The \code{\link[ergm]{ergm_state}} object resulting
#'    from the simulation.  Attributes are attached as for other output types.}
#'  Note that when using \code{simulate_formula.networkDynamic} with either
#'    \code{"final"} or \code{"ergm_state"} for \code{output}, the nodes
#'    included in these objects are those produced by \code{network.collapse}
#'    at the start time.
#'
#' @examples
#' \donttest{
#' data(samplk)
#' 
#' # Fit a transition from Time 1 to Time 2
#' samplk12 <- tergm(list(samplk1, samplk2)~
#'                   Form(~edges+mutual+transitiveties+cyclicalties)+
#'                   Diss(~edges+mutual+transitiveties+cyclicalties),
#'                   estimate="CMLE")
#'
#' # direct simulation from tergm object
#' sim1 <- simulate(samplk12, nw.start="last")
#'
#' # equivalent simulation from formula with network LHS;
#' # must pass dynamic=TRUE for tergm simulation
#' sim2 <- simulate(samplk2 ~ Form(~edges+mutual+transitiveties+cyclicalties) +
#'                            Diss(~edges+mutual+transitiveties+cyclicalties),
#'                            coef = coef(samplk12),
#'                            dynamic=TRUE)
#'
#' # the default simulate output is a networkDynamic, and we can simulate
#' # with a networkDynamic LHS as well
#' sim3 <- simulate(sim2 ~ Form(~edges+mutual+transitiveties+cyclicalties) +
#'                         Diss(~edges+mutual+transitiveties+cyclicalties),
#'                         coef = coef(samplk12),
#'                         dynamic=TRUE)
#' }
#' @importFrom stats simulate
#' @export
simulate.tergm<-function(object, nsim=1, seed=NULL,
                          coef=coefficients(object),
                          constraints = object$constraints,
                          monitor = object$targets,
                          time.slices = 1, time.start=NULL, time.burnin=0, time.interval=1,
                          control=control.simulate.tergm(),
                          output=c("networkDynamic", "stats", "changes", "final", "ergm_state"),
                          nw.start = NULL,
                          stats = FALSE,
                          verbose=FALSE, ...){
  if(is(control, "control.simulate.stergm") || is(control, "control.simulate.network")) {
    control$MCMC.prop <- control$MCMC.prop.form
    control$MCMC.prop.weights <- control$MCMC.prop.weights.form
    control$MCMC.prop.args <- control$MCMC.prop.args.form
    control <- set.control.class("control.simulate.tergm")
  }  
  
  check.control.class(c("simulate.tergm","simulate.formula.tergm"), "simulate.tergm")
  
  control.transfer <- list(MCMC.prop.weights="MCMC.prop.weights",
                           MCMC.prop.args="MCMC.prop.args",
                           MCMC.prop="MCMC.prop",
                           MCMC.packagenames="MCMC.packagenames",
                           MCMC.maxedges="MCMC.maxedges",
                           MCMC.maxchanges="MCMC.maxchanges",
                           EGMME.MCMC.burnin.min="MCMC.burnin.min",
                           EGMME.MCMC.burnin.max="MCMC.burnin.max",
                           EGMME.MCMC.burnin.pval="MCMC.burnin.pval",
                           EGMME.MCMC.burnin.add="MCMC.burnin.add")
  for(arg in names(control.transfer))
    if(is.null(control[[control.transfer[[arg]]]]))
      control[control.transfer[[arg]]] <- list(object$control[[arg]])

  control <- set.control.class("control.simulate.formula.tergm")

  nw <- object$network
  if(is.null(nw.start)){
    if(is(object, "tergm_CMLE")) stop('Simulating from TERGM CMLE fit requires the starting network to be specified in the nw.start argument: "first", "last", a numeric index of the network in the series (with "first"==1), or a network (NOT networkDynamic at this time).')
    nw.start <- nw
  }else if(is.numeric(nw.start)){
    nwl <- uncombine_network(nw)
    if(nw.start == 1) nw.start <- (nwl[[1]] %n% ".PrevNets")[[1]]
    else nw.start <- nwl[[nw.start - 1]]
    if(!is.network(nw.start)) stop("Invalid starting network specification.")
  }else if(is.character(nw.start)){
    nwl <- uncombine_network(nw)
    nw.start <- switch(nw.start,
                       first = (nwl[[1]] %n% ".PrevNets")[[1]],
                       last = nwl[[length(nwl)]],
                       stop("Invalid starting network specification."))
    if(!is.network(nw.start)) stop("Invalid starting network specification.")                   
  }else if(is.networkDynamic(nw.start)){
    stop("Using a networkDynamic to start a simulation from a TERGM is not supported at this time.")
  }

  simulate_formula.network(object=object$formula, basis=nw.start,nsim=nsim,coef=coef, constraints=constraints, monitor=monitor, time.start=time.start, time.slices=time.slices, time.burnin=time.burnin, time.interval=time.interval,control=control, output=match.arg(output), stats=stats, verbose=verbose, dynamic=TRUE, ...)
}



#' @rdname simulate.tergm
#' @export
simulate_formula.network <- function(object, nsim=1, seed=NULL,
                             coef = NULL,
                             constraints = ~.,
                             monitor = NULL,
                             time.slices = 1, time.start=NULL, time.burnin=0, time.interval=1, time.offset=1,
                             control=control.simulate.formula.tergm(),
                             output=c("networkDynamic", "stats", "changes", "final", "ergm_state"),
                             stats = FALSE,
                             verbose=FALSE,
                             ..., basis=ergm.getnetwork(object), dynamic=FALSE) {
                             
  if(!dynamic) {
    if(missing(dynamic) && "lasttoggle" %in% list.network.attributes(basis)) warning("For dynamic simulation in ", sQuote("tergm"), " you must pass ", sQuote("dynamic=TRUE"), ".  Attempting ", sQuote("ergm"), " simulation instead...")

    mc <- match.call()
    mc[[1]] <- ergm_simulate_formula_network
    return(eval.parent(mc))
  }

  # reassign names for consistency with existing code
  formula <- object
  nw <- basis

  if(length(list(...))) stop("Unknown arguments: ",names(list(...)))
  check.control.class("simulate.formula.tergm", "TERGM simulate_formula.network")

  output <- match.arg(output)
  
  if(!is.null(seed)) set.seed(as.integer(seed))
  # output is a "main call" parameter, since it affects what to
  # compute rather than just how to compute it, but it's convenient to
  # have it as a part of the control data structure.
  if((time.burnin!=0 || time.interval!=1) && output!="stats" && output!="final"){
    stop("Only network statistics or the final network can be returned when time.burnin!=0 or time.interval!=1.")
  }
  
  control$changes <- output != "stats"

  if(output=="stats" && is.null(monitor) && !stats){
    warning("Requesting a statistics matrix as output, but there are no statistics to return (monitor is NULL, and stats is FALSE). Output will be empty.")
  }

  if(is.character(monitor)) {
    formula_pieces <- .extract.fd.formulae(formula)
    
    monitor <- switch(monitor,
                      formation = formula_pieces$form,
                      dissolution = formula_pieces$pers,
                      all = formula_pieces$all
                      )
  }

  proposal <- ergm_proposal(constraints, arguments = control$MCMC.prop.args, nw = nw,
                            weights = control$MCMC.prop.weights, hints = control$MCMC.prop, class="t")

  model <- ergm_model(formula, nw, dynamic=TRUE, term.options=control$term.options, extra.aux=list(proposal=proposal$auxiliaries, system=~.lasttoggle))
  proposal$aux.slots <- model$slots.extra.aux$proposal

  if(!missing(coef) && nparam(model)!=length(coef)) stop("coef has ", length(coef), " elements, while the model requires ",nparam(model)," parameters.")

  model.mon <- if(!is.null(monitor)) ergm_model(monitor, nw, dynamic=TRUE, term.options=control$term.options) else NULL
  
  if(missing(coef)) {
    coef <- rep(0,nparam(model, canonical=TRUE))
    warning("No parameter values given, using Bernoulli model.")
  }
    
  eta <- ergm.eta(coef, model$etamap)
  
  control$time.burnin <- time.burnin
  control$time.interval <- time.interval
  control$time.samplesize <- time.slices
  control$collect <- stats
  
  out <- replicate(nsim, {
    nw <- .set.default.net.obs.period(nw, time.start)
    nw %n% "time" <- start <- .get.last.obs.time(nw, time.start)
    z <- tergm_MCMC_sample(nw, model, model.mon,
                              proposal,
                              eta, control=control, verbose=verbose)
    
    stats.gen <- if(control$collect) mcmc(sweep(z$statsmatrix.gen,2,summary(model, nw=nw),"+"),start=time.burnin+1,thin=time.interval)
    stats.mon <- if(!is.null(model.mon)) mcmc(sweep(z$statsmatrix.mon,2,summary(model.mon, nw=nw),"+"),start=time.burnin+1,thin=time.interval)

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
                                    list(stats = list(model = stats.gen, monitor = stats.mon),
                                         coef = coef,
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
                                        list(stats = list(model = stats.gen, monitor = stats.mon),
                                             coef = coef,
                                             start = nw%n%"time" + 0,
                                             end = nw%n%"time" + time.slices))
               changes
             },
             stats = {
               list(model = stats.gen, monitor = stats.mon)
             },
             final = {
               changes <- z$changed
               # update the times on the list of changes returned to match the update time requested by time.offset
               # TODO: this is a horrible hack, should really make the change inside simulate.stergm
               # assume that simulate.stergm has added +1 to all the time values, so subtract 1 for an offset of 0
               changes[,1]<-changes[,1]-1+time.offset
               newnw <- as.network(z$newnetwork)
               newnw <- .add.net.obs.period.spell(newnw, start-1+time.offset, time.slices)
               attributes(newnw) <- c(attributes(newnw), # Don't clobber existing attributes!
                                      list(stats = list(model = stats.gen, monitor = stats.mon),
                                           coef = coef,
                                           start = nw%n%"time" + 0,
                                           end = nw%n%"time" + time.slices,
                                           changes = changes))
               newnw
             },
             ergm_state = {
               changes <- z$changed
               # update the times on the list of changes returned to match the update time requested by time.offset
               # TODO: this is a horrible hack, should really make the change inside simulate.stergm
               # assume that simulate.stergm has added +1 to all the time values, so subtract 1 for an offset of 0
               changes[,1]<-changes[,1]-1+time.offset
               newnw <- z$newnetwork
               attributes(newnw) <- c(attributes(newnw), # Don't clobber existing attributes!
                                      list(stats = list(model = stats.gen, monitor = stats.mon),
                                           coef = coef,
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
           },
           ergm_state = {
             out
           })
  }
}
#' @rdname simulate.tergm
#' @export
simulate_formula.networkDynamic <- function(object, nsim=1, seed=NULL,
                                    coef = attr(basis, "coef"),
                                    constraints = ~.,
                                    monitor = NULL,
                                    time.slices = 1, time.start=NULL, time.burnin=0, time.interval=1, time.offset=1,
                                    control=control.simulate.formula.tergm(),
                                    output=c("networkDynamic", "stats", "changes", "final", "ergm_state"),
                                    stats = FALSE,
                                    verbose=FALSE, ..., basis=eval_lhs.formula(object), dynamic=FALSE){

  if(!dynamic) {
    if(missing(dynamic)) warning("For dynamic simulation in ", sQuote("tergm"), " you must pass ", sQuote("dynamic=TRUE"), ".  Attempting ", sQuote("ergm"), " simulation instead...")

    mc <- match.call()
    mc[[1]] <- ergm_simulate_formula_network
    return(eval.parent(mc))
  }

  # reassign names for consistency with existing code  
  formula <- object
  object <- basis
    
  if(nsim>1) stop("Simulating more than one chain of networks is not supported at this time. If you want to simulate over multiple time steps, use the time.slices argument.")

  # Resolve the starting time by setting the initial (implicit) net.obs.period.
  object <- .set.default.net.obs.period(object, time.start)
  start <- .get.last.obs.time(object, time.start)
  
  # if the network does not have a vertex pid, create one
  if(is.null(object%n%'vertex.pid')){
    set.vertex.attribute(object,'tergm_pid',1:network.size(object))
    object%n%'vertex.pid'<-'tergm_pid'
  }
  
  if(verbose) message("extracting state of networkDynamic at time  ", start)
  
  # extract nwd to nw
  
  nw <- network.extract.with.lasttoggle(object, start)

  output <- match.arg(output)
  # get back a 'changes' matrix for the next sim step with columns 'time','tail','head','to'
  sim <- simulate_formula.network(object=formula, basis=nw, nsim=1, seed=NULL,
                          coef=coef,
                          constraints = constraints,
                          monitor = monitor,
                          time.slices=time.slices, time.start=time.start, time.burnin=time.burnin, time.interval=time.interval, time.offset=time.offset,
                          control=control,          
                          output=switch(output, networkDynamic = "changes", output),
                          stats = stats,
                          verbose=verbose, dynamic=dynamic, ...)
  
  ## Map the vertex IDs back to the original dynamic network from the static sim:
  #  using the pid methods to be safe for when pids are non-numeric
  # pids are either pre-existing, or set by network.extract.with.lasttoggle
  #sim[,"tail"] <- vActiveIDs[sim[,"tail"]]
  #sim[,"head"] <- vActiveIDs[sim[,"head"]]

  if(output %in% c("networkDynamic", "changes")) {
    sim[,"tail"] <- get.vertex.id(object,get.vertex.pid(nw,sim[,"tail"]))
    sim[,"head"] <- get.vertex.id(object,get.vertex.pid(nw,sim[,"head"]))
  }
  
  ## if user does not want a networkDynamic returned, we are done
  if(output!="networkDynamic") return(sim)

  if(verbose) message("Updating networkDynamic ", appendLF = FALSE)
  
  object  <- networkDynamic.apply.changes(object, sim)
  # set up net.obs.period list to describe time period simulated
  object <- .add.net.obs.period.spell(object, start+time.offset-1, time.slices)
  
  if(verbose){
    obs<-(object%n%'net.obs.period')$observations
    message("with simulated time: ( ", obs[[length(obs)]], " ).")
  }
  
  attributes(object) <- c(attributes(object), # Don't clobber existing attributes!
                          list(stats = attr(sim, "stats"),
                               coef = coef,
                               changes = rbind(attr(object,"changes"),matrix(c(sim), nrow=nrow(sim),ncol=ncol(sim),dimnames=list(rownames(sim),colnames(sim))))
                               ))
  object
}
