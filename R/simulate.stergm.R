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
                          verbose=FALSE, ...){
  check.control.class(c("simulate.stergm","simulate.network"))
  
  control.transfer <- c("MCMC.burnin", "MCMC.burnin.mul", "MCMC.prop.weights", "MCMC.prop.args", "MCMC.packagenames", "MCMC.init.maxedges", "MCMC.init.maxchanges")
  for(arg in control.transfer)
    if(is.null(control[[arg]]))
      control[arg] <- list(object$control[[arg]])

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
  
  simulate.network(nw.start,formation=object$formation,dissolution=object$dissolution,nsim=nsim,coef.form=coef.form, coef.diss=coef.diss, constraints=constraints, monitor=monitor, time.start=time.start, time.slices=time.slices, time.burnin=time.burnin, time.interval=time.interval,control=control, statsonly=statsonly, output=output, stats.form = stats.form, stats.diss = stats.diss, verbose=verbose,...)
}



# Note that we are overriding simulate.network here, since the first argument is a network.
simulate.network <- function(object, nsim=1, seed=NULL,
                             formation, dissolution,
                             coef.form,coef.diss,
                             constraints = ~.,
                             monitor = NULL,
                             time.slices = 1, time.start=NULL, time.burnin=0, time.interval=1,
                             control=control.simulate.network(),
                             statsonly=NULL,
                             output=c("networkDynamic", "stats", "changes", "final"),
                             stats.form = FALSE,
                             stats.diss = FALSE,
                             verbose=FALSE, ...) {
  if(length(list(...))) stop("Unknown arguments: ",names(list(...)))
  check.control.class(myname="STERGM simulate.network")

  if(!is.null(statsonly)){
    warning("Argument `statsonly' for STERGM simulate() is deprecated and may be removed in a future version. Use `output' instead.")
    output <- if(statsonly) "stats" else "networkDynamic"
  }

  output <- match.arg(output)
  
  if(!is.null(seed)) set.seed(as.integer(seed))
  # output is a "main call" parameter, since it affects what to
  # compute rather than just how to compute it, but it's convenient to
  # have it as a part of the control data structure.
  if((time.burnin!=0 || time.interval!=1) && output!="stats"){
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

  formation<-ergm.update.formula(formation,nw~.)
  dissolution<-ergm.update.formula(dissolution,nw~.)

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
                      all = append.rhs.formula(~nw, unique(lapply(c(term.list.formula(formation[[3]]),term.list.formula(dissolution[[3]])), unset.offset.call)))
                      )
  }
  
  if(!is.null(monitor)) monitor<-ergm.update.formula(monitor,nw~.)
  
  model.form <- ergm.getmodel(formation, nw, role="formation")
  if(!missing(coef.form) && coef.length.model(model.form)!=length(coef.form)) stop("coef.form has ", length(coef.form), " elements, while the model requires ",coef.length.model(model.form)," parameters.")

  model.diss <- ergm.getmodel(dissolution, nw, role="dissolution")
  if(!missing(coef.diss) && coef.length.model(model.diss)!=length(coef.diss)) stop("coef.diss has ", length(coef.diss), " elements, while the model requires ",coef.length.model(model.diss)," parameters.")

  model.mon <- if(!is.null(monitor)) ergm.getmodel(monitor, nw, role="target") else NULL
  
  verbose <- match(verbose,
                c("FALSE","TRUE", "very"), nomatch=1)-1
  if(missing(coef.form)) {
    coef.form <- rep(0,length(model.form$coef.names))
    warning("No parameter values given, using Bernouli formation.\nThis means that every time step, half the non-tie dyads will gain a tie!")
  }

  if(missing(coef.diss)) {
    coef.diss <- rep(0,length(model.diss$coef.names))
    warning("No parameter values given, using Bernoulli dissolution.\nThis means that every time step, half the ties get dissolved!\n")
  }
    
  MHproposal.form <- MHproposal(constraints,control$MCMC.prop.args.form,nw,
                                weights=control$MCMC.prop.weights.form,class="f")
  MHproposal.diss <- MHproposal(constraints,control$MCMC.prop.args.diss,nw,
                                weights=control$MCMC.prop.weights.diss,class="d")

  eta.form <- ergm.eta(coef.form, model.form$etamap)
  eta.diss <- ergm.eta(coef.diss, model.diss$etamap)
  
  control$time.burnin <- time.burnin
  control$time.interval <- time.interval
  control$time.samplesize <- time.slices
  control$collect.form <- stats.form
  control$collect.diss <- stats.diss
  
  out <- replicate(nsim, {
    if(is.null(nw %n% "lasttoggle")) nw %n% "lasttoggle" <- rep(round(-.Machine$integer.max/2), network.dyadcount(nw))
    nwtime <- nw %n% "time"

    nw %n% "time" <- if(is.null(nwtime)) NVL(time.start,0) else{
      if(!is.null(time.start)){
        if(time.start!=nwtime) warning("Argument time.start specified for a network that already has a time stamp. Overriding the time stamp.")
        time.start
      }else nwtime
    }
    
    
    z <- stergm.getMCMCsample(nw, model.form, model.diss, model.mon,
                              MHproposal.form, MHproposal.diss,
                              eta.form, eta.diss, control, verbose)
    
    stats.form <- if(control$collect.form) mcmc(sweep(z$statsmatrix.form,2,summary(formation),"+"),start=time.burnin+1,thin=time.interval)
    stats.diss <- if(control$collect.diss) mcmc(sweep(z$statsmatrix.diss,2,summary(dissolution),"+"),start=time.burnin+1,thin=time.interval)
    stats <- if(!is.null(model.mon)) mcmc(sweep(z$statsmatrix.mon,2,summary(monitor),"+"),start=time.burnin+1,thin=time.interval)

    out <-
      switch(output,
             networkDynamic = {
               library(networkDynamic)
               nwd <- to.networkDynamic.lasttoggle(nw)
               nwd <- networkDynamic.apply.changes(nwd,z$changed)
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
                                         start = nw%n%"time" + 0,
                                         end = nw%n%"time" + time.slices,
                                         changes = z$changed))
               nwd
             },
             changes = {
               changes <- z$changed
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
                                           changes = z$changed))
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

simulate.networkDynamic <- function(object, nsim=1, seed=NULL,
                                    formation = attr(object, "formation"), dissolution = attr(object, "dissolution"),
                                    coef.form = attr(object, "coef.form"), coef.diss = attr(object, "coef.diss"),
                                    constraints = NVL(attr(object, "constraints"),~.),
                                    monitor = attr(object, "monitor"),
                                    time.slices = 1, time.start=NULL, time.burnin=0, time.interval=1,
                                    control=control.simulate.network(),
                                    statsonly=NULL,
                                    output=c("networkDynamic", "stats", "changes"),
                                    stats.form = FALSE,
                                    stats.diss = FALSE,
                                    verbose=FALSE, ...){
  if(nsim>1) stop("Simulating more than one chain of networks is not supported at this time. If you want to simulate over multiple time steps, use the time.slices argument.")

  if(!is.null(statsonly)){
    warning("Argument `statsonly' for STERGM simulate() is deprecated and may be removed in a future version. Use `output' instead.")
    output <- if(statsonly) "stats" else "networkDynamic"
  }

  nwend <- attr(object,"end")
  start <- if(is.null(nwend)) NVL(time.start,0) else{
    if(!is.null(time.start)){
      if(time.start!=nwend) warning("Argument time.start specified for a network that already has a time stamp. Overriding the time stamp.")
      time.start
    }else nwend
  }

  nw <- network.extract.with.lasttoggle(object, start)
  vActiveIDs <- nw %v% ".networkDynamicID"
  delete.vertex.attribute(nw, ".networkDynamicID")
 
  output <- match.arg(output)
  
  sim <- simulate.network(nw, nsim=1, seed=NULL,
                          formation=formation, dissolution=dissolution,
                          coef.form=coef.form,coef.diss=coef.diss,
                          constraints = constraints,
                          monitor = monitor,
                          time.slices=time.slices, time.burnin=time.burnin, time.interval=time.interval,
                          control=control,          
                          output=switch(output, networkDynamic = "changes", output),
                          stats.form = stats.form,
                          stats.diss = stats.diss,
                          verbose=verbose, ...)
  
  ## If all the user wants is statistics or a list of toggles, we are done.
  if(output!="networkDynamic") return(sim)

  if(verbose) cat("Updating networkDynamic.\n")
  
  ## Map the vertex IDs:
  sim[,"tail"] <- vActiveIDs[sim[,"tail"]]
  sim[,"head"] <- vActiveIDs[sim[,"head"]]

  object  <- networkDynamic.apply.changes(object, sim)

  attributes(object) <- c(attributes(object), # Don't clobber existing attributes!
                          list(formation = ergm.update.formula(formation,nw~.),
                               dissolution = ergm.update.formula(dissolution,nw~.),
                               stats.form = rbind(if(isTRUE(attr(sim,"formation")==attr(object,"formation"))) attr(object,"stats.form"),attr(sim,"stats.form")),
                               stats.diss = rbind(if(isTRUE(attr(sim,"dissolution")==attr(object,"dissolution"))) attr(object,"stats.diss"),attr(sim,"stats.diss")),
                               monitor = attr(sim,"monitor"),
                               stats = rbind(if(isTRUE(attr(sim,"monitor")==attr(object,"monitor"))) attr(object,"stats"),attr(sim,"stats")),
                               coef.form=coef.form,
                               coef.diss=coef.diss,
                               constraints=constraints,
                               start = attr(object,"start") + 0,
                               end = start + time.slices,
                               changes = rbind(attr(object,"changes"),matrix(c(sim), nrow=nrow(sim),ncol=ncol(sim),dimnames=list(rownames(sim),colnames(sim))))
                               ))
  object
}
