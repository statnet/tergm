stergm.CMLE <- function(nw, formation, dissolution, constraints, times, offset.coef.form, offset.coef.diss,
                        eval.loglik,
                        estimate,
                        control,
                        verbose) {

  if(is.null(times)){
    if(inherits(nw, "network.list") || is.list(nw)){
      times  <- c(1,2)
      warning("Time points not specified for a list. Modeling transition from the first to the second network. This behavior may change in the future.")
    }else if(inherits(nw,"networkDynamic")){
      times  <- c(0,1)
      warning("Time points not specified for a networkDynamic. Modeling transition from time 0 to 1.")
    }
  }
  
  if(length(times)<2) stop("Time points whose transition is to be modeled was not specified.")
  if(length(times)>2) stop("Only two time points (one transition) are supported at this time.")

  if(inherits(nw, "network.list") || is.list(nw)){
    y0 <- nw[[times[1]]]
    y1 <- nw[[times[2]]]
  }else if(inherits(nw,"networkDynamic")){
    require(networkDynamic) # This is needed for the "%t%.network" function
    y0 <- nw %t% times[1]
    y1 <- nw %t% times[2]
  }

  if(!is.network(y0) || !is.network(y1)) stop("nw must be a networkDynamic, a network.list, or a list of networks.")
  
  
  # Construct the formation and dissolution networks; the
  # network.update cannot be used to copy attributes from y0 to y.form and
  # y.diss, since NA edges will be lost.
  
  y.form <- y0 | y1
  y.form <- nvattr.copy.network(y.form, y0)
  formation <- ergm.update.formula(formation, y.form~.)

  y.diss <- y0 & y1
  y.diss <- nvattr.copy.network(y.diss, y0)
  dissolution <- ergm.update.formula(dissolution, y.diss~.)

  # Construct new constraints

  constraints.form <- if(constraints==~.) ~atleast(y0) else ergm.update.formula(constraints, ~.+atleast(y0))
  constraints.diss <- if(constraints==~.) ~atmost(y0) else ergm.update.formula(constraints, ~.+atmost(y0))
  
  # Apply initial values passed to control.stergm() the separate controls, if necessary.
  if(is.null(control$CMLE.control.form$init)) control$CMLE.control.form$init <- control$init.form
  if(is.null(control$CMLE.control.diss$init)) control$CMLE.control.diss$init <- control$init.diss
 
  model.form<-ergm.getmodel(formation, y.form, initialfit=TRUE)
  model.diss<-ergm.getmodel(dissolution, y.diss, initialfit=TRUE)

  if(!is.null(control$CMLE.control.form$init)){
    # Check length of control$CMLE.control.form$init.
    if(length(control$CMLE.control.form$init)!=length(model.form$etamap$offsettheta)) {
      if(verbose) cat("control$CMLE.control.form$init is", control$CMLE.control.form$init, "\n", "number of statistics is",length(model.form$coef.names), "\n")
      stop(paste("Invalid starting formation parameter vector control$CMLE.control.form$init:",
                 "wrong number of parameters."))
    }
  }else control$CMLE.control.form$init <- rep(NA, length(model.form$etamap$offsettheta)) # Set the default value of control$CMLE.control.form$init.
  if(!is.null(offset.coef.form)) control$CMLE.control.form$init[model.form$etamap$offsettheta]<-offset.coef.form
  names(control$CMLE.control.form$init) <- model.form$coef.names

  if(!is.null(control$CMLE.control.diss$init)){
    # Check length of control$CMLE.control.diss$init.
    if(length(control$CMLE.control.diss$init)!=length(model.diss$etamap$offsettheta)) {
      if(verbose) cat("control$CMLE.control.diss$init is", control$CMLE.control.diss$init, "\n", "number of statistics is",length(model.diss$coef.names), "\n")
      stop(paste("Invalid starting dissolution parameter vector control$CMLE.control.diss$init:",
                 "wrong number of parameters."))
    }
  }else control$CMLE.control.diss$init <- rep(NA, length(model.diss$etamap$offsettheta)) # Set the default value of control$CMLE.control.diss$init.  
  if(!is.null(offset.coef.diss)) control$CMLE.control.diss$init[model.diss$etamap$offsettheta]<-offset.coef.diss
  names(control$CMLE.control.diss$init) <- model.diss$coef.names


  # Translate the "estimate" from the stergm() argument to the ergm() argument.
  ergm.estimate <- switch(estimate,
                     CMLE = "MLE",
                     CMPLE = "MPLE")
  
  # Now, call the ergm()s:
  cat("Fitting formation...\n")
  fit.form <- ergm(formation, constraints=constraints.form, offset.coef=offset.coef.form, eval.loglik=eval.loglik, estimate=ergm.estimate, control=control$CMLE.control.form, verbose=verbose)
  cat("Fitting dissolution...\n")
  fit.diss <- ergm(dissolution, constraints=constraints.diss, offset.coef=offset.coef.diss, eval.loglik=eval.loglik, estimate=ergm.estimate, control=control$CMLE.control.diss, verbose=verbose)

  # Construct the output list. Conveniently, this is mainly a list consisting of two ergms.
  
  list(network=nw, times=times, formation=formation, dissolution=dissolution, formation.fit=fit.form, dissolution.fit=fit.diss, estimate=estimate)
}
