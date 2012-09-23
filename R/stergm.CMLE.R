stergm.CMLE <- function(nw, formation, dissolution, constraints, times, offset.coef.form, offset.coef.diss,
                        eval.loglik,
                        estimate,
                        control,
                        verbose) {

  if(is.null(times)){
    if(inherits(nw, "network.list") || is.list(nw)){
      times  <- seq_along(nw)
      warning("Time points not specified for a list. Modeling transition from the between successive networks jointly. This behavior may change in the future.")
    }else if(inherits(nw,"networkDynamic")){
      times  <- c(0,1)
      warning("Time points not specified for a networkDynamic. Modeling transition from time 0 to 1.")
    }
  }
  
  if(length(times)<2) stop("Time points whose transition is to be modeled was not specified.")  


  # Construct a list of "from" networks and a list of "to" networks.
   if(inherits(nw,"networkDynamic")){
    require(networkDynamic) # This is needed for the "%t%.network" function
    y0s <- lapply(times[-length(times)], function(t) nw %t% t)
    y1s <- lapply(times[-1], function(t) nw %t% t)
  }else if(inherits(nw, "network.list") || is.list(nw)){
    y0s <- nw[times[-length(times)]]
    y1s <- nw[times[-1]]
    
    if(!all(sapply(c(y0s,y1s),is.network))) stop("nw must be a networkDynamic, a network.list, or a list of networks.")
  }

  y0s <- lapply(y0s,standardize.network)
  y1s <- lapply(y1s,standardize.network)

  y0s.NA <- sapply(y0s, network.naedgecount)>0
  if(any(y0s.NA)){    
    y0s <- switch(control$CMLE.NA.impute,
                  stop = stop("Transitioned-from network(s) at time(s) ", paste.and(times[y0s.NA]), " has missing dyads. Fix or choose an imputation option via CMLE.NA.impute control parameter."),
                  previous = {
                    if(y0s.NA[1]) stop("Imputation option `previous' cannot impute dyads of the first network in the series.")
                    for(t in seq_along(y0s)[-1])
                      if(y0s.NA[t]){
                        # Workaround for a bug in network (Ticket #80 in Trac)
                        na.el <-as.edgelist(is.na(y0s[[t]]))
                        na.eids <- apply(na.el, 1, function(e) get.edgeIDs(y0s[[t]], e[1],e[2], na.omit=FALSE))
                        y0s[[t]] <- delete.edges(y0s[[t]], na.eids)
                        y0s[[t]][na.el] <- y0s[[t-1]][na.el]
                      }
                    y0s
                  },
                  majority = {
                    lapply(y0s, function(y){
                      # Workaround for a bug in network (Ticket #80 in Trac)
                      na.el <-as.edgelist(is.na(y))
                      na.eids <- apply(na.el, 1, function(e) get.edgeIDs(y, e[1],e[2], na.omit=FALSE))
                      impute <- if(network.edgecount(y,na.omit=TRUE)/network.dyadcount(y,na.omit=TRUE)>0.5) 1 else 0
                      y <- delete.edges(y, na.eids)
                      y[na.el] <- impute
                      y
                    })
                  },
                  `0` = {
                    lapply(y0s, function(y){
                      # Workaround for a bug in network (Ticket #80 in Trac)
                      na.el <-as.edgelist(is.na(y))
                      na.eids <- apply(na.el, 1, function(e) get.edgeIDs(y, e[1],e[2], na.omit=FALSE))
                      y <- delete.edges(y, na.eids)
                    })
                  },
                  `1` = {
                    lapply(y0s, function(y){
                      # Workaround for a bug in network (Ticket #80 in Trac)
                      na.el <-as.edgelist(is.na(y))
                      na.eids <- apply(na.el, 1, function(e) get.edgeIDs(y, e[1],e[2], na.omit=FALSE))
                      y <- delete.edges(y, na.eids)
                      y[na.el] <- 1
                      y
                    })
                  }
                  
                  )
  }

  if(length(times)>2){
    y0 <- combine.networks(y0s, standardized=TRUE,blockname=".stergm.CMLE.time.index")
    y1 <- combine.networks(y1s, standardized=TRUE,blockname=".stergm.CMLE.time.index")

    # Check that these networks can be combined for this model.
    bad.stat <-
      !((apply(rbind(sapply(y0s, function(nw) summary(ergm.update.formula(formation,nw~.)))),1,sum)==summary(ergm.update.formula(formation,y0~.))) &
        (apply(rbind(sapply(y0s, function(nw) summary(ergm.update.formula(dissolution,nw~.)))),1,sum)==summary(ergm.update.formula(dissolution,y0~.))) &
        (apply(rbind(sapply(y1s, function(nw) summary(ergm.update.formula(formation,nw~.)))),1,sum)==summary(ergm.update.formula(formation,y1~.))) &
        (apply(rbind(sapply(y1s, function(nw) summary(ergm.update.formula(dissolution,nw~.)))),1,sum)==summary(ergm.update.formula(dissolution,y1~.))))
    if(any(bad.stat)) stop("Fitting the terms ", paste.and(names(bad.stat)[bad.stat]), " over multiple network transitions is not supported at this time.")
  }else{
    # We are about to do logical operations on networks, so make sure
    # tail-head orderings match up.
    y0 <- standardize.network(y0s[[1]])
    y1 <- standardize.network(y1s[[1]])
  }
  
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
  if(length(times)>2) constraints.form <- ergm.update.formula(constraints.form, ~.+blockdiag(".stergm.CMLE.time.index"))
  # TODO: Some unlucky variable names can break this. We need to figure out a way around this.
  environment(constraints.form) <- environment()
  
  constraints.diss <- if(constraints==~.) ~atmost(y0) else ergm.update.formula(constraints, ~.+atmost(y0))
  if(length(times)>2) constraints.diss <- ergm.update.formula(constraints.diss, ~.+blockdiag(".stergm.CMLE.time.index"))
  # TODO: Some unlucky variable names can break this. We need to figure out a way around this.
  environment(constraints.diss) <- environment()
  
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
