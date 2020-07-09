#  File R/tergm.EGMME.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2019 Statnet Commons
#######################################################################

tergm.EGMME <- function(formula, constraints, offset.coef,
                   targets, target.stats, SAN.offsets, estimate,
                 control,
                 verbose) {
  nw <- eval_lhs.formula(formula)

  if(!is.network(nw)) stop("Argument nw must be a network.")

  if(is.null(nw %n% "time")) nw %n% "time" <- 0

  # EGMME requires targets, or there will be an error
  if (is.null(targets)) stop('tergm.EGMME requires targets parameter be specified')
  
  if(is.character(targets)) {    
    targets <- switch(targets,
                      formation = .extract.fd.formulae(formula)$form,
                      dissolution = .extract.fd.formulae(formula)$diss)
  }
  
  if(length(targets)==3) {
    warning("Targets formula has an LHS, which will be ignored in favor of nw.")
    targets <- targets[c(1,3)] # in a formula f<-y~x, f[1]=~, f[2]=y, and f[3]=x
  }

  targets <- nonsimp_update.formula(targets,nw~., from.new="nw")
  formula <- nonsimp_update.formula(formula,nw~., from.new="nw")
  SAN.formula <- targets # including any offsets

  if (any(ergm_model(targets, nw)$etamap$offset)) {
    message("Targets contains offset terms;
                they will only be used during the SAN run.")
    targets <- statnet.common::filter_rhs.formula(targets, function(x) !inherits(x, "call") || !(x[[1]] == "offset"))
  }

  control.transfer <- list(EGMME.MCMC.burnin.min="MCMC.burnin.min",
                           EGMME.MCMC.burnin.max="MCMC.burnin.max",
                           EGMME.MCMC.burnin.pval="MCMC.burnin.pval",
                           EGMME.MCMC.burnin.add="MCMC.burnin.add")
  for(arg in names(control.transfer))
      if(is.null(control[[control.transfer[[arg]]]]))
          control[control.transfer[[arg]]] <- list(control[[arg]])

  if (verbose) message("Initializing Metropolis-Hastings proposal.")
  proposal <- ergm_proposal(constraints, weights=control$MCMC.prop.weights, control$MCMC.prop.args, nw, class="t")
  proposal.SAN <- ergm_proposal(constraints, weights=control$SAN.control$SAN.prop.weights, control$SAN.control$SAN.prop.args, nw, class="c")
  
  model <- ergm_model(formula, nw, term.options=control$term.options, extra.aux=list(proposal=proposal$auxiliaries, system=~.lasttoggle))
  model.SAN <- ergm_model(SAN.formula, nw, term.options=control$SAN.control$term.options, extra.aux=list(proposal=proposal.SAN$auxiliaries))  
  
  proposal$aux.slots <- model$slots.extra.aux$proposal
  proposal.SAN$aux.slots <- model.SAN$slots.extra.aux$proposal
  
  model.mon <- ergm_model(targets, nw, term.options=control$term.options)

  if(any(model$etamap$canonical==0) || any(model.mon$etamap$canonical==0) || any(model.SAN$etamap$canonical==0)) stop("Equilibrium GMME for models based on curved ERGMs is not supported at this time.")

  p.free<-sum(!model$etamap$offsettheta)
  if(p.free==0) stop("Model specification has no free parameters (all are offsets).")
  q<-length(model.mon$etamap$offsettheta)
  if(p.free>q) stop("Fitting ",p.free," free parameters on ",q," target statistics. The specification is underidentified.")

  nw.stats<-summary(model.mon, nw=nw)
  if(!is.null(target.stats)){
    if(length(nw.stats)!=length(target.stats))
      stop("Incorrect length of the target.stats vector: should be ", length(nw.stats), " but is ",length(target.stats),".")
        
    if(verbose) message("Constructing an approximate response network.")
    ## If target.stats are given, overwrite the given network and targets
    ## with SAN-ed network and targets.
    
    nw <- TARGET_STATS <-
        san(model.SAN, basis=nw, target.stats=target.stats,
            constraints=proposal.SAN,
            control=control$SAN.control,
            only.last=TRUE,
            verbose=verbose,
            offset.coef=SAN.offsets)

    targets<-nonsimp_update.formula(targets,TARGET_STATS~., from.new="TARGET_STATS")
    formula <- nonsimp_update.formula(formula,TARGET_STATS~., from.new="TARGET_STATS")
    nw.stats <- summary(model.mon, nw)

    if(verbose){
      message("SAN summary statistics:")
      message_print(nw.stats)
      message("Meanstats Goal:")
      message_print(target.stats)
      message("Difference: SAN target.stats - Goal target.stats =")
      message_print(round(nw.stats-target.stats,0))
    }
  }

  model.mon$nw.stats <- nw.stats
  model.mon$target.stats <- if(!is.null(target.stats)) vector.namesmatch(target.stats, names(model.mon$nw.stats)) else model.mon$nw.stats

  # If some control$init is specified...
  
  if(!is.null(control$init)){
    # Check length of control$init.
    if(length(control$init)!=length(model$etamap$offsettheta)) {
      if(verbose) message("control$init is ", paste(control$init, collapse = " "), "\n", " number of statistics is ", length(model$coef.names))
      stop(paste("Invalid starting formation parameter vector control$init:",
                 "wrong number of parameters."))
    }
  }else control$init <- rep(NA, length(model$etamap$offsettheta)) # Set the default value of control$init.
  if(!is.null(offset.coef)) control$init[model$etamap$offsettheta]<-offset.coef
  names(control$init) <- model$coef.names

  initialfit <- tergm.EGMME.initialfit(control$init, nw, model, formula, model.mon, targets, control, verbose)
  
  if(verbose) message("Fitting TERGM Equilibrium GMME.")

  if(control$parallel){
    ergm.getCluster(control, verbose=verbose)
    if(verbose && !is.null(ergm.getCluster(control))) message("Using parallel cluster.")
  }
  
  Cout <- switch(control$EGMME.main.method,
                 "Gradient-Descent" = tergm.EGMME.GD(initialfit$coef,
                   nw, model, model.mon,
                   control=control, proposal=proposal,
                  verbose),
                 stop("Method ", control$EGMME.main.method, " is not implemented.")
                )

  out <- list(network = nw,
              formula = formula,
              coef = Cout$eta,              
              targets = targets, 
              target.stats=model.mon$target.stats, 
              estimate=estimate, 
              covar = Cout$covar, 
              opt.history=Cout$opt.history, 
              sample=Cout$sample, 
              sample.obs=NULL, 
              control=control, 
              reference = ~Bernoulli, 
              mc.se = Cout$mc.se, 
              constraints = constraints,
              fit = with(Cout, list(network=nw, 
                                    formula=formula, 
                                    coef = eta, 
                                    covar=covar, 
                                    etamap = model$etamap, 
                                    offset = model$etamap$offsettheta, 
                                    constraints=constraints, 
                                    estimate=estimate, 
                                    control=control, 
                                    reference = ~Bernoulli, 
                                    mc.se = mc.se, 
                                    ergm_version = packageVersion("ergm"))))
  class(out$fit)<-"ergm"
  
  out
}
