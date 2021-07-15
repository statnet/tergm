#  File R/tergm.EGMME.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2021 Statnet Commons
################################################################################

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
                      dissolution = .extract.fd.formulae(formula)$pers)
  }

  if(length(targets)==3) {
    warning("Targets formula has an LHS, which will be ignored in favor of nw.")
    targets <- targets[c(1,3)] # in a formula f<-y~x, f[1]=~, f[2]=y, and f[3]=x
  }

  targets <- nonsimp_update.formula(targets,nw~., from.new="nw")
  formula <- nonsimp_update.formula(formula,nw~., from.new="nw")
  SAN.formula <- targets # including any offsets

  target_model <- ergm_model(targets, nw, dynamic=TRUE, term.options = control$term.options)
  if(any(target_model$etamap$offsetmap)) {
    message("Targets contains offset statistics; they will only be used during the SAN run, and removal of the offset statistics will be attempted for the EGMME targets.")

    non_offsets <- !target_model$etamap$offsetmap
    targets <- trim_env(~.SubsetStatistics(targets, non_offsets), keep = c("targets", "non_offsets"))
    
    updated_target_model <- ergm_model(targets, nw, dynamic=TRUE, term.options = control$term.options)
    if(any(updated_target_model$etamap$offsetmap) || 
       sum(!target_model$etamap$offsetmap) != nparam(updated_target_model, canonical = TRUE)) {
         stop("Failed to remove offsets from targets formula; please specify targets formula without offsets.")
    }
  }

  control.transfer <- list(EGMME.MCMC.burnin.min="MCMC.burnin.min",
                           EGMME.MCMC.burnin.max="MCMC.burnin.max",
                           EGMME.MCMC.burnin.pval="MCMC.burnin.pval",
                           EGMME.MCMC.burnin.add="MCMC.burnin.add")
  for(arg in names(control.transfer))
      if(is.null(control[[control.transfer[[arg]]]]))
          control[control.transfer[[arg]]] <- list(control[[arg]])

  if (verbose) message("Initializing Metropolis-Hastings proposal.")
  proposal <- ergm_proposal(constraints, weights = control$MCMC.prop.weights, arguments = control$MCMC.prop.args, nw = nw, hints = control$MCMC.prop, class="t")
  proposal.SAN <- ergm_proposal(constraints, weights = control$SAN$SAN.prop.weights, arguments = control$SAN$SAN.prop.args, nw = nw, hints = control$SAN$SAN.prop, class="c")
  
  model <- ergm_model(formula, nw, dynamic=TRUE, term.options=control$term.options, extra.aux=list(proposal=proposal$auxiliaries, system=~.lasttoggle))
  model.SAN <- ergm_model(SAN.formula, nw, dynamic=TRUE, term.options=control$SAN$term.options, extra.aux=list(proposal=proposal.SAN$auxiliaries))  
  
  proposal$aux.slots <- model$slots.extra.aux$proposal
  proposal.SAN$aux.slots <- model.SAN$slots.extra.aux$proposal
  
  model.mon <- ergm_model(targets, nw, dynamic=TRUE, term.options=control$term.options)

  if(is.curved(model) || is.curved(model.mon) || (control$SAN$SAN.maxit > 0 && is.curved(model.SAN))) stop("Equilibrium GMME for models based on curved ERGMs is not supported at this time.")

  p.free<-sum(!model$etamap$offsettheta)
  if(p.free==0) stop("Model specification has no free parameters (all are offsets).")
  q<-length(model.mon$etamap$offsettheta)
  if(p.free>q) stop("Fitting ",p.free," free parameters on ",q," target statistics. The specification is underidentified.")

  nw.stats<-summary(model.mon, nw=nw, dynamic=TRUE)
  if(!is.null(target.stats)){
    if(length(nw.stats)!=length(target.stats))
      stop("Incorrect length of the target.stats vector: should be ", length(nw.stats), " but is ",length(target.stats),".")
        
    if(verbose) message("Constructing an approximate response network.")
    ## If target.stats are given, overwrite the given network and targets
    ## with SAN-ed network and targets.
    
    if(control$SAN$SAN.maxit > 0) {
      if(sum(model.SAN$etamap$offsettheta) != length(SAN.offsets)) {
        stop("Incorrect number of offset coefficients specified for SAN: expected ", sum(model.SAN$etamap$offsettheta), "; got ", length(SAN.offsets), ".");
      }
      
      nw <- san(model.SAN, 
                basis = nw, 
                target.stats = target.stats,
                constraints = proposal.SAN,
                control = control$SAN,
                only.last = TRUE,
                verbose = verbose,
                offset.coef = SAN.offsets)
    }
    
    TARGET_STATS <- nw
    
    targets<-nonsimp_update.formula(targets,TARGET_STATS~., from.new="TARGET_STATS")
    formula <- nonsimp_update.formula(formula,TARGET_STATS~., from.new="TARGET_STATS")
    nw.stats <- summary(model.mon, nw, dynamic=TRUE)

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
      if(verbose) message("control$init is ", paste(control$init, collapse = " "), "\n", " number of statistics is ", length(param_names(model, canonical = TRUE)))
      stop(paste("Invalid starting formation parameter vector control$init:",
                 "wrong number of parameters."))
    }
  }else control$init <- rep(NA, length(model$etamap$offsettheta)) # Set the default value of control$init.
  
  if(!is.null(offset.coef)) {
    if(length(control$init[model$etamap$offsettheta])!=length(offset.coef)) {
      stop("Invalid offset parameter vector offset.coef: ",
           "wrong number of parameters: expected ",
           length(control$init[model$etamap$offsettheta]),
           " got ", length(offset.coef), ".")
    }
    control$init[model$etamap$offsettheta] <- offset.coef
  }

  # Make sure any offset elements are given in control$init.
  if(any(is.na(control$init) & model$etamap$offsettheta)) stop("The model contains offset terms whose parameter values have not been specified:", paste.and(param_names(model)[is.na(control$init)&model$offsettheta]), ".", sep="")
    
  names(control$init) <- param_names(model)

  initialfit <- tergm.EGMME.initialfit(control$init, nw, model, formula, model.mon, targets, control, verbose)
  
  if(verbose) message("Fitting TERGM Equilibrium GMME.")

  if(control$parallel){
    ergm.getCluster(control, verbose=verbose)
    if(verbose && !is.null(ergm.getCluster(control))) message("Using parallel cluster.")
  }
  
  Cout <- switch(control$EGMME.main.method,
                 "Gradient-Descent" = tergm.EGMME.GD(coef(initialfit),
                   nw, model, model.mon,
                   control=control, proposal=proposal,
                  verbose),
                 stop("Method ", control$EGMME.main.method, " is not implemented.")
                )

#' @importFrom utils getS3method
  out <- c(Cout,
           list(network = nw,
                coef = Cout$eta,
                targets = targets,
                target.stats=model.mon$target.stats,
                estimate=estimate,
                sample.obs=NULL,
                control=control,
                reference = ~Bernoulli,
                constraints = constraints,
                etamap = model$etamap,
                offset = model$etamap$offsettheta,
                mle.lik = NA,
                MPLE_is_MLE = FALSE,
                ergm_version = packageVersion("ergm")))

  class(out) <- c("tergm_EGMME", "tergm", "ergm")
  out
}
