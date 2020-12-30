#  File R/stergm.EGMME.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################
################################################################################
# The <stergm> function fits stergms from a specified formation and dissolution
# formula returning approximate MLE's based on MCMC estimation.
#
# --PARAMETERS--
#   formation   : the formation formula, as 'nw ~ term(s)'
#   dissolution : the dissolution formula, as 'nw ~ term(s)'
#   theta.form0 : the intial theta formation coefficients, or optionally if
#                 these are to be estimates, the string "MPLE";
#                 default="MPLE"
#   theta.diss  : the initial theta dissolution coefficients
#   seed        : an integer starting value for the random number generator;
#                 default=NULL
#   MH.burnin   : the number of proposals used in each MCMC step; this is ignored
#                 unless 'control$main.method'="Robbins-Monro"; any other style or
#                 the default style will not recognize this parameter;
#                 default=1000
#   constraints : a one-sided formula of the constraint terms; options are
#                      bd        degrees        nodegrees
#                      edges     degreedist     indegreedist
#                      observed  outdegreedist
#                default="~ ."; these may not work currently.
#   target.stats   :  a vector of the mean value parameters;
#                  default=the observed statistic from the 'nw' in formula
#   control     :  a list of control parameters returned from <control.stergm>;
#                  default=<control.stergm>()
#   verbose     :  whether ergm should be verbose (T or F); default=FALSE
#
# --RETURNED--
#   because a stergm object is the return type of several functions, and
#   because this is a rather lengthy list, and because the returned items
#   of this function borrow from the other stergm.* functions, this list 
#   provides the returned items for all funtions returning a stergm.
#   The symbol preceding each component indicates which function returns it,
#   but remember that, <stergm> will additionally return the items from
#   one of the other stergm functions as well:
#       <stergm>             = $
#       <stergm.RM>          = @
#
#   the components include:
#
#     @   coef.form   : the estimated formation model coefficients
#     @   coef.diss   : the estimated dissolution model coefficients
#      &  eta         : the estimated formation ?? coefficients
#   $     offset      : a logical vector whose ith entry tells whether the
#                        ith curved theta coeffient was offset/fixed
#   $     etamap      :  the list constituting the theta->eta mapping for the
#                        formation model; for details of its components,
#                        see <ergm.etamap>
#   $     MH.burnin   :  the number of proposals made in each MCMC step
#   $     formation   : the formation formula, as 'nw ~ term(s)'
#   $     dissolution : the dissolution formula, as 'nw ~ term(s)'
#   $     constraints : the constraints formula
#     @&  newnetwork  :  the final network sampled
#     @&  network    :  the 'nw' inputted to <ergm> via the 'formula'
#   $     prop.args.form     :  the MHP formation arguments passed to the
#                               InitErgmProposal rountines
#   $     prop.args.diss     :  the MHP dissolution arguments passed to the
#                               InitErgmProposal rountines
#   $     prop.weights.form  :  the method used to allocate probabilities of
#                               being proposed to dyads in the formation stage,
#                               as "TNT", "random", "nonobserved", or "default"
#      &  theta.original     :  the theta values at the start of the MCMC 
#                               sampling
#     @   theta.form.original:  the formation theta values at the start of the
#                               MCMC sampling
#   $     prop.weights.diss  :  as 'prop.weights.form', but for the dissolution
#                               model
#
################################################################################

stergm.EGMME <- function(nw, formation, dissolution, constraints, offset.coef.form, offset.coef.diss,
                   targets, target.stats, estimate,
                 control,
                 verbose) {

  if(!is.network(nw)) stop("Argument nw must be a network.")

#  if(is.null(nw %n% "lasttoggle")) nw %n% "lasttoggle" <- rep(round(-.Machine$integer.max/2), network.dyadcount(nw))
	
  if(is.null(nw %n% "time")) nw %n% "time" <- 0

  # EGMME requires targets, or there will be an error
  if (is.null(targets)) stop('stergm.EGMME requires targets parameter be specified')
 
  # Allow the user to specify targets as copied from formation or dissolution formula.
  if(is.character(targets)){
    targets <- switch(targets,
                      formation = formation,
                      dissolution = dissolution)
  }
  
  if(length(targets)==3){
    warning("Targets formula has an LHS, which will be ignored in favor of nw.")
    targets <- targets[c(1,3)] # in a formula f<-y~x, f[1]=~, f[2]=y, and f[3]=x
  }

  targets <- nonsimp_update.formula(targets,nw~., from.new="nw")
  formation <- nonsimp_update.formula(formation,nw~., from.new="nw")
  dissolution <- nonsimp_update.formula(dissolution,nw~., from.new="nw")

  # target formula should not have offsets. removing them
  if (any(offset.info.formula(targets)$term)) {
    message("Targets formula should not contain offset terms;
                they have been been removed.")
    targets <- filter_rhs.formula(targets, function(x) (if(is.call(x)) x[[1]] else x)!="offset")
  }

  control.transfer <- list(EGMME.MCMC.burnin.min="MCMC.burnin.min",
                           EGMME.MCMC.burnin.max="MCMC.burnin.max",
                           EGMME.MCMC.burnin.pval="MCMC.burnin.pval",
                           EGMME.MCMC.burnin.add="MCMC.burnin.add")
  for(arg in names(control.transfer))
      if(is.null(control[[control.transfer[[arg]]]]))
          control[control.transfer[[arg]]] <- list(control[[arg]])


  model.form <- ergm_model(formation, nw, expanded=TRUE, role="formation", term.options=control$term.options)
  model.diss <- ergm_model(dissolution, nw, expanded=TRUE, role="dissolution", term.options=control$term.options)
  model.mon <- ergm_model(targets, nw, expanded=TRUE, role="target", term.options=control$term.options)

  if(any(model.form$etamap$canonical==0) || any(model.diss$etamap$canonical==0) || any(model.mon$etamap$canonical==0)) stop("Equilibrium GMME for models based on curved ERGMs is not supported at this time.")

  p.free<-sum(!model.form$etamap$offsettheta)+sum(!model.diss$etamap$offsettheta)
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
        san(model.mon, basis=nw, target.stats=target.stats,
            constraints=constraints,
            control=control$SAN,
            only.last=TRUE,
            verbose=verbose)

    targets<-nonsimp_update.formula(targets,TARGET_STATS~., from.new="TARGET_STATS")
    formation <- nonsimp_update.formula(formation,TARGET_STATS~., from.new="TARGET_STATS")
    dissolution <- nonsimp_update.formula(dissolution,TARGET_STATS~., from.new="TARGET_STATS")
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

  if (verbose) message("Initializing Metropolis-Hastings proposals.")
  proposal.form <- ergm_proposal(constraints, weights=control$MCMC.prop.weights.form, control$MCMC.prop.args.form, nw, class="f")
  proposal.diss <- ergm_proposal(constraints, weights=control$MCMC.prop.weights.diss, control$MCMC.prop.args.diss, nw, class="d")

  if(!is.dyad.independent(proposal.form$arguments$constraints) || !is.dyad.independent(proposal.diss$arguments$constraints)){
    warning("Dyad-dependent constraint imposed. Note that the constraint is applied to the post-formation and post-dissolution networks y+ and y-, not the next time-step's network. This behavior may change in the future.")
  }

  # If some control$init is specified...
  
  if(!is.null(control$init.form)){
    # Check length of control$init.form.
    if(length(control$init.form)!=length(model.form$etamap$offsettheta)) {
      if(verbose) message("control$init.form is ", paste(control$init.form, collapse = " "), "\n", " number of statistics is ", length(model.form$coef.names))
      stop(paste("Invalid starting formation parameter vector control$init.form:",
                 "wrong number of parameters."))
    }
  }else control$init.form <- rep(NA, length(model.form$etamap$offsettheta)) # Set the default value of control$init.form.
  if(!is.null(offset.coef.form)) control$init.form[model.form$etamap$offsettheta]<-offset.coef.form
  names(control$init.form) <- model.form$coef.names

  if(!is.null(control$init.diss)){
    # Check length of control$init.diss.
    if(length(control$init.diss)!=length(model.diss$etamap$offsettheta)) {
      if(verbose) message("control$init.diss is ", paste(control$init.diss, collapse = " "), "\n", " number of statistics is ", length(model.diss$coef.names))
      stop(paste("Invalid starting dissolution parameter vector control$init.diss:",
                 "wrong number of parameters."))
    }
  }else control$init.diss <- rep(NA, length(model.diss$etamap$offsettheta)) # Set the default value of control$init.diss.  
  if(!is.null(offset.coef.diss)) control$init.diss[model.diss$etamap$offsettheta]<-offset.coef.diss
  names(control$init.diss) <- model.diss$coef.names

  initialfit <- stergm.EGMME.initialfit(formation, dissolution, targets, control$init.form, control$init.diss, nw, model.form, model.diss, model.mon, control, verbose)
  
  if(verbose) message("Fitting STERGM Equilibrium GMME.")

  if(control$parallel){
    ergm.getCluster(control, verbose=verbose)
    if(verbose && !is.null(ergm.getCluster(control))) message("Using parallel cluster.")
  }
  
  Cout <- switch(control$EGMME.main.method,
                 "Gradient-Descent" = stergm.EGMME.GD(initialfit$formation.fit$coef,
                   initialfit$dissolution.fit$coef, nw, model.form, model.diss, model.mon,
                   control=control, proposal.form=proposal.form,
                   proposal.diss=proposal.diss,
                  verbose),
                 stop("Method ", control$EGMME.main.method, " is not implemented.")
                )

  out <- list(network = nw, formation = formation, dissolution = dissolution, targets = targets, target.stats=model.mon$target.stats, estimate=estimate, covar = Cout$covar, opt.history=Cout$opt.history, sample=Cout$sample, sample.obs=NULL, control=control, reference = ~Bernoulli, mc.se = Cout$mc.se, constraints = constraints,
              formation.fit = with(Cout, list(network=nw, formula=formation, coef = eta.form, covar=covar.form, etamap = model.form$etamap, offset = model.form$etamap$offsettheta, constraints=constraints, estimate=estimate, control=control, reference = ~Bernoulli, mc.se = mc.se.form)),
              dissolution.fit = with(Cout, list(network=nw, formula=dissolution, coef = eta.diss, covar=covar.diss, etamap = model.diss$etamap, offset = model.diss$etamap$offsettheta, constraints=constraints, estimate=estimate, control=control, reference = ~Bernoulli, mc.se = mc.se.diss))
              )
  class(out$formation.fit)<-class(out$dissolution.fit)<-"ergm"
  
  out
}
