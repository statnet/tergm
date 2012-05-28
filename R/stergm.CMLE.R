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
#                               InitMHP rountines
#   $     prop.args.diss     :  the MHP dissolution arguments passed to the
#                               InitMHP rountines
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

stergm.CMLE <- function(nw, formation, dissolution, times, offset.coef.form, offset.coef.diss,
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

  CMPLE.is.CMLE <- (ergm.independencemodel(model.form)
                    && ergm.independencemodel(model.diss))

  
  # Get the initial fit:
  initialfit <- stergm.CMLE.initialfit(init.form=control$CMLE.control.form$init,
                                       init.diss=control$CMLE.control.diss$init,
                                       nw0=y0, nw.form=y.form, nw.diss=y.diss,
                                       model.form=model.form,
                                       model.diss=model.diss,
                                       control=control, verbose=verbose)

  if(estimate=="CMPLE" || (CMPLE.is.CMLE && !control$force.main)){
    initialfit <- c(initialfit, list(network=nw, times=times, estimate=estimate))
    initialfit$formation.fit$formula <- formation
    initialfit$dissolution.fit$formula <- dissolution
    initialfit$formation.fit$constraints <- ~atleast(y0)
    initialfit$dissolution.fit$constraints <- ~atmost(y0)
    initialfit$formation.fit$target.stats <- summary(formation)
    initialfit$dissolution.fit$target.stats <- summary(dissolution)
    class(initialfit$formation.fit) <- class(initialfit$dissolution.fit) <- "ergm"

    initialfit$estimate <- estimate
    initialfit$control<-control
    return(initialfit)
  }else{
    # Copy the initial parameters to ergm controls.
    control$CMLE.control.form$init <- initialfit$formation.fit$coef
    control$CMLE.control.diss$init <- initialfit$dissolution.fit$coef
  }

  # Translate the "estimate" from the stergm() argument to the ergm() argument.
  ergm.estimate <- switch(estimate,
                     CMLE = "MLE",
                     CMPLE = "MPLE")
  
  # Now, call the ergm()s:
  cat("Fitting formation:\n")
  fit.form <- ergm(formation, constraints=~atleast(y0), offset.coef=offset.coef.form, eval.loglik=eval.loglik, estimate=ergm.estimate, control=control$CMLE.control.form, verbose=verbose)
  cat("Fitting dissolution:\n")
  fit.diss <- ergm(dissolution, constraints=~atmost(y0), offset.coef=offset.coef.diss, eval.loglik=eval.loglik, estimate=ergm.estimate, control=control$CMLE.control.diss, verbose=verbose)

  # Construct the output list. Conveniently, this is mainly a list consisting of two ergms.
  
  list(network=nw, times=times, formation=formation, dissolution=dissolution, formation.fit=fit.form, dissolution.fit=fit.diss, estimate=estimate, initialfit = initialfit)
}
