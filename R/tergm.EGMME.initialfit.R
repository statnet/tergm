#  File R/tergm.EGMME.initialfit.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2019 Statnet Commons
#######################################################################
tergm.EGMME.initialfit<-function(init, nw, model, formula, model.mon, formula.mon, control, verbose=FALSE){

  # Remove offset() from coefficient names.
  .do <- function(x) sub('offset\\((.+)\\)', '\\1', x)
  
  if(!is.null(control$init.method) && control$init.method == "zeros") {
    init[is.na(init)] <- 0
  } else if(!any(is.na(init))) {
    # Don't need to do anything.
  } else {
    fd.formulae <- .extract.fd.formulae(formula)
    
    form <- fd.formulae$form
    diss <- fd.formulae$diss
    
    model.form <- ergm_model(form, nw=nw, term.options=control$term.options)
    model.diss <- ergm_model(diss, nw=nw, term.options=control$term.options)
    
    wf <- which(grepl("Form\\(.+\\)", model$coef.names))
    wd <- which(grepl("Diss\\(.+\\)", model$coef.names))
    
    init.form <- init[wf]
    init.diss <- init[wd]
    
    names(init.form) <- model.form$coef.names
    names(init.diss) <- model.diss$coef.names
    
    form <- nonsimp_update.formula(form, nw ~ ., from.new = "nw")
    model.form$formula <- form
  
    if(all(model.form$coef.names[!model.form$etamap$offsettheta] %in% .do(model.mon$coef.names))
             && (
                  all(model.diss$etamap$offsettheta)
                  || (
                       length(model.diss$coef.names) == 1
                       && .do(model.diss$coef.names) == "edges"
                       && "mean.age" %in% .do(model.mon$coef.names)
                     )
                )
             && all(.do(model.diss$coef.names) %in% model.form$coef.names)
             && is.dyad.independent(model.diss)) {
      if(verbose) message("Formation statistics are analogous to targeted statistics, dissolution is fixed or is edges with a mean.age target, dissolution terms appear to have formation analogs, and dissolution process is dyad-independent, so using edges dissolution approximation  (Carnegie et al.).")
      
      if(!all(model.diss$etamap$offsettheta)){ # This must mean that the two provisos above are satisfied.
        mean.age <- model.mon$target.stats[.do(model.mon$coef.names)=="mean.age"]
        init.diss <- log(mean.age+1)
        names(init.diss) <- "edges"
      }
      
      form.offset.coef <- if(any(model.form$etamap$offsetmap)) init.form[model.form$etamap$offsetmap] else NULL
      
      # Fit an ERGM to the formation terms:
      form.targets <- model.mon$target.stats[match(model.form$coef.names,.do(model.mon$coef.names))]
      form.targets <- form.targets[!model.form$etamap$offsettheta]
      
      ergm_control <- NVL(control$EGMME.initialfit.control, control.ergm())
      ergm_control$init <- init.form
      
      init.form<-coef(ergm(model.form$formula,control=ergm_control, target.stats=form.targets, offset.coef = form.offset.coef, eval.loglik=FALSE))
      # Now, match up non-offset formation terms with dissolution terms.
      # In case it's not obvious (it's not to me) what the following
      # does, it takes non-offset elements of init.form, then, from
      # those, it takes those elements that correspond (by name) to the
      # dissolution coefficient names and decrements them by init.diss.
      #
      # Yes, I am also a little surprised that assigning to a
      # double-index works.
      init.form[!model.form$etamap$offsettheta][match(.do(names(init.diss)),names(init.form[!model.form$etamap$offsettheta]))] <-
        init.form[!model.form$etamap$offsettheta][match(.do(names(init.diss)),names(init.form[!model.form$etamap$offsettheta]))] - init.diss
        
      init[wd] <- init.diss
      init[wf] <- init.form
    }else{
      stop("No initial parameter method for specified model and targets combination is implemented. Specify via control$init.")
    }
  }
  out <- list(formula=formula, coef = init, targets = formula.mon, target.stats=model.mon$target.stats, nw = nw, control = control, fit = list(coef=init, etamap = model$etamap))
  class(out)<-"tergm"
  out
}
