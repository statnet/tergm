stergm.CMLE.initialfit <- function(init.form, init.diss, nw0, nw.form, nw.diss, model.form, model.diss, control, verbose=FALSE){
  if(!is.null(control$init.method) && control$init.method == "zeros"){
    init.form[is.na(init.form)]<-0
    init.diss[is.na(init.diss)]<-0
    form.fit <- list(coef = init.form)
    diss.fit <- list(coef = init.diss)
  }else if(!any(is.na(init.form)) && !any(is.na(init.diss))){
    # Don't need to do anything.
    form.fit <- list(coef = init.form)
    diss.fit <- list(coef = init.diss)
  }else{
    # Fit formation MPLE:
    # If it's either in the initial network or missing in the formation network, it's fixed:
    Clist.fixed <- ergm.Cprepare(nw0 | is.na(nw.form), model.form) 
    Clist <- ergm.Cprepare(nw.form, model.form)
    form.fit <- ergm.mple(Clist, Clist.fixed, model.form, init=init.form, control=control$CMLE.control.form, verbose=verbose)
    
    # Fit dissolution MPLE:
    # If it's either absent in the initial network or missing in the dissolution network, it's fixed:
    Clist.fixed <- ergm.Cprepare((!nw0) | is.na(nw.diss), model.diss) 
    Clist <- ergm.Cprepare(nw.diss, model.diss)
    diss.fit <- ergm.mple(Clist, Clist.fixed, model.diss, init=init.diss, control=control$CMLE.control.diss, verbose=verbose)    
  }
  
  return(list(formation.fit=c(form.fit, offset=model.form$etamap$offsettheta, network = nw.form, reference="Bernoulli", estimate="MPLE", control=control$CMLE.control.form),
              dissolution.fit=c(diss.fit, offset=model.diss$etamap$offsettheta, network = nw.diss, reference="Bernoulli", estimate="MPLE", control=control$CMLE.control.diss)))
}


