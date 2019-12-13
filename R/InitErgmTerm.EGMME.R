InitErgmTerm.FormE <- function(nw, arglist, response=NULL,  ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  # get the network and formula
  f <- a$formula

  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)
  
  m <- ergm_model(f, nw, response=response,...)
  inputs <- to_ergm_Cdouble(m)

  nw0 <- nw
  nw0[,] <- 0 # Delete edges but not lasttoggles.
  gs <- summary(m, nw0)
  
  c(list(name="on_union_lt_net_Network",
         coef.names = paste0("Form",'(',param_names(m, canonical=TRUE),')'),
         auxiliaries = ~.union.lt.net(),
         inputs=inputs,
         emptynwstats=gs,
         dependence=!is.dyad.independent(m),
         duration=TRUE),
    passthrough.curved.ergm_model(m, function(x) paste0('Form(',x,')')))
}

InitErgmTerm..union.lt.net<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_union_lt_net_Network",
       coef.names=c(),
       duration=TRUE,
       dependence=FALSE)
}


InitErgmTerm.DissE <- function(nw, arglist, response=NULL,  ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

    # get the network and formula
  f <- a$formula

  if(length(f)==2) f <- nonsimp_update.formula(f, nw~.)
  else nw <- ergm.getnetwork(f)
  
  m <- ergm_model(f, nw, response=response,...)
  inputs <- to_ergm_Cdouble(m)
  
  gs <- summary(m)
  
  c(list(name="on_intersect_lt_net_Network",
         coef.names = paste0("Diss",'(',param_names(m, canonical=TRUE),')'),
         auxiliaries = ~.intersect.lt.net(),
         inputs=inputs,
         emptynwstats=gs,
         dependence=!is.dyad.independent(m),
         duration=TRUE),
    passthrough.curved.ergm_model(m, function(x) paste0('Diss(',x,')')))
}

InitErgmTerm..intersect.lt.net<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_intersect_lt_net_Network",
       coef.names=c(),
       duration=TRUE,
       dependence=FALSE)
}

InitErgmTerm..lasttoggle <- function(nw, arglist, ...){
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_lasttoggle",
       coef.names=c(),
       duration=TRUE,
       dependence=FALSE,
       ext.encode = function(el, nw0) list(time=as.integer(nw0 %n% "time"), lasttoggle=as.integer(nw0 %n% "lasttoggle")),
       ext.decode = function(ext.state, el, nw0){
         nw0 %n% "time" <- ext.state$time
         nw0 %n% "lasttoggle" <- matrix(ext.state$lasttoggle, ncol=3)
         list(el=el, nw0=nw0)
       }
       )
}
