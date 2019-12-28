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

  nw0 <- nw
  nw0[,] <- 0 # Delete edges but not lasttoggles.
  wm <- wrap.ergm_model(m, nw0, response=response, function(x) paste0('Form(',x,')'))
  ext.encode <- wm$ext.encode
  wm$ext.encode <-
    if(!is.null(ext.encode)){
      function(el, nw0)
        c(ext.encode(el, nw0), list(list(time=as.integer(nw0 %n% "time"), lasttoggle=as.integer(nw0 %n% "lasttoggle"))))
    }else{
      function(el, nw0)
        list(list(time=as.integer(nw0 %n% "time"), lasttoggle=as.integer(nw0 %n% "lasttoggle")))
    }

  ## TODO: Ideally, this term could grab the extended state from the
  ## auxiliary, avoiding duplication. Fortunately, these functions
  ## won't be called very often.
  wm$emptynwstats <- function(ext.state){
    mine <- ext.state[[length(ext.state)]]
    t <- mine$time
    lt <- matrix(mine$lasttoggle, ncol=3)
    changed <-  lt[,3]==t # If an edge just got toggled, then it must be in the formation network.   
    s <- ergm_state(nw0, el=lt[changed,1:2], model=m, ext.state=ext.state[-length(ext.state)])
    summary(s)
  }

  c(list(name="on_union_lt_net_Network",
         auxiliaries = ~.union.lt.net(),
         submodel = m,
         duration=TRUE),
    wm)
}

InitErgmTerm..union.lt.net<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_union_lt_net_Network",
       coef.names=c(),
       auxiliaries = ~ .lasttoggle,
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

  nw0 <- nw
  nw0[] <- FALSE
  c(list(name="on_intersect_lt_net_Network",
         auxiliaries = ~.intersect.lt.net(),
         submodel = m,
         duration=TRUE),
    wrap.ergm_model(m, nw0, response=response, function(x) paste0('Diss(',x,')')))
}

InitErgmTerm..intersect.lt.net<-function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_intersect_lt_net_Network",
       coef.names=c(),
       auxiliaries = ~ .lasttoggle,
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
