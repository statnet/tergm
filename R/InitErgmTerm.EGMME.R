#  File R/InitErgmTerm.EGMME.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################

stopifnot_dynamic <- function(nw, ..., dynamic=FALSE, .netseries.OK=FALSE){
  if(!dynamic && ! "lasttoggle" %in% list.network.attributes(nw)){
    msg <- paste0("This term requires either dynamic data (",
                  if(.netseries.OK) "network series or ",
                  "network dynamic or last toggle information) ",
                  "or must be set to dynamic mode (typically by passing ",
                  sQuote("dynamic=TRUE)"), " to ", sQuote("summary()"),
                  ", ", sQuote("simulate()"), ", etc..")
    ergm_Init_abort(msg)
  }
}

`InitErgmTerm.Form (dynamic)` <- function(nw, arglist,  ...) {
  stopifnot_dynamic(nw, .netseries.OK=TRUE, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  c(list(name="on_union_lt_net_Network",
         auxiliaries = ~.union.lt.net + .lasttoggle + .previous.lt.net,
         submodel = m,
         duration=TRUE),
    ergm_propagate_ext.encode(m),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("Form")))
}

InitErgmTerm..union.lt.net<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
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


`InitErgmTerm.Persist (dynamic)` <- function(nw, arglist,  ...) {
  stopifnot_dynamic(nw, .netseries.OK=TRUE, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  c(list(name="on_intersect_lt_net_Network",
         auxiliaries = ~.intersect.lt.net() + .lasttoggle,
         submodel = m,
         duration=TRUE),
    ergm_propagate_ext.encode(m),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("Persist")))
}

`InitErgmTerm.Diss (dynamic)` <- function(nw, arglist,  ...) {
  stopifnot_dynamic(nw, .netseries.OK=TRUE, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  c(list(name="negate_on_intersect_lt_net_Network",
         auxiliaries = ~.intersect.lt.net() + .lasttoggle,
         submodel = m,
         duration=TRUE),
    ergm_propagate_ext.encode(m),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("Diss")))
}

InitErgmTerm..intersect.lt.net<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
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


`InitErgmTerm.Change (dynamic)` <- function(nw, arglist,  ...) {
  stopifnot_dynamic(nw, .netseries.OK=TRUE, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)

  c(list(name="on_discord_lt_net_Network",
         auxiliaries = ~.discord.lt.net() + .lasttoggle + .previous.lt.net,
         submodel = m,
         duration=TRUE),
    ergm_propagate_ext.encode(m),
    wrap.ergm_model(m, nw, ergm_mk_std_op_namewrap("Change")))
}


InitErgmTerm..discord.lt.net<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_discord_lt_net_Network",
       coef.names=c(),
       auxiliaries = ~ .lasttoggle,
       duration=TRUE,
       dependence=FALSE)
}


InitErgmTerm..previous.lt.net<-function(nw, arglist, ...) {
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_previous_lt_net_Network",
       coef.names=c(),
       auxiliaries = ~ .lasttoggle,
       duration=TRUE,
       dependence=FALSE)
}

InitErgmTerm..lasttoggle <- function(nw, arglist, ...){
  stopifnot_dynamic(nw, ...)
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c(),
                      vartypes = c(),
                      defaultvalues = list(),
                      required = c())

  list(name="_lasttoggle",
       coef.names=c(),
       duration=TRUE,
       dependence=FALSE,
       ext.encode = function(el, nw0) list(time=as.integer(NVL(nw0 %n% "time", 0)), lasttoggle=as.integer(nw0 %n% "lasttoggle")),
       ext.decode = function(ext.state, el, nw0){
         nw0 %n% "time" <- ext.state$time
         nw0 %n% "lasttoggle" <- matrix(ext.state$lasttoggle, ncol=3)
         list(el=el, nw0=nw0)
       }
       )
}
