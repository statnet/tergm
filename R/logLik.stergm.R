#  File R/logLik.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2018 Statnet Commons
#######################################################################



#' A \code{\link{logLik}} method for \code{\link[=ergm.object]{stergm}}.
#' 
#' Functions to return the log-likelihood associated with a
#' [`stergm`] CMLE fit, evaluating it if necessary.
#' See \code{\link[ergm]{logLik.ergm}} documentation for details and caveats.
#' 
#' If the log-likelihood was not computed for \code{object}, produces an error
#' unless \code{eval.loglik=TRUE}
#'
#' @param object A [`stergm`] fit, returned by
#' \code{\link{stergm}}, for \code{estimate="CMLE"}.
#' @param add Logical: If TRUE, instead of returning the log-likelihood, return
#' \code{object} with log-likelihood value set.
#' @param force.reeval Logical: If TRUE, reestimate the log-likelihood even if
#' \code{object} already has an estiamte.
#' @param eval.loglik Logical: If TRUE, evaluate the log-likelihood if not set
#' on \code{object}.
#' @param control A list of control parameters for algorithm tuning.
#' Constructed using \code{\link{control.logLik.ergm}}.
#' @param \dots Other arguments to the likelihood functions.
#' @return For \code{logLik.stergm}, \code{add=FALSE} (the default), a
#' \code{\link{logLik}} object. If \code{add=TRUE} (the default), an
#' \code{\link[=ergm]{ergm}} object or a
#' \code{\link[=stergm]{stergm}} object with the log-likelihood set. For
#' \code{logLikNull.stergm}, a \code{\link{logLik}} object.
#' @seealso \code{\link{logLik}}, \code{\link{ergm.bridge.llr}},
#' \code{\link{ergm.bridge.dindstart.llk}}
#' @references Hunter, D. R. and Handcock, M. S. (2006) \emph{Inference in
#' curved exponential family models for networks}, Journal of Computational and
#' Graphical Statistics.
#' @keywords models
#' @export
logLik.stergm<-function(object, add=FALSE, force.reeval=FALSE, eval.loglik=add || force.reeval, control=control.logLik.stergm(), ...){
  check.control.class("logLik.stergm","logLik.stergm")
  if(object$estimate=="EGMME") stop("Log-likelihood for ",object$estimate," is not meaningful.")
  
  if(!is.null(control$seed))  set.seed(as.integer(control$seed))

  if(add){
    object$formation.fit <- logLik(object$formation.fit, add=add, force.reeval=force.reeval, eval.loglik = add || force.reeval, control=control$control.form)
    object$dissolution.fit <- logLik(object$dissolution.fit, add=add, force.reeval=force.reeval, eval.loglik = add || force.reeval, control=control$control.diss)
    
    object
  }else{
    llk.form <- logLik(object$formation.fit, add=add, force.reeval=force.reeval, eval.loglik = eval.loglik, control=control$control.form)
    llk.diss <- logLik(object$dissolution.fit, add=add, force.reeval=force.reeval, eval.loglik = eval.loglik, control=control$control.diss)

    llk <- llk.form + llk.diss
    class(llk) <- "logLik"
    attr(llk,"df") <- attr(llk.form,"df") + attr(llk.diss,"df")
    attr(llk,"nobs") <- attr(llk.form,"nobs") + attr(llk.diss,"nobs")
    
    llk
  }
}

#' @rdname logLik.stergm
#' @description `logLikNull` method computes the null model
#'   likelihood. See [ergm::logLikNull()].
#' @importFrom ergm logLikNull
#' @export
logLikNull.stergm <- function(object, control=control.logLik.stergm(), ...){
    llk.form <- logLikNull(object$formation.fit, control=control$control.form)
    llk.diss <- logLikNull(object$dissolution.fit, control=control$control.diss)

    llk <- llk.form + llk.diss
    class(llk) <- "logLik"
    attr(llk,"df") <- attr(llk.form,"df") + attr(llk.diss,"df")
    attr(llk,"nobs") <- attr(llk.form,"nobs") + attr(llk.diss,"nobs")
    
    llk
}
