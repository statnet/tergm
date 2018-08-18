#  File R/gof.stergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2018 Statnet Commons
#######################################################################


#' Goodness-of-fit methods for STERGM CMLE and CMPLE fits
#' 
#' For now, these are simple wrappers around \code{\link{gof.ergm}},
#' \code{\link{print.gof}}, \code{\link{summary.gof}}, and
#' \code{\link{plot.gof}}, respectively, to run goodness-of-fit for
#' formation and dissolution models separately. This may change in the future.
#' 
#' 
#' @param object For \code{gof.stergm}, \code{\link{stergm}}
#'   conditional MLE (CMLE) or conditional MPLE (CMPLE) fit. For the
#'   others, a \code{gof.stergm} object returned by \code{gof.stergm}.
#' @param main Gives the title of the goodness-of-fit plots, which
#'   will have "Formation:" and "Dissolution:" prepended to it.
#' @param \dots Additional arguments passed through to the respective
#'   functions in the \code{\link[=ergm-package]{ergm}} package.
#' @return For \code{gof.stergm}, an object of class
#'   \code{gof.stergm}, which is simply a list with two named
#'   elements: \code{formation} and \code{dissolution}, each of them a
#'   \code{gof} returned by \code{\link{gof.ergm}}.
#' 
#' For the others, nothing.
#' @seealso [stergm()], [ergm()], [simulate.stergm()], [ergm::print.gof()], [ergm::plot.gof()],
#' summary.gof, mcmc.diagnostics.ergm
#' @keywords models
#' @examples
#' 
#' \donttest{
#' data(samplk)
#' 
#' # Fit a transition from Time 1 to Time 2
#' samplk12 <- stergm(list(samplk1, samplk2),
#'                    formation=~edges+mutual+transitiveties+cyclicalties,
#'                    dissolution=~edges+mutual+transitiveties+cyclicalties,
#'                    estimate="CMLE")
#' 
#' samplk12.gof <- gof(samplk12)
#' 
#' samplk12.gof
#' 
#' summary(samplk12.gof)
#' 
#' plot(samplk12.gof)
#' 
#' plot(samplk12.gof, plotlogodds=TRUE)
#' }
#' 
#' @export
gof.stergm <- function (object, ...){
  if(object$estimate=="EGMME") stop("Goodness of fit for STERGM EGMME is not implemented at this time.")
  out<-list(formation=gof(object$formation.fit,...),
            dissolution=gof(object$dissolution.fit,...))
  class(out)<-"gof.stergm"
  out
}

#' @rdname gof.stergm
#' @param x A \code{gof.stergm} object returned by
#'   \code{gof.stergm}.
#' @method print gof.stergm
#' @export
print.gof.stergm <- function(x, ...){
  cat("\n================================\n")
    cat("Formation model goodness of fit:\n")
    cat("================================\n")

  print(x$formation, ...)
  
  cat("\n==================================\n")
    cat("Dissolution model goodness of fit:\n")
    cat("==================================\n")
  
  print(x$dissolution, ...)
}

#' @rdname gof.stergm
#' @method summary gof.stergm
#' @export
summary.gof.stergm <- function(object, ...) {
  print.gof.stergm(object, ...) # Nothing better for now
}

#' @rdname gof.stergm
#' @method plot gof.stergm
#' @importFrom graphics plot
#' @export
plot.gof.stergm <- function(x, ..., main="Goodness-of-fit diagnostics"){
  plot(x$formation, ..., main=paste("Formation:", main))
  
  plot(x$dissolution, ..., main=paste("Dissolution:", main))
}
