#  File R/gof.tergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2019 Statnet Commons
#######################################################################


#' Goodness-of-fit methods for TERGM CMLE and CMPLE fits
#' 
#' For now, these are simple wrappers around \code{\link{gof.ergm}},
#' \code{\link{print.gof}}, \code{\link{summary.gof}}, and
#' \code{\link{plot.gof}}, respectively.
#' 
#' 
#' @param object For \code{gof.tergm}, \code{\link{tergm}}
#'   conditional MLE (CMLE) or conditional MPLE (CMPLE) fit. For the
#'   others, a \code{gof.tergm} object returned by \code{gof.tergm}.
#' @param main Gives the title of the goodness-of-fit plot.
#' @param \dots Additional arguments passed through to the respective
#'   functions in the \code{\link[=ergm-package]{ergm}} package.
#' @return For \code{gof.tergm}, an object of class
#'   \code{gof.tergm}.
#' 
#' For the others, nothing.
#' @seealso [tergm()], [ergm()], [simulate.tergm()], [ergm::print.gof()], [ergm::plot.gof()],
#' summary.gof, mcmc.diagnostics.ergm
#' @keywords models
#' @examples
#' 
#' data(samplk)
#' 
#' # Fit a transition from Time 1 to Time 2
#' samplk12 <- tergm(list(samplk1, samplk2) ~ Form(~edges+mutual+transitiveties+cyclicalties) +
#'                                            Diss(~edges+mutual+transitiveties+cyclicalties),
#'                                            estimate="CMLE")
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
#' 
#' @export
gof.tergm <- function (object, ...){
  if(object$estimate=="EGMME") stop("Goodness of fit for TERGM EGMME is not implemented at this time.")
  out <- gof(object$fit,...)
  class(out)<-c("gof.tergm", class(out))
  out
}

#' @rdname gof.tergm
#' @param x A \code{gof.tergm} object returned by
#'   \code{gof.tergm}.
#' @method print gof.tergm
#' @export
print.gof.tergm <- function(x, ...){
  cat("\n================================\n")
    cat("Model goodness of fit:\n")
    cat("================================\n")

  NextMethod()
}

#' @rdname gof.tergm
#' @method summary gof.tergm
#' @export
summary.gof.tergm <- function(object, ...) {
  print(object, ...) # Nothing better for now
}

#' @rdname gof.tergm
#' @method plot gof.tergm
#' @importFrom graphics plot
#' @export
plot.gof.tergm <- function(x, ..., main="Goodness-of-fit diagnostics"){
  NextMethod()
}
