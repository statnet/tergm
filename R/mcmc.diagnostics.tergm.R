#  File R/mcmc.diagnostics.tergm.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2020 Statnet Commons
#######################################################################





#' Conduct MCMC diagnostics on an ergm or stergm fit
#' 
#' This function prints diagnistic information and creates simple diagnostic
#' plots for the MCMC sampled statistics produced from a \code{\link{stergm}}
#' fit.
#' 
#' The plots produced are a trace of the sampled output and a density estimate
#' for each variable in the chain.  The diagnostics printed include
#' correlations and convergence diagnostics.
#' 
#' In fact, an \code{object} contains the matrix of statistics from the MCMC
#' run as component \code{$sample}.  This matrix is actually an object of class
#' \code{mcmc} and can be used directly in the \code{coda} package to assess
#' MCMC convergence. \emph{Hence all MCMC diagnostic methods available in
#' \code{coda} are available directly.} See the examples and
#' \url{http://www.mrc-bsu.cam.ac.uk/software/bugs/the-bugs-project-winbugs/coda-readme/}.
#' 
#' More information can be found by looking at the documentation of
#' \code{\link{stergm}}.
#' 
#' @param object A stergm object.  See documentation for
#'   \code{\link{stergm}}.
#' @param center Logical: If TRUE, ; center the samples on the
#'   observed statistics.
#' @param esteq Logical: If TRUE, summarize the estimating equation values
#'   (evaluated at the MLE of any non-linear parameters), rather than
#'   their canonical components.
#' @param vars.per.page Number of rows (one variable per row) per
#'   plotting page.  Ignored if \code{latticeExtra} package is not
#'   installed.
#' @param \dots Additional arguments, to be passed to plotting
#'   functions.
#' @return \code{\link{mcmc.diagnostics.ergm}} returns some degeneracy
#'   information, if it is included in the original object.  The
#'   function is mainly used for its side effect, which is to produce
#'   plots and summary output based on those plots.
#' @seealso \code{\link{ergm}}, \code{\link{stergm}},\code{network}
#'   package, \code{coda} package, \code{\link{summary.ergm}}
#' @references Raftery, A.E. and Lewis, S.M. (1995).  The number of
#'   iterations, convergence diagnostics and generic Metropolis
#'   algorithms.  In Practical Markov Chain Monte Carlo (W.R. Gilks,
#'   D.J. Spiegelhalter and S. Richardson, eds.).  London, U.K.:
#'   Chapman and Hall.
#' 
#' This function is based on the \code{coda} package It is based on the the R
#' function \code{raftery.diag} in \code{coda}.  \code{raftery.diag}, in turn,
#' is based on the FORTRAN program \code{gibbsit} written by Steven Lewis which
#' is available from the Statlib archive.
#' @keywords models
#' @importFrom ergm mcmc.diagnostics

#' @export
mcmc.diagnostics.tergm <- function(object, 
                                    center=TRUE,
                                    esteq=TRUE,
                                    vars.per.page=3, ...){
  if(!is.null(object$fit$sample)){
    cat("\n==========================\n")
    cat("Model fit diagnostics\n")
    cat("==========================\n\n")
    object$fit$control$MCMLE.termination <- NVL(object$fit$control$MCMLE.termination,"NA")
    mcmc.diagnostics(object$fit, center=center, esteq=esteq, vars.per.page=vars.per.page, ...)
  }
  if(!is.null(object$sample)){
    cat("\n==========================\n")
    cat("EGMME diagnostics\n")
    cat("==========================\n\n")
    object$control$MCMLE.termination <- NVL(object$control$MCMLE.termination,"NA")    
    getS3method("mcmc.diagnostics","ergm")(object, center=center, esteq=esteq, vars.per.page=vars.per.page, ...)
  }
}

