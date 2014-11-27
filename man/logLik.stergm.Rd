%  File man/logLik.stergm.Rd in package tergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2014 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{logLik.stergm}
\alias{logLik.stergm}
\alias{logLikNull.stergm}
\title{
  A \code{\link{logLik}} method for \code{\link[=ergm.object]{stergm}}.
}
\description{
Functions to return the log-likelihood associated with a
\code{\link[=stergm.object]{stergm}} CMLE fit, evaluating it if
necessary. See \code{\link[ergm]{logLik.ergm}} documentation for details
and caveats.
}
\usage{

\method{logLik}{stergm}(object, 
                      add=FALSE, 
                      force.reeval=FALSE, 
                      eval.loglik=add || force.reeval, 
                      control=control.logLik.stergm(), 
                      \dots)

\method{logLikNull}{stergm}(object, control=control.logLik.stergm(),
                    ...)

}
\arguments{
  \item{object}{
    A \code{\link[=stergm.object]{stergm}} fit, returned
                      by \code{\link{stergm}}, for \code{estimate="CMLE"}.
  }
  \item{add}{Logical:  If TRUE, instead of 
    returning the log-likelihood, return \code{object} with
    log-likelihood value set.
  }
  \item{force.reeval}{Logical:  If TRUE,
    reestimate the log-likelihood even if \code{object} already
    has an estiamte. 
  }
  \item{eval.loglik}{Logical:  If TRUE,
    evaluate the log-likelihood if not set on \code{object}.
  }
  \item{control}{
    A list of control parameters for algorithm
    tuning. Constructed using \code{\link{control.logLik.ergm}}. 
  }
  \item{\dots}{
    Other arguments to the likelihood functions.
  }
}
\details{If the log-likelihood was not computed for \code{object},
  produces an error unless \code{eval.loglik=TRUE}
}
\value{
  For \code{logLik.stergm}, \code{add=FALSE} (the default), a \code{\link{logLik}} object. If
  \code{add=TRUE} (the default), an \code{\link[=ergm.object]{ergm}}
  object or a \code{\link[=stergm.object]{stergm}} object with the
  log-likelihood set. For \code{logLikNull.stergm}, a
  \code{\link{logLik}} object.
}
\references{
  Hunter, D. R. and Handcock, M. S. (2006)
  \emph{Inference in curved exponential family models for networks},
  Journal of Computational and Graphical Statistics.
}

\seealso{
  \code{\link{logLik}}, \code{\link{ergm.bridge.llr}}, \code{\link{ergm.bridge.dindstart.llk}}
}
\keyword{ models  }
