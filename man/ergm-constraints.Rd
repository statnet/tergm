\name{ergm-constraints}
\alias{ergm-constraints}
\alias{constraints-ergm}
\alias{ergm.constraints}
\alias{constraints.ergm}
\alias{tergm-constraints}
\alias{constraints-tergm}
\alias{tergm.constraints}
\alias{constraints.tergm}

\alias{InitConstraint.atleat}
\alias{InitConstraint.atmost}

\docType{methods}

\title{Formation and Dissolution Constraints for Exponential Family Random Graph Models}
\description{
  \code{\link{ergm}} is used to fit linear exponential random graph
  models (ERGMs), in which
  the probability of a given network, \eqn{y}, on a set of nodes is 
  \eqn{h(y) \exp\{\eta(\theta) \cdot
    g(y)\}/c(\theta)}, where
  \eqn{h(y)} is the reference measure (usually \eqn{h(y)=1}),
  \eqn{g(y)} is a vector of network statistics for \eqn{y},
  \eqn{\eta(\theta)} is a natural parameter vector of the same 
  length (with \eqn{\eta(\theta)=\theta} for most terms), and \eqn{c(\theta)} is the
  normalizing constant for the distribution.

  This page describes the constraints (the networks \eqn{y} for which \eqn{h(y)>0}) 
  that are included with the
  \code{\link[=tergm-package]{tergm}} package.
}
\section{Possible network constraints}{
  The constraint terms currently implemented are 
  \describe{
    \item{\code{atleast(nw)}}{\emph{The Formation Constraint:}
      Preserve all ties in network \code{nw}. Only dyads that are not
      ties in \code{nw} may be changed.
    }
    \item{\code{atmost(nw)}}{\emph{The Dissolution Constraint:}
      Prevent all nonties in network \code{nw}. Only dyads that have
      ties in \code{nw} may be changed.
    }   
  }
}

\references{

  Krivitsky PN, Handcock MS (2010). A Separable Model for Dynamic
  Networks. \url{http://arxiv.org/abs/1011.1937}
  
  Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).
  A \pkg{statnet} Tutorial.
  \emph{Journal of Statistical Software}, 24(8).
  \url{http://www.jstatsoft.org/v24/i08/}.
  
  Hunter, D. R. and Handcock, M. S. (2006)
  \emph{Inference in curved exponential family models for networks},
  Journal of Computational and Graphical Statistics.
  
  Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
  \pkg{ergm}: A Package to Fit, Simulate and Diagnose
  Exponential-Family Models for Networks.
  \emph{Journal of Statistical Software}, 24(3).
  \url{http://www.jstatsoft.org/v24/i03/}.
  
  Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
  Networks. \emph{Electronic Journal of Statistics}, 2012, 6,
  1100-1128. \href{http://dx.doi.org/10.1214/12-EJS696}{\code{doi:10.1214/12-EJS696}}
  
  Morris M, Handcock MS, Hunter DR (2008).
  Specification of Exponential-Family Random Graph Models:
  Terms and Computational Aspects.
  \emph{Journal of Statistical Software}, 24(4).
  \url{http://www.jstatsoft.org/v24/i04/}.
  
}

\keyword{models}
