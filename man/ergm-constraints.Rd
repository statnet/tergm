%  File man/ergm-constraints.Rd in package tergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2008-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-constraints}
\alias{ergm-constraints}
\alias{constraints-ergm}
\alias{ergm.constraints}
\alias{constraints.ergm}

\alias{InitConstraint.atleast}
\alias{InitConstraint.atmost}

\docType{methods}

\title{Formation and Dissolution Constraints for Exponential Family Random Graph Models}
\description{
  This page describes the network sample space constraints
  that are included with the
  \code{\link[=tergm-package]{tergm}} package. For more information, and
  instructions for using constraints, see
  \code{\link[ergm]{ergm-constraints}} and \code{\link{ergm}}.
}
\section{Constraints implemented in  the \code{\link[=tergm-package]{tergm}} package}{
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
