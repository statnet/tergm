%  File man/ergm-hints.Rd in package tergm, part of the Statnet suite
%  of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution
%
%  Copyright 2008-2020 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-hints}
\alias{ergm-hints}
\alias{hints-ergm}
\alias{ergm.hints}
\alias{hints.ergm}

\alias{discord}
\alias{discordTNT}
\alias{InitErgmConstraint.discord}
\alias{InitErgmConstraint.discordTNT}

\docType{methods}

\title{Hints for Temporal Exponential Family Random Graph Models}
\description{
  This page describes the hints and network sample
  space constraints that are included with the
  \code{\link[=tergm-package]{tergm}} package. For more information, and
  instructions for using hints and constraints, see
  \code{\link[ergm]{ergm-hints}} and \code{\link{ergm}}.
}
\section{Hints implemented in the \code{\link[=tergm-package]{tergm}} package}{
  \describe{
    \item{\code{discord}}{\emph{The Dynamic Discordance Hint:}
      Propose toggling discordant dyads with greater
      frequency (typically about 50 percent).  May be used in
      dynamic fitting and simulation.
    }
  }
}

\references{

  Krivitsky PN and Handcock MS (2014) A Separable Model for Dynamic Networks. \emph{Journal of the Royal Statistical Society, Series B}, 76(1): 29-46. \doi{10.1111/rssb.12014}

  Goodreau SM, Handcock MS, Hunter DR, Butts CT, Morris M (2008a).
  A \pkg{statnet} Tutorial.
  \emph{Journal of Statistical Software}, 24(8).
  \url{https://www.jstatsoft.org/v24/i08/}.

  Hunter, D. R. and Handcock, M. S. (2006)
  \emph{Inference in curved exponential family models for networks},
  Journal of Computational and Graphical Statistics.

  Hunter DR, Handcock MS, Butts CT, Goodreau SM, Morris M (2008b).
  \pkg{ergm}: A Package to Fit, Simulate and Diagnose
  Exponential-Family Models for Networks.
  \emph{Journal of Statistical Software}, 24(3).
  \url{https://www.jstatsoft.org/v24/i03/}.

  Krivitsky PN (2012). Exponential-Family Random Graph Models for Valued
  Networks. \emph{Electronic Journal of Statistics}, 2012, 6,
  1100-1128. \href{https://dx.doi.org/10.1214/12-EJS696}{\code{doi:10.1214/12-EJS696}}

  Morris M, Handcock MS, Hunter DR (2008).
  Specification of Exponential-Family Random Graph Models:
  Terms and Computational Aspects.
  \emph{Journal of Statistical Software}, 24(4).
  \url{https://www.jstatsoft.org/v24/i04/}.

}

\keyword{models}
