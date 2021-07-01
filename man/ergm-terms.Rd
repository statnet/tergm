%  File man/ergm-terms.Rd in package tergm, part of the
%  Statnet suite of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution .
%
%  Copyright 2008-2021 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-terms}
\alias{ergm-terms}
\alias{ergm.terms}
\alias{terms-ergm}
\alias{terms.ergm}
\alias{nodemix.mean.age}
\alias{InitErgmTerm.edgecov}

\docType{methods}
\title{Temporally-Sensitive Operator and Durational Terms used in Exponential Family Random Graph Models}
\description{
  The terms described here are unique to temporal networks:  each summarizes
  some type of change or durational information.
  
  The operator terms include: \code{Form()}, \code{Persist()}, \code{Diss()}, 
  \code{Cross()} and \code{Change()}.  These are used to specify how the 
  \code{\link[ergm]{ergm-terms}} in a formula are evaluated across a network
  time-series.  Note, you cannot use an operator within another operator, so
  \code{Cross(~Form(~edges))} is not a valid specification.
  
  The durational terms are distinguished either by their name, \code{mean.age}, or
  their name extensions: \code{<name>.ages}, 
  \code{<name>.mean.age}, and \code{<name>.age.interval}.  In contrast to their named
  \code{\link[ergm]{ergm-terms}} counterparts, these durational terms take into account
  the elapsed time since each (term-relevant) dyad in the network was last toggled.  
  
  As currently implemented, the package does not support use of durational 
  terms during estimation.  But durational terms may be used as targets, monitors,
  or summary statistics.  The ability to 
  use these terms in the estimation of models is under development.

  All terms listed here currently work with binary-valued ties only.
}

\section{Terms to represent network statistics included in the
  \code{\link[=tergm-package]{tergm}} package}{
}
\references{
  \itemize{
    \item Handcock M. S., Hunter D. R., Butts C. T., Goodreau S. G.,
      Krivitsky P. N. and Morris M.
(2012). _Fit, Simulate and Diagnose Exponential-Family Models for
Networks_. Version 3.1.  Project home
page at <URL: https://statnet.org>, <URL:
CRAN.R-project.org/package=ergm>.
  
\item  Krivitsky, P.N. (2012). Modeling of Dynamic Networks based on
  Egocentric Data with Durational Information. \emph{Pennsylvania State
  University Department of Statistics Technical Report},
  2012(2012-01). \url{https://web.archive.org/web/20170830053722/https://stat.psu.edu/research/technical-report-files/2012-technical-reports/TR1201A.pdf}

\item  Krivitsky, P.N. (2012). Modeling Tie Duration in ERGM-Based Dynamic
  Network Models. \emph{Pennsylvania State University Department of
  Statistics Technical Report},
  2012(2012-02).
  }
}
\seealso{\code{\link[ergm:ergm-terms]{ergm-terms}} (from the \code{\link[=ergm-package]{ergm}} package), \code{\link{ergm}}, \code{\link{network}}, \code{\link{\%v\%}}, \code{\link{\%n\%}}}
\keyword{models}
