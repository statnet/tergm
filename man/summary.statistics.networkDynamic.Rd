%  File man/summary.statistics.networkDynamic.Rd in package tergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2013 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{summary.statistics.networkDynamic}
\alias{summary.statistics.networkDynamic}
\alias{summary.formula}
\alias{summary}
\alias{summary.statistics}
% Note that some of these aliases overlap with summary.statistics.Rd in
% ergm. This is intentional.
\title{
  Calculation of networkDynamic statistics.
}
\description{
  A method for \code{\link[ergm]{summary.statistics}} to calculate the specified 
  statistics for an observed \code{\link{networkDynamic}}.
  See \code{\link{ergm-terms}} for more information on the statistics that may
  be specified.
}
\usage{
\method{summary.statistics}{networkDynamic}(object, at, \dots, basis = NULL)
}
\arguments{
  \item{object}{
    An \R \code{\link{formula}} object with a
    \code{\link{networkDynamic}} as its LHS. (See
    \code{\link[ergm]{summary.statistics}} for more details.)
}
  \item{at}{
    A vector of time points at which to calculate the statistics.
}
 
  \item{basis}{
    An optional \code{\link{networkDynamic}} object relative to which
    the statistics should be calculated.}
 \item{\dots}{
    Further arguments passed to or used by methods.
}
}
\value{
  A matrix with \code{length(at)} rows, one for each time point,
  containing the statistics measured on the network.
}
\seealso{
  ergm, networkDynamic, ergm-terms, summary.statistics.network
}
\keyword{ models }
