%  File man/ergm-terms.Rd in package tergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2008-2017 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-terms}
\alias{ergm-terms}
\alias{ergm.terms}
\alias{terms-ergm}
\alias{terms.ergm}
\alias{edge.ages}
\alias{mean.age}
\alias{edges.ageinterval}
\alias{edgecov.ages}
\alias{edgecov.mean.age}
\alias{degree.mean.age}
\alias{degrange.mean.age}
\alias{InitErgmTerm.edge.ages}
\alias{InitErgmTerm.edgecov}
\alias{InitErgmTerm.edgecov.ages}
\alias{InitErgmTerm.edges.ageinterval}
\alias{InitErgmTerm.mean.age}
\alias{InitErgmTerm.edgecov.mean.age}
\alias{InitErgmTerm.degrange.mean.age}
\alias{InitErgmTerm.degree.mean.age}

\docType{methods}
\title{Temporally-Sensitive Terms used in Exponential Family Random Graph Models}
\description{
  Unlike ordinary \code{\link[ergm]{ergm-terms}}, which take only a single
  network as an argument, the terms documented here also take into
  account the "ages" of extant ties in the network: the time elapsed since
  their formation.

  As implemented, many of these terms cannot be used to "drive" the process
  of network evolution, but they can be used as target statistics to
  infer the terms that do. More concretely, they may appear in
  \code{targets=} or \code{monitor=} formulas of
  \code{\link{stergm}}, \code{\link{simulate.stergm}}, or
  \code{\link[ergm]{summary.formula}} (with an ERGM formula), but they may
  not appear in their \code{\link[=stergm]{formation=}} and
  \code{\link[=stergm]{dissolution=}} formulas. These terms are marked with
  "(target-only)".

  All terms listed here are binary.
}

\section{Terms to represent network statistics included in the
  \code{\link[=tergm-package]{tergm}} package}{
  \describe{

    \item{\code{degrange.mean.age(from, to=+Inf, byarg=NULL, emptyval=0)}
      (target-only)}{\emph{Average age of ties incident on nodes having
      degree in a given range:} 
      The \code{from} and \code{to} arguments are vectors of distinct
      integers or \code{+Inf}, for \code{to}. If one of the vectors has
      length 1, it is recycled to the length of the other. Otherwise, they
      must have the same length. This term adds one
      network statistic to the model for each element of \code{from} (or \code{to}); the \eqn{i}th
      such statistic equals  the average, among all ties incident on nodes
      with degree greater than or equal to 
      \code{from[i]} but strictly less than \code{to[i]}, of the amount of time elapsed
      since the tie's formation. The optional argument
      \code{by} is a character string giving the name of an attribute in the
      network's vertex attribute list. If specified, then separate degree 
      statistics are calculated for nodes having each separate
      value of the attribute.
      
      Because this average is undefined for a network that does not have
      any actors with degree in the specified range, the
      argument \code{emptyval} can be used to specify the value returned
      if this is the case. This is, technically, an arbitrary value, but it should
      not have a substantial effect unless a non-negligible fraction of
      networks at the parameter configuration of interest has no actors
      with specified degree. }

    \item{\code{degree.mean.age(d, by=NULL,
	emptyval=0)} (target-only)}{\emph{Average age of ties incident
	on nodes having a given degree:}
      The \code{d} argument is a vector of distinct integers. This term adds one
      network statistic to the model for each element in \code{d}; the \eqn{i}th
      such statistic equals the average, among all ties incident on nodes
      with degree exactly \code{d[i]}, of the amount of time elapsed
      since the tie's formation. The optional argument
      \code{by} is a character string giving the name of an attribute in the
      network's vertex attribute list. If specified, then separate degree 
      statistics are calculated for nodes having each separate
      value of the attribute.
      
      Because this average is undefined for a network that does not have
      any actors with degree \code{d[i]}, the
      argument \code{emptyval} can be used to specify the value returned
      if this is the case. This is, technically, an arbitrary value, but it should
      not have a substantial effect unless a non-negligible fraction of
      networks at the parameter configuration of interest has no actors
      with specified degree.  }

    
    \item{\code{edges.ageinterval(from, to=+Inf)} (dissolution- and target-only)}{\emph{Number of edges
      with age falling into a specified range:}
      This term counts the number of edges in the network for
      which the time elapsed since formation is greater than or equal to
      \code{from} but strictly less than \code{to}. In other words, it
      is in the semiopen interval \code{[from, to)}. \code{from} and
      \code{to} may be scalars, vectors of the same length, or one of
      them must have length one, in which case it is recycled.

      When used in the dissolution formula of a STERGM, it can be used
      to model a non-Markovian dissolution process, controlling the
      hazard function in the interval directly.
    }

    \item{\code{edge.ages} (target-only)}{\emph{Sum of ages of extant
	ties:} 
      This term adds one statistic equaling sum, over all ties
      present in the network, of the amount of time elapsed since
      formation.

      Unlike \code{\link{mean.age}}, this statistic is well-defined on
      an empty network. However, if used as a target, it appears to
      produce highly biased dissolution parameter estimates if the goal
      is to get an intended average duration.
    }

    \item{\code{edgecov.ages(x, attrname=NULL)} (target-only)}{\emph{Weighted sum of ages of extant
	ties:} 
      This term adds one statistic equaling sum, over all ties
      present in the network, of the amount of time elapsed since
      formation, multiplied by a dyadic covariate. See the help for the
      \code{\link[ergm:ergm-terms]{edgecov}} term for details for
      specifying the covariate.
      
      "Weights" can be negative.
      
      Unlike \code{\link{edgecov.mean.age}}, this statistic is well-defined on
      an empty network. However, if used as a target, it appears to
      produce highly biased dissolution parameter estimates if the goal
      is to get an intended average duration.
    }
    
    \item{\code{edgecov.mean.age(x, attrname=NULL,
	emptyval=0)} (target-only)}{\emph{Weighted average age of an extant
	tie:} 
      This term adds one statistic equaling the average, over all ties
      present in the network, of the amount of time elapsed since
      formation, weighted by a (nonnegative) dyadic covariate. See the help for the
      \code{\link[ergm:ergm-terms]{edgecov}} term for details for
      specifying the covariate.
      
      The behavior when there are negative weights is undefined.
      
      Because this average is undefined for an empty network (or a
      network all of whose extant edges have been weighted 0), the
      argument \code{emptyval} can be used to specify the value returned
      if this is the case. This is, technically, an arbitrary value, but it should
      not have a substantial effect unless a non-negligible fraction of
      networks at the parameter configuration of interest is empty
      and/or if only a few dyads have nonzero weights.  }
    
    
    \item{\code{mean.age(emptyval=0)} (target-only)}{\emph{Average age of an extant
	tie:} 
      This term adds one statistic equaling the average, over all ties
      present in the network, of the amount of time elapsed since
      formation.
      
      Because this average is undefined for an empty network, the
      argument \code{emptyval} can be used to specify the value returned
      if it is. This is, technically, an arbitrary value, but it should
      not have a substantial effect unless a non-negligible fraction of
      networks at the parameter configuration of interest is empty.  }

    
  }
}
\references{
  \itemize{
    \item Handcock M. S., Hunter D. R., Butts C. T., Goodreau S. G.,
      Krivitsky P. N. and Morris M.
(2012). _Fit, Simulate and Diagnose Exponential-Family Models for
Networks_. Version 3.1.  Project home
page at <URL: http://www.statnet.org>, <URL:
CRAN.R-project.org/package=ergm>.
  
\item  Krivitsky, P.N. (2012). Modeling of Dynamic Networks based on
  Egocentric Data with Durational Information. \emph{Pennsylvania State
  University Department of Statistics Technical Report},
  2012(2012-01). \url{http://stat.psu.edu/research/technical-report-files/2012-technical-reports/modeling-of-dynamic-networks-based-on-egocentric-data-with-durational-information}

\item  Krivitsky, P.N. (2012). Modeling Tie Duration in ERGM-Based Dynamic
  Network Models. \emph{Pennsylvania State University Department of
  Statistics Technical Report},
  2012(2012-02). \url{http://stat.psu.edu/research/technical-report-files/2012-technical-reports/TR1201A.pdf}
  }
}
\seealso{\code{\link[ergm:ergm-terms]{ergm-terms}} (from the \code{\link[=ergm-package]{ergm}} package), \code{\link{ergm}}, \code{\link{network}}, \code{\link{\%v\%}}, \code{\link{\%n\%}}}
\keyword{models}
