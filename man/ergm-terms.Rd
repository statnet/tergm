%  File man/ergm-terms.Rd in package tergm, part of the Statnet suite
%  of packages for network analysis, https://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  https://statnet.org/attribution
%
%  Copyright 2008-2020 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{ergm-terms}
\alias{ergm-terms}
\alias{ergm.terms}
\alias{terms-ergm}
\alias{terms.ergm}
\alias{Form}
\alias{Diss}
\alias{edge.ages}
\alias{mean.age}
\alias{nodefactor.mean.age}
\alias{nodemix.mean.age}
\alias{edges.ageinterval}
\alias{edgecov.ages}
\alias{edgecov.mean.age}
\alias{degree.mean.age}
\alias{degrange.mean.age}
\alias{InitErgmTerm.Form}
\alias{InitErgmTerm.Diss}
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
  network as an argument, the terms documented here may also take into account
  the elapsed time since each dyad in the network was last toggled.  As currently 
  implemented, such "durational" terms may not be used in the formation
  or dissolution models (whether specified separably with \code{formation} and 
  \code{dissolution} arguments in the old style, or inside \code{Form} and 
  \code{Diss} operator terms in the new style), but may be used as targets, monitors,
  or summary statistics, and may also appear as non-separable effects 
  in the generative model (necessarily specified in the new style).  The ability to 
  use these terms in the formation and dissolution models may be added in the future.

  All terms listed here are binary.
}

\section{Operator Terms included in the
  \code{\link[=tergm-package]{tergm}} package}{
  \describe{
    \item{\code{Form(formula)}}{
      \emph{The Formation Operator Term:} This term accepts a model formula
      \code{formula} and produces the corresponding formation model.
    }
    
    \item{\code{Diss(formula)}}{
      \emph{The Dissolution Operator Term:} This term accepts a model formula
      \code{formula} and produces the corresponding dissolution model.
    }

    \item{\code{Cross(formula)}}{
      \emph{The Crossection Operator Term:} This term accepts a model formula
      \code{formula} and produces the corresponding model for the cross-sectional
      network. It is mainly useful for CMLE estimation, and has no effect (i.e.,
      \code{Cross(~TERM) == ~TERM}) for EGMME and dynamic simulation.
    }
  }
}

\section{Terms to represent network statistics included in the
  \code{\link[=tergm-package]{tergm}} package}{
  \describe{

    \item{\code{degrange.mean.age(from, to=+Inf, byarg=NULL, emptyval=0)}
      }{\emph{Average age of ties incident on nodes having
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
      \code{byarg} specifies a vertex attribute (see 
      \link[ergm:nodal_attributes]{Specifying Vertex Attributes and Levels}
      for details). If specified, then separate degree 
      statistics are calculated for nodes having each separate
      value of the attribute.
      
      Because this average is undefined for a network that does not have
      any actors with degree in the specified range, the
      argument \code{emptyval} can be used to specify the value returned
      if this is the case. This is, technically, an arbitrary value, but it should
      not have a substantial effect unless a non-negligible fraction of
      networks at the parameter configuration of interest has no actors
      with specified degree. }

    \item{\code{degree.mean.age(d, byarg=NULL,
	emptyval=0)}}{\emph{Average age of ties incident
	on nodes having a given degree:}
      The \code{d} argument is a vector of distinct integers. This term adds one
      network statistic to the model for each element in \code{d}; the \eqn{i}th
      such statistic equals the average, among all ties incident on nodes
      with degree exactly \code{d[i]}, of the amount of time elapsed
      since the tie's formation. The optional argument
      \code{byarg} specifies a vertex attribute (see 
      \link[ergm:nodal_attributes]{Specifying Vertex Attributes and Levels}
      for details).  If specified, then separate degree 
      statistics are calculated for nodes having each separate
      value of the attribute.
      
      Because this average is undefined for a network that does not have
      any actors with degree \code{d[i]}, the
      argument \code{emptyval} can be used to specify the value returned
      if this is the case. This is, technically, an arbitrary value, but it should
      not have a substantial effect unless a non-negligible fraction of
      networks at the parameter configuration of interest has no actors
      with specified degree.  }

    
    \item{\code{edges.ageinterval(from, to=+Inf)}}{\emph{Number of edges
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

    \item{\code{edge.ages}}{\emph{Sum of ages of extant
	ties:} 
      This term adds one statistic equaling sum, over all ties
      present in the network, of the amount of time elapsed since
      formation.

      Unlike \code{\link{mean.age}}, this statistic is well-defined on
      an empty network. However, if used as a target, it appears to
      produce highly biased dissolution parameter estimates if the goal
      is to get an intended average duration.
    }

    \item{\code{edgecov.ages(x, attrname=NULL)}}{\emph{Weighted sum of ages of extant
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
	emptyval=0)}}{\emph{Weighted average age of an extant
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
    
    
    \item{\code{mean.age(emptyval=0, log=FALSE)}}{\emph{Average age of an extant tie:} 
      This term adds one statistic equaling the average, over all ties
      present in the network, of the amount of time elapsed since
      formation.
      
      Because this average is undefined for an empty network, the
      argument \code{emptyval} can be used to specify the value returned
      if it is. This is, technically, an arbitrary value, but it should
      not have a substantial effect unless a non-negligible fraction of
      networks at the parameter configuration of interest is empty.
      
      To get mean log age instead of mean age, set \code{log=TRUE}.          
    }

    \item{\code{nodefactor.mean.age(attr, levels=NULL, emptyval=0, log=FALSE)}}{\emph{Average ages of extant half-ties incident on nodes of specified attribute levels:} 
      This term adds one statistic for each level of \code{attr}, 
      equaling the average, over all half-ties incident on nodes of that level,
      of the amount of time elapsed since formation.

      To control what levels are included, use the \code{levels} argument.  Note that the default
      \code{levels} value for \code{nodefactor.mean.age} retains all levels, unlike the default
      for \code{nodefactor}, which omits the first level.
      
      See \link[ergm:nodal_attributes]{Specifying Vertex Attributes and Levels} 
      for details on specifying vertex attributes and levels.
      
      The argument \code{emptyval} functions like it does for \code{mean.age},
      except that a different value may be specified for each level of \code{attr}.
      The length of \code{emptyval} should either be 1 (in which case that value is
      used for every level of \code{attr}) or should be equal to the number of retained
      levels of \code{attr}, in which case the \code{i}th value in \code{emptyval} is used
      for the \code{i}th retained level of \code{attr}.      
      
      To get mean log age instead of mean age, set \code{log=TRUE}; this is applied to all levels.
    }
      
    \item{\code{nodemix.mean.age(attr, b1levels=NULL, b2levels=NULL, levels=NULL, levels2=NULL, emptyval=0, log=FALSE)}}{\emph{Average ages of extant ties of specified mixing types:} 
      This term adds one statistic for each mixing type of \code{attr}, 
      equaling the average, over all ties of that mixing type,
      of the amount of time elapsed since formation.
      
      The \code{levels}-related arguments function just like they do for the ordinary \code{nodemix} term.
      See \link[ergm:nodal_attributes]{Specifying Vertex Attributes and Levels} 
      for details on specifying vertex attributes and levels.      
      
      The argument \code{emptyval} functions like it does for \code{mean.age},
      except that a different value may be specified for each mixing type of \code{attr}.
      The length of \code{emptyval} should either be 1 (in which case that value is
      used for every mixing type of \code{attr}) or should be equal to the number of retained
      mixing types of \code{attr}, in which case the \code{i}th value in \code{emptyval} is used
      for the \code{i}th retained mixing type of \code{attr}.
            
      To get mean log age instead of mean age, set \code{log=TRUE}; this is applied to all mixing types.
    }
  }
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
