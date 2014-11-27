%  File man/impute.network.list.Rd in package tergm, part of the Statnet suite
%  of packages for network analysis, http://statnet.org .
%
%  This software is distributed under the GPL-3 license.  It is free,
%  open source, and has the attribution requirements (GPL Section 7) at
%  http://statnet.org/attribution
%
%  Copyright 2003-2014 Statnet Commons
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
\name{impute.network.list}
\alias{impute.network.list}
\title{
  Impute missing dyads in a series of networks
}
\description{
  This function takes a list of networks with missing dyads and returns
  a list of networks with missing dyads imputed according to a list of
  imputation directives.
}
\usage{
impute.network.list(nwl, imputers = c(), nwl.prepend = list(), nwl.append = list())
}
\arguments{
  \item{nwl}{
    A list of \code{\link{network}} objects or a
    \code{\link{network.list}} object.
  }
  \item{imputers}{
    A character vector giving one or more methods to impute missing
    dyads. Currenly implemented methods are as follows:
    \describe{
      \item{\code{"next"}}{Impute the state of the same dyad in the next
	network in the list (or later, if that one is also missing).
	This imputation method is likely to
	lead to an underestimation of the formation and dissolution
	rates. The last network in the list cannot be imputed this way.}
      \item{\code{"previous"}}{Impute the state of the same dyad in the
	previous network in the list (or earlier, if that one is also
	missing). The first network in the list cannot be imputed
	this way.}
      \item{\code{"majority"}}{Impute the missing dyad with the value of
	the majority among the non-missing dyads in that time step's
	network. A network that has exactly the same number of ties as
	non-missing non-ties cannot be imputed this way.}
      \item{\code{"0"}}{Assume missing dyads are all non-ties.}
      \item{\code{"1"}}{Assume missing dyads are all ties.}
    }
    If \code{length(imputers)>1} the specified
    imputation methods will be applied in succession. For example,
    \code{imputers=c("next","previous","majority","0")} would
    first try to impute a missing dyad with the next time step's
    value. If it, and all of the later values for that dyad are missing,
    it will try to impute it with the previous time step's value. If it,
    and all of the earlier values for that dyad are missing as well, it
    will try to impute it with the value of the majority of non-missing
    dyads for that time step. If there is an exact tie, it will impute
    0.
  }
  \item{nwl.prepend}{
    An optional list of networks to treat as preceding those in
    \code{nwl}. They will not be imputed or returned, but they can be
    useful for imputing dyads in the first network in \code{nwl}, when
    using \code{"previous"} imputer.
  }
  \item{nwl.append}{
    An optional list of networks to treat as following those in
    \code{nwl}. They will not be imputed or returned, but they can be
    useful for imputing dyads in the last network in \code{nwl}, when
    using \code{"next"} imputer.
  }
}
\value{
  A list of networks with missing dyads imputed.
}
\seealso{
  \code{\link{network}}, \code{\link{is.na}}
}
\keyword{manip}

