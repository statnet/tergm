% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/InitErgmTerm.duration.R
\name{degree.mean.age-ergmTerm}
\alias{degree.mean.age-ergmTerm}
\alias{InitErgmTerm.degree.mean.age}
\title{Average age of ties incident on nodes having a given degree}
\usage{
# binary: degree.mean.age(d, byarg=NULL, emptyval=0)
}
\arguments{
\item{d}{a vector of distinct integers}

\item{byarg}{specifies a vertex attribute (see Specifying Vertex attributes and Levels (\code{?nodal_attributes}) for details.). If specified, then separate degree
statistics are calculated for nodes having each separate
value of the attribute.}

\item{emptyval}{can be used to specify the value returned if the network does not have any actors
with degree in the specified range. This is, technically, an arbitrary value, but it should
not have a substantial effect unless a non-negligible fraction of
networks at the parameter configuration of interest has no actors
with specified degree.}
}
\description{
Average age of ties incident on nodes having a given degree
}
\details{
This term adds one
network statistic to the model for each element in \code{d} ; the \eqn{i} th
such statistic equals the average, among all ties incident on nodes
with degree exactly \code{d[i]} , of the amount of time elapsed
since the tie's formation. The optional argument
\code{byarg} specifies a vertex attribute (see
Specifying Vertex Attributes and Levels
for details). If specified, then separate degree
statistics are calculated for nodes having each separate
value of the attribute.
}
\seealso{
\code{\link[ergm:ergmTerm]{ergmTerm}} for index of model terms currently visible to the package.

\Sexpr[results=rd,stage=render]{ergm:::.formatTermKeywords("ergmTerm", "degree.mean.age", "subsection")}
}
\concept{durational}
