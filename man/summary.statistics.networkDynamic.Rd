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
  statistics for an observed \code{\link{networkDynamic}} at the specified time point(s).
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
  A matrix with \code{length(at)} rows, one for each time point in \code{at},
  and columns for each term of the formula, containing the corresponding statistics measured on the network.
}
\seealso{
  ergm, networkDynamic, ergm-terms, summary.statistics.network
}
\examples{
# create a toy dynamic network
my.nD <- network.initialize(100,directed=FALSE)
activate.vertices(my.nD, onset=0, terminus = 10)
add.edges.active(my.nD,tail=1:2,head=2:3,onset=5,terminus=8)

# use a summary formula to display number of isolates and edges
# at discrete time points
summary(my.nD~isolates+edges, at=1:10)
}
\keyword{ models }
