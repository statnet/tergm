\name{stergm.getMCMCsample}
\alias{stergm.getMCMCsample}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Collects a sample of networks and returns the formation and dissolution statistics of each sample
}
\description{
\code{stergm.getMCMCsample} is a low-level internal function not intended to be called directly by end users. It collects a sample of networks and
 returns the formation and dissolution statistics of each sample, along with
 a toggle matrix of the changes needed from the original network to each
 in the sample
}
\usage{
stergm.getMCMCsample(nw, model.form, model.diss, model.mon, 
                    MHproposal.form, MHproposal.diss, 
                    eta.form, eta.diss, control, verbose)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nw}{
 a \code{\link{network}} object
}
  \item{model.form}{
a formation model, as returned by \code{\link{ergm.getmodel}}
}
  \item{model.diss}{
a dissolution model, as returned by \code{\link{ergm.getmodel}}
}
  \item{model.mon}{
%%     ~~Describe \code{model.mon} here~~
}
  \item{MHproposal.form}{
a list of parameters needed for MHproposals of the
  formations

}
  \item{MHproposal.diss}{
 a list of parameters needed for MHproposals of the
    dissolutions
}
  \item{eta.form}{
%%     ~~Describe \code{eta.form} here~~
}
  \item{eta.diss}{
%%     ~~Describe \code{eta.diss} here~~
}
  \item{control}{
list of control paramters, probably from \code{\link{control.stergm}}
}
  \item{verbose}{
logical; whether this and other functions should be verbose
}
}
\details{
This function is normally called inside \code{\link{simulate.stergm}} to prepare inputs for the C sampling code and return its results
}
\seealso{
\code{\link{simulate.stergm}}
}
\value{
   returns the MCMC sample as a list containing:
   \itemize{
     \item{statsmatrix.form: the matrix of sampled statistics for 'model.form' RELATIVE TO INITIAL NETWORK}
     \item{statsmatrix.diss: the matrix of sampled statistics for 'model.form' RELATIVE TO INITIAL NETWORK}
     \item{statsmatrix.mon: the matrix of sampled statistics for 'model.mon' RELATIVE TO INITIAL NETWORK}
     \item{newnetwork      : the final network from the sampling process}
     \item{changed         : a toggle matrix, where the first column is
                       the timestamp of the toggle and the 2nd and 3rd
                       columns are the head & tail of the toggle; this
                       is only returned if 'control'$changes is not NULL}
     \item{maxchanges      : the "MCMC Dyn workspace"; see 'maxchanges' in the
                       input param list}
    }
}

\keyword{ internal }% __ONLY ONE__ keyword per line
