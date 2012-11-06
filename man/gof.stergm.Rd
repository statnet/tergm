\name{gof.stergm}
\alias{gof.stergm}
\alias{print.gof.stergm}
\alias{summary.gof.stergm}
\alias{plot.gof.stergm}

\title{
  Goodness-of-fit methods for STERGM CMLE and CMPLE fits
}
\description{
  For now, these are simple wrappers around \code{\link{gof.ergm}},
  \code{\link{print.gofobject}}, \code{\link{summary.gofobject}}, and
  \code{\link{plot.gofobject}}, respectively, to
  run goodness-of-fit for formation and dissolution models
  separately. This may change in the future.
}
\usage{
\method{gof}{stergm}(object, ...)
\method{print}{gof.stergm}(x, ...)
\method{summary}{gof.stergm}(object, ...)
\method{plot}{gof.stergm}(x, ..., main="Goodness-of-fit diagnostics")
}
\arguments{
  \item{object, x}{
    For \code{gof.stergm}, \code{\link{stergm}} conditional MLE (CMLE)
    or conditional MPLE (CMPLE) fit. For the others, a \code{gof.stergm}
    object returned by \code{gof.stergm}.
  }
  \item{main}{
    Gives the title of the goodness-of-fit plots, which will have
    "Formation:" and "Dissolution:" prepended to it.
  }
  \item{\dots}{
    Additional arguments passed through to the respective functions in the
    \code{\link[=ergm-package]{ergm}} package.
}
}
\value{
  For \code{gof.stergm}, an object of class \code{gof.stergm}, which is
  simply a list with two named elements: \code{formation} and
  \code{dissolution}, each of them a \code{gofobject} returned by
  \code{\link{gof.ergm}}.

  For the others, nothing.
}
\seealso{
  stergm, ergm, simulate.stergm, print.gofobject, plot.gofobject,
  summary.gofobject, mcmc.diagnostics.ergm
}
\examples{
\donttest{
data(samplk)

# Fit a transition from Time 1 to Time 2
samplk12 <- stergm(list(samplk1, samplk2),
                   formation=~edges+mutual+transitiveties+cyclicalties,
                   dissolution=~edges+mutual+transitiveties+cyclicalties,
                   estimate="CMLE")

samplk12.gof <- gof(samplk12)

samplk12.gof

summary(samplk12.gof)

plot(samplk12.gof)

plot(samplk12.gof, plotlogodds=TRUE)
}
}
\keyword{ models }
