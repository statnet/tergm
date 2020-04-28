
#' @describeIn tergm Print the summary of the model fit.
#' @param object A \code{tergm} object.
#' @export
summary.tergm <- function(object, ...) {
  out<-list(fit=summary(object$fit,...))
  class(out)<-"summary.tergm"
  out
}


#' @noRd
#' @importFrom utils getS3method
#' @export
print.summary.tergm <- function(x, ...){
  cat("\n==============================\n")
    cat("Summary of model fit \n")
    cat("==============================\n\n")

  cat("Formula:   ")
  f <- x$fit$formula
  if(length(f)==3) f<-f[c(1,3)]
  print(f, showEnv=FALSE)
  cat("\n")

  getS3method("print","summary.ergm")(x$fit, ..., print.header=FALSE, print.formula=FALSE, print.degeneracy=FALSE, print.drop=FALSE, print.deviances=x$fit$estimate!="EGMME")

}
