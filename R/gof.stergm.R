gof.stergm <- function (object, ...){
  if(object$estimate=="EGMME") stop("Goodness of fit for STERGM EGMME is not implemented at this time.")
  out<-list(formation=gof(object$formation.fit,...),
            dissolution=gof(object$dissolution.fit,...))
  class(out)<-"gof.stergm"
  out
}

print.gof.stergm <- function(x, ...){
  cat("\n================================\n")
    cat("Formation model goodness of fit:\n")
    cat("================================\n")

  print(x$formation, ...)
  
  cat("\n==================================\n")
    cat("Dissolution model goodness of fit:\n")
    cat("==================================\n")
  
  print(x$dissolution, ...)
}

summary.gof.stergm <- function(object, ...) {
  print.gof.stergm(object, ...) # Nothing better for now
}

plot.gof.stergm <- function(x, ..., main="Goodness-of-fit diagnostics"){
  plot(x$formation, ..., main=paste("Formation:", main))
  
  plot(x$dissolution, ..., main=paste("Dissolution:", main))
}
