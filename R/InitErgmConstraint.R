
InitErgmConstraint.maxedges <- function(nw, arglist, ...) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("n"),
                      vartypes = c("numeric"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  if(length(a$n) != 1 || as.double(a$n) <= 0) {
    ergm_Init_abort("Argument `n` to `maxedges` constraint must be a single positive number.")
  }

  list(n = as.double(a$n))
}
