library(tergm)

opttest({
  data(florentine)
  net <- flobusiness
  set.seed(1)
  
  mod1 <- stergm(flobusiness, formation= ~edges + degree(3), 
                 dissolution= ~offset(edges),
                 offset.coef.diss=log(9), 
                 targets="formation",
                 estimate="EGMME",
                 control=control.stergm(parallel=3, parallel.type="PSOCK")
  )

}, testname='tergm_parallel')

opttest({
  data(florentine)
  net <- flobusiness
  set.seed(1)
  
  mod1 <- stergm(flobusiness, formation= ~edges + degree(3), 
                 dissolution= ~offset(edges),
                 offset.coef.diss=log(9), 
                 targets="formation",
                 estimate="EGMME",
                 control=control.stergm(parallel=3, parallel.type="MPI")
  )
  
}, testname='tergm_parallel_MPI', testvar="ENABLE_MPI_TESTS")