library(tergm)

opttest({
  data(florentine)
  net <- flobusiness
  set.seed(1)
  
  mod1 <- stergm(flobusiness, formation= ~edges + offset(degree(3)), 
                 offset.coef.form=0.8, 
                 dissolution= ~offset(edges),
                 offset.coef.diss=log(9), 
                 targets="formation",
                 estimate="EGMME",
                 control=control.stergm(parallel=3, parallel.type="PSOCK")
  )

}, testname='tergm_parallel')