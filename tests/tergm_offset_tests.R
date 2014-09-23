library(tergm)

opttest({
  data(florentine)
  net <- flobusiness
  summary(net ~ edges+degree(3))
  
  # default initialization
  mod1 <- stergm(flobusiness, formation= ~edges + offset(degree(3)), 
                 offset.coef.form=0.8, 
                 dissolution= ~offset(edges),
                 offset.coef.diss=log(9), 
                 targets="formation",
                 estimate="EGMME"
  )
  
  # init.method set to zeros works
  mod2 <- stergm(flobusiness, formation= ~edges + offset(degree(3)), 
                 offset.coef.form=0.8, 
                 dissolution= ~offset(edges),
                 offset.coef.diss=log(9), 
                 targets="formation",
                 estimate="EGMME",control=control.stergm(init.method='zeros'))
  
  # this works, auto defaulting SAN coefs if they are different from stergm init.form
  mod3 <- stergm(flobusiness, formation= ~edges + offset(degree(3)), 
                 offset.coef.form=0.8, 
                 dissolution= ~offset(edges),
                 offset.coef.diss=log(9), 
                 targets='formation', estimate="EGMME",   
                 control=control.stergm(init.form=c(-3,0.8))
  )
  
  # we can explicitly specify the target and target stats
  mod4 <- stergm(flobusiness, formation= ~edges + offset(degree(3)), 
                 offset.coef.form=0.8, 
                 dissolution= ~offset(edges),
                 offset.coef.diss=log(9), 
                 targets=~edges, target.stats = 15,
                 estimate="EGMME",   
                 control=control.stergm(init.form=c(-3,0.8))
  )
  
  sapply(list(mod1, mod2, mod3, mod4), function(x) coef(x)$formation)
  
  for (mod in list(mod1, mod2, mod3, mod4)) {
    print(apply(simulate(mod, monitor=~edges+degree(3), output="stats", time.slices=200), 2, mean))
  }
}, testname='target_offset')
