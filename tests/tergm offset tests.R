library(tergm)

opttest({
  net <- network.initialize(10,directed=F)
  net[1,2] <- net[2,3] <- net[3,4] <- net[4,5] <- net[5,6] <- 1 
  
  # this works, specifying initial coefs; note the difference between init.form and san
  mod1 <- stergm(net, formation= ~edges + offset(degree(3)), 
                 offset.coef.form=-4, 
                 target.stats = 10, edapprox=T,
                 dissolution= ~offset(edges),
                 offset.coef.diss=3, 
                 targets=~edges, estimate="EGMME",   
                 control=control.stergm(init.form=c(-2,-4), SAN.control=control.san(coef=-2))
  )
  
  # init.method set to zeros works
  mod1 <- stergm(net, formation= ~edges + offset(degree(3)), 
                 offset.coef.form=-4, 
                 target.stats = 10, edapprox=T,
                 dissolution= ~offset(edges),
                 offset.coef.diss=3, 
                 targets=~edges, estimate="EGMME", control=control.stergm(init.method='zeros'))
  
  # this works, auto defaulting SAN coefs if they are different from stergm init.form
  mod1 <- stergm(net, formation= ~edges + offset(degree(3)), 
                 offset.coef.form=-4, 
                 target.stats = 10, edapprox=T,
                 dissolution= ~offset(edges),
                 offset.coef.diss=3, 
                 targets=~edges, estimate="EGMME",   
                 control=control.stergm(init.form=c(-2,-4))
  )
}, testname='target_offset')
