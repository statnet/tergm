# This tests the internal function combine.networks.

n <- 10
T <- 10

library(tergm)
yl <- replicate(T,
                {
                  y <- network.initialize(n,dir=FALSE)
                  y <- simulate(y~edges, coef=-1)
                  y %n% "x" <- matrix(runif(n*n),n,n)
                  y %v% "v" <- runif(n)
                  y %e% "e" <- runif(network.edgecount(y))
                  y
                },
                simplify=FALSE)

yc <- tergm:::combine.networks(yl)
ym <- as.matrix(yc)

for(t in seq_len(T)){
  J <- I <- (t-1)*10 + seq_len(n)
  stopifnot(all(as.matrix(yc)[I,J]==as.matrix(yl[[t]]))) # Check ties.
  stopifnot(all((yc %n% "x")[I,J]==(yl[[t]] %n% "x"))) # Check dyadic attributes.
  stopifnot(all((yc %v% "v")[I]==(yl[[t]] %v% "v"))) # Check vertex attributes.
  stopifnot(all(as.matrix(yc,attr="e")[I,J]==as.matrix(yl[[t]],attr="e"))) # Check edge attributes.

  ym[I,J] <- 0
}
stopifnot(all(ym==0))


