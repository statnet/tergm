#  File R/combine.networks.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2003-2014 Statnet Commons
#######################################################################
all.same <- function(x){
  if(length(x)==0) return(TRUE)
  v0 <- x[1]
  for(v in x[-1]) if(!identical(v0,v)) return(FALSE)
  return(TRUE)
}

# create a single block-diagonal network by combining multible networks
combine.networks <- function(nwl, ignore.nattr=c("bipartite","directed","hyper","loops","mnext","multiple","n"), ignore.vattr=c(), ignore.eattr=c(), blockname="NetworkID", detect.edgecov=TRUE, standardized=FALSE){
  if(any(sapply(nwl, is.bipartite))) stop("Bipartite networks are not supported at this time.")
  if(any(diff(sapply(nwl, is.directed)))) stop("All networks must have the same directedness.")
  

  if(!standardized) nwl <- lapply(nwl, standardize.network)
  ns <- sapply(nwl, network.size)
  blks <- c(0, cumsum(ns))

  out <- network.initialize(sum(ns), directed=is.directed(nwl[[1]]))


  # Concatenate network attributes. If you run into what looks like a covariate matrix, combine it correctly.
  
  for(a in setdiff(Reduce(intersect,lapply(nwl, list.network.attributes)),
                          ignore.nattr)){ # I.e., iterate through common attributes.
    vl <- lapply(nwl, get.network.attribute, a, unlist=FALSE)

    # Here, try to autodetect covariate matrices and combine them.
    if(detect.edgecov
       && all(sapply(vl, is.matrix))
       && all(sapply(vl, nrow)==ns)
       && all(sapply(vl, ncol)==ns)
       && all.same(sapply(vl, mode))){
      
      m <- matrix(NA, sum(ns), sum(ns))
      mode(m) <- mode(vl[[1]])
      
      for(b in seq_along(vl)){
        inds <- (blks[b]+1):blks[b+1]
        m[inds, inds] <- vl[[b]]
      }

      vl <- m
    }
    
    out <- set.network.attribute(out, a, vl)
  }

  # Concatenate vertex attributes.
  
  for(a in setdiff(Reduce(intersect,lapply(nwl, list.vertex.attributes)),
                          ignore.vattr)){ # I.e., iterate through common attributes.
    out <- set.vertex.attribute(out, a,
                                do.call(c, lapply(nwl, get.vertex.attribute, a, unlist=FALSE))
                                )
  }

  # Add ties and attributes

  for(b in seq_along(nwl)){
    el <- rbind(as.edgelist(nwl[[b]]),as.edgelist(is.na(nwl[[b]])))
    eids <- apply(el, 1, function(e) get.edgeIDs(nwl[[b]], e[1], e[2], na.omit=FALSE))

    vals <- lapply(nwl[[b]]$mel,"[[","atl")[eids]
    names <- lapply(vals, names)

    out <- add.edges(out, el[,1]+blks[b], el[,2]+blks[b], names.eval=names, vals.eval=vals)
  }

  # Finally, add a vertex attribute specifying the blocks

  out <- set.vertex.attribute(out, blockname, rep(seq_along(ns),ns))

  out
}
