# Keeps a list of "named" graphic devices.
#
# Usage: get.dev(name)
#
# If a graphic device with a given name exists, switch to it. If not,
# tries to grab an "unnamed" device and gives it a name. If none are
# available, creates a new device, switches to it, and "remembers" the
# name.

get.dev <- local({
  devs <- list()
  function(name){    
    if(is.null(devs[[name]]) || dev.set(devs[[name]])!=devs[[name]]){
      # Try to find an "unnamed" device to take over.
      free.devs <- setdiff(dev.list(),unlist(devs))
      if(length(free.devs))
        dev.set(free.devs[1])
      else
        dev.new()
      
      devs[[name]] <<- dev.cur()
    }else dev.set(devs[[name]])
    return(devs[[name]])
  }
})

# A customized level plot that uses cyan for positive values, red for
# negative, white for close to 0, black for exactly 0, and ensures
# that the scales on both sides are the same.
.my.levelplot <- function(m,levels=80,...){
  bound <- max(na.omit(c(abs(m))))

  levelplot(m, at=unique(c(seq(-bound,-.Machine$double.eps,length.out=levels/2+1),
                 seq(.Machine$double.eps,+bound,length.out=levels/2+1))),
            col.regions=c(hsv(h=0,s=seq(1,2/levels,length.out=levels/2),v=1),rgb(0,0,0),
              hsv(h=.5,s=seq(2/levels,1,length.out=levels/2),v=1)),
            ...)
}


# A wrapper around network.extract
network.extract.with.lasttoggle <- function(nwd, at){
  nw <- network.extract(nwd, at=at)
  nw %v% ".networkDynamicID" <- which(is.active(nwd, at=at, v=seq_len(network.size(nwd))))
  
  # There is probably a more efficient way to do this:

  # Note that nw is still a networkDynamic, and still has vertex and
  # edge activity.
  lttails <- lapply(nw$mel, "[[", "outl")
  ltheads <- lapply(nw$mel, "[[", "inl")
  # I.e., from the edge attribute list, grab the "active" matrix, and from it, extract the highest (most recent) timestamp that's prior to or at at.
  # Note that if x[x<=at] is empty, max(x[x<=at]) returns -Inf, which is what we want.
  ltlts <- lapply(lapply(lapply(nw$mel, "[[", "atl"), "[[", "active"), function(x) suppressWarnings(max(x[x<=at])))

  # Now, strip the networkDynamic information from it.
  delete.vertex.attribute(nw, "active")
  delete.edge.attribute(nw, "active")
  class(nw) <- "network"
  
  ltm <-
    if(is.bipartite(nw)) m <- matrix(-Inf, nw%n%"bipartite", network.size(nw) - nw%n%"bipartite")
    else m <- matrix(-Inf, network.size(nw), network.size(nw))

  for(i in seq_along(ltlts))
    if(ltlts[[i]]!=-Inf){
      e<-c(lttails[[i]],ltheads[[i]])
      if(!all(e)) next
      if(!is.directed(nw)) e <- c(min(e),max(e))
      if(is.bipartite(nw)) e[2] <- e[2] - nw %n% "bipartite"
      # The -1 is important: lasttoggle is shifted by -1 relative to
      # networkDynamic (at least for now).
      m[e[1],e[2]] <- ltlts[[i]] - 1 
    }
  m[m==-Inf] <- round(-.Machine$integer.max/2)
  
  nw %n% "time" <- at
  nw %n% "lasttoggle" <- to.lasttoggle.matrix(m, is.directed(nw), is.bipartite(nw))

  nw
}

to.networkDynamic.lasttoggle <- function(nw){
  nwd <- nw
  if(!is.null(nw %n% "lasttoggle")){

    lt.edges <- ergm.el.lasttoggle(nw)

    lt.edges <- lt.edges[lt.edges[,3]>round(-.Machine$integer.max/2),,drop=FALSE] 

    # The +1 after lt.edges[,3] is important: lasttoggle is shifted by -1 relative to
    # networkDynamic (at least for now).
    if(nrow(lt.edges)) nwd <- deactivate.edges(nwd, onset=-Inf, terminus=lt.edges[,3]+1, e=apply(lt.edges[,1:2,drop=FALSE],1,function(e) get.edgeIDs(nw, e[1], e[2])))
  }
  nwd<-delete.network.attribute(nwd, "time")
  nwd<-delete.network.attribute(nwd, "lasttoggle")
  class(nwd) <- c("networkDynamic","network")
  attr(nwd,"end") <- nw %n% "time"
  
  nwd
}

networkDynamic.apply.changes <- function(nwd, changes){
  ## Add edges that were never present in the initial network.
  extant.edges <- as.edgelist(nwd)
  changed.edges <- unique(changes[,c("tail","head"),drop=FALSE])
  new.edges <- changed.edges[!(paste(changed.edges[,1],changed.edges[,2]) %in% paste(extant.edges[,1],extant.edges[,2])),,drop=FALSE]
  nwd <- add.edges(nwd,as.list(new.edges[,1]),as.list(new.edges[,2]))
  
  for(i in seq_len(nrow(changes))){
    eID <- get.edgeIDs(nwd,changes[i,"tail"],changes[i,"head"])
    # The following two lines are the "correct" way to do this, but
    # since we know the direction, and since we can assume that
    # everything after attr(nwd,"end") is "censored", we can do it
    # faster.
    #
    # if(changes[i,"to"]==0) nwd <- deactivate.edges(nwd, onset=changes[i,"time"], terminus=+Inf, e=eID)
    # if(changes[i,"to"]==1) nwd <- activate.edges(nwd, onset=changes[i,"time"], terminus=+Inf, e=eID)
    
    if(changes[i,"to"]==0){
      # If we are dissolving a tie, we are changing the bottom-right
      # cell in the spell matrix. However, we can't assume that a
      # spell matrix exists for a given edge, so we need to check, and
      # add a (-Inf,time) row if it doesn't.
      am <- nwd$mel[[eID]]$atl$active
      if(is.null(am)) am <- rbind(c(-Inf,+Inf))
      am[nrow(am),2] <- changes[i,"time"]
      nwd[["mel"]][[eID]]$atl$active <- am
    }else if(changes[i,"to"]==1){
      # If we are forming a tie, we are adding a new row.  If (and
      # this shouldn't happen), there is a spell ending at the same
      # time as we are trying to start it, merge it in instead.
      am <- nwd[["mel"]][[eID]]$atl$active
      if(!is.null(am) && am[nrow(am),2]==changes[i,"time"]) am[nrow(am),2] <- +Inf
      else am <- rbind(am, c(changes[i,"time"],+Inf))
    }
    nwd[["mel"]][[eID]]$atl$active <- am
  }

  nwd
}

### A potentially faster version, that might be worth refining later.

## networkDynamic.apply.changes <- function(nwd, changes){
##   ## Add edges that were never present in the initial network.
##   extant.edges <- as.edgelist(nwd)
##   changed.edges <- unique(changes[,c("tail","head"),drop=FALSE])
##   new.edges <- changed.edges[!(paste(changed.edges[,1],changed.edges[,2]) %in% paste(extant.edges[,1],extant.edges[,2])),,drop=FALSE]
##   nwd <- add.edges(nwd,as.list(new.edges[,1]),as.list(new.edges[,2]))

##   changes <- changes[order(changes[,"tail"], changes[,"head"], changes[,"time"], changes[,"to"]),,drop=FALSE]
##   edge.changes <- aggregate(as.data.frame(changes[,c("time","to")]),by=list(tail=changes[,"tail"],head=changes[,"head"]),FUN=identity)
  
##   for(i in seq_len(nrow(edge.changes))){
##     eID <- get.edgeIDs(nwd,edge.changes[i,"tail"],edge.changes[i,"head"])
##     times <- edge.changes[i,"time"][[1]]
##     tos <- edge.changes[i,"to"][[1]]
    
##     if(!all(abs(diff(tos))==1)) stop("Problem with change matrix.")

##     am <- nwd$mel[[eID]]$atl$active # Extant spells.
    
##     if(tos[1]==0){ # First change is a dissolution.
##       # No spell matrix:
##       if(is.null(am)) am <- rbind(c(-Inf,+Inf))

##       am.new <- matrix(c(am[nrow(am),1],
##                          times,
##                          if(tos[length(tos)]==1) Inf # If ending with a formation, spell continues forever.
##                          ),ncol=2,byrow=TRUE)
      
##       nwd$mel[[eID]]$atl$active <- rbind(am[-nrow(am),],am.new)
##     }else if(tos[1]==1){
##       am.new <- matrix(c(times,
##                          if(tos[length(tos)]==1) Inf # If ending with a formation, spell continues forever.
##                          ),ncol=2,byrow=TRUE)
##       nwd$mel[[eID]]$atl$active <- rbind(am,am.new)
##     }
##   }

##   nwd
## }
