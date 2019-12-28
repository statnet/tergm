#  File R/stergm.utils.R in package tergm, part of the Statnet suite
#  of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution
#
#  Copyright 2008-2019 Statnet Commons
#######################################################################
# Keeps a list of "named" graphic devices.
#
# Usage: get.dev(name)
#
# If a graphic device with a given name exists, switch to it. If not,
# tries to grab an "unnamed" device and gives it a name. If none are
# available, creates a new device, switches to it, and "remembers" the
# name.

#' @import grDevices
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

  requireNamespace('lattice', quietly=TRUE)
  lattice::levelplot(m, at=unique(c(seq(-bound,-.Machine$double.eps,length.out=levels/2+1),
                 seq(.Machine$double.eps,+bound,length.out=levels/2+1))),
            col.regions=c(hsv(h=0,s=seq(1,2/levels,length.out=levels/2),v=1),rgb(0,0,0),
              hsv(h=.5,s=seq(2/levels,1,length.out=levels/2),v=1)),
            ...)
}


# A wrapper around network.extract.
#
# Extracts the network at the specified time point `at` and attaches
# a network attribute representing the time point `at` and
# a numeric matrix named "lasttoggle" representing the last toggle 
# time of every edge that is extant at `at` as well as the last toggle
# time for non-edges that were toggled off at time at.
#
# updates: lasttoggle is NULL when duration.dependent is FALSE
network.extract.with.lasttoggle <- function(nwd, at, duration.dependent){
  nw <- network.extract(nwd, onset=at-1, terminus=at+1)

  if(duration.dependent) {
    lttails <- unlist(lapply(nw$mel, "[[", "outl"))
    ltheads <- unlist(lapply(nw$mel, "[[", "inl"))
    ltlts <- unlist(lapply(lapply(lapply(nw$mel, "[[", "atl"), "[[", 
                    "active"), function(x) suppressWarnings(max(x[x <= at]))))

    ltlts[ltlts == -Inf] <- round(-.Machine$integer.max/2)
    ltlts <- as.integer(ltlts)

    w <- lttails & ltheads

    lttails <- lttails[w]
    ltheads <- ltheads[w]
    ltlts <- ltlts[w]
    
    if(!is.directed(nw)) {
      lttmp <- lttails
      w2 <- lttails > ltheads
      lttails[w2] <- ltheads[w2]
      ltheads[w2] <- lttmp[w2]
    }

    lasttoggle <- matrix(c(lttails, ltheads, ltlts), nrow=length(lttails), ncol=3)
  } else {  # non-duration dependent model
    lasttoggle <- NULL
  }

  nw <- network.collapse(nwd, at=at) #  Convert to a network network.
  nw %n% "time" <- at
  nw %n% "lasttoggle" <- lasttoggle
  nw
}


to.networkDynamic.lasttoggle <- function(nw){
  nwd <- nw
  if(!is.null(nw %n% "lasttoggle")){
    
    lt.edges <- edgelist_with_lasttoggle(nw)
    
    lt.edges <- lt.edges[lt.edges[,3]>round(-.Machine$integer.max/2),,drop=FALSE] 

    # removed the +1 after lt.edges[,3]
    if(nrow(lt.edges)) nwd <- deactivate.edges(nwd, onset=-Inf, terminus=lt.edges[,3], e=apply(lt.edges[,1:2,drop=FALSE],1,function(e) get.edgeIDs(nw, e[1], e[2])))
  }
  nwd<-delete.network.attribute(nwd, "time")
  nwd<-delete.network.attribute(nwd, "lasttoggle")
  class(nwd) <- c("networkDynamic","network")
  #attr(nwd,"end") <- nw %n% "time"
  
  nwd
}

networkDynamic.apply.changes <- function(nwd, changes){
  
  # if there are no changes, just return the existing network
  if(nrow(changes)==0){
    return(nwd)
  }
  
  ## Add edges that were never present in the initial network.
  extant.edges <- as.edgelist(nwd)
  changed.edges <- unique(changes[,c("tail","head"),drop=FALSE])
  new.edges <- changed.edges[!(paste(changed.edges[,1],changed.edges[,2]) %in% paste(extant.edges[,1],extant.edges[,2])),,drop=FALSE]
  nwd <- add.edges(nwd,as.list(new.edges[,1]),as.list(new.edges[,2]))

  changes <- changes[order(changes[,"tail"], changes[,"head"], changes[,"time"], changes[,"to"]),,drop=FALSE]

  # Group changes by tail and head, resulting in a "data frame" with
  # columns tail, head, time, and to, with time and to being columns
  # of lists of changes for that dyad.
  changes <- aggregate(as.data.frame(changes[,c("time","to"),drop=FALSE]),by=list(tail=changes[,"tail"],head=changes[,"head"]),FUN=list)

  # Now, for each row in this "data frame", construct a list of lists,
  # each containing elements eID, times, and tos.
  changes <- apply(changes, 1, function(r)
                   list(eID=get.edgeIDs(nwd,r[["tail"]],r[["head"]]),
                        times=r[["time"]], tos=r[["to"]]))
  
  for(e in changes){
    tos <- e$tos
    times <- e$times
    eID <- e$eID
    
    if(!all(abs(diff(tos))==1)) stop("Problem with change matrix.")

    am <- nwd$mel[[eID]]$atl$active # Extant spells.
    
    if(tos[1]==0){ # First change is a dissolution.
      # No spell matrix:
      if(is.null(am)) am <- rbind(c(-Inf,+Inf))

      # If the last formation toggle is at the same time as the new
      # dissolution toggle, the whole spell gets dissolved away, with
      # the last row of am getting dropped below and not replaced by
      # anything. (This should not, normally, happen for the
      # simulate() functions.)
      #
      # Otherwise, prepend the onset of the extant tie.      
      if(am[nrow(am),1]==times[1]) times <- times[-1]
      else times <- c(am[nrow(am),1],times)

      # If ending with a formation, spell continues forever.
      if(tos[length(tos)]==1) times <- c(times, +Inf)

      # Construct a new spell matrix. (If times is empty, it's NULL, which still works.)
      am.new <- if(length(times)) matrix(times,ncol=2,byrow=TRUE)
      
      nwd[["mel"]][[eID]]$atl$active <- rbind(am[-nrow(am),],am.new)
    }else if(tos[1]==1){ # First change is a formation.

      # If the last dissolution toggle is at the same time as the new
      # formation toggle, the spell resumes as if never
      # dissolved. (This should not, normally, happen for the
      # simulate() functions.)
      if(!is.null(am) && am[nrow(am),2]==times[1]){
        times[1] <- am[nrow(am),1]
        am <- am[-nrow(am),,drop=FALSE]
      }

      # If ending with a formation, spell continues forever.
      if(tos[length(tos)]==1) times <- c(times, +Inf)
    
      # Construct a new spell matrix. (If times is empty, it's NULL, which still works.)
      am.new <- if(length(times)) matrix(times,ncol=2,byrow=TRUE)

      nwd[["mel"]][[eID]]$atl$active <- rbind(am,am.new)
    }
  }

  nwd
}

# extract edgelist with lasttoggle times for edges
edgelist_with_lasttoggle <- function(nw) {
  rv <- as.edgelist(nw)
  
  # handle no edges this way to avoid warning from cbind
  if(NROW(rv) == 0) return(matrix(0L,0,3))
  
  rv <- cbind(rv, round(-.Machine$integer.max/2)) # default time
  
  lt <- if(is(nw, "ergm_state")) nw$nw0 %n% "lasttoggle" else nw %n% "lasttoggle"
  
  # if a non-default time exists, use it instead
  for(i in seq_len(NROW(rv))) {
    for(j in seq_len(NROW(lt))) {
      if(rv[i,1] == lt[j,1] && rv[i,2] == lt[j,2]) {
        rv[i,3] <- lt[j,3]
        break
      }
    }
  }
  
  rv
}

