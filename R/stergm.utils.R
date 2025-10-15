#  File R/stergm.utils.R in package tergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################
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
network.extract.with.lasttoggle <- function(nwd, at){
  # first obtain all tails, heads, and lasttoggle times <= at in the nwd
  valid_eids <- valid.eids(nwd)

  tails <- unlist(lapply(nwd$mel, "[[", "outl")[valid_eids])
  heads <- unlist(lapply(nwd$mel, "[[", "inl")[valid_eids])
  lts <- unlist(lapply(lapply(lapply(nwd$mel, "[[", "atl"), "[[", 
                  "active"), function(x) suppressWarnings(max(x[x <= at]))))[valid_eids]

  # which edges are active at time at?
  active_edges <- is.active(nwd, at = at, e = valid_eids)
  
  # which edges have non-default lasttoggle times? (handling NULL case separately)
  has_lasttoggle_time <- NVL2(lts, lts != -Inf, NULL)
  
  # which edges to keep in lasttoggle?
  which_to_keep <- active_edges & has_lasttoggle_time

  lttails <- tails[which_to_keep]
  ltheads <- heads[which_to_keep]
  ltlts <- lts[which_to_keep]

  ## now *additionally* get lasttoggle times (of at) for all edges that were toggled off at time at
  relevant_nonedges <- is.active(nwd, at = at - 1, e = valid_eids) & !is.active(nwd, at = at, e = valid_eids)
  
  lttails <- c(lttails, tails[relevant_nonedges])
  ltheads <- c(ltheads, heads[relevant_nonedges])
  ltlts <- c(ltlts, rep(at, sum(relevant_nonedges)))

  ## it is apparently allowed that tail or head may be 0 which means we should omit that edge....
  w <- lttails & ltheads

  lttails <- lttails[w]
  ltheads <- ltheads[w]
  ltlts <- ltlts[w]
  
  ## now ensure things have the correct type and tail < head if undirected
  ltlts <- as.integer(ltlts)
  
  if(!is.directed(nwd)) {
    lttmp <- lttails
    w2 <- lttails > ltheads
    lttails[w2] <- ltheads[w2]
    ltheads[w2] <- lttmp[w2]
  }

  lasttoggle <- matrix(c(lttails, ltheads, ltlts), nrow=length(lttails), ncol=3)
  ## this is the matrix we want except that the nodal indices are relative
  ## to the entire networkDynamic; we want them just for the network collapsed
  ## to time at, so we now adjust the indices to account for any inactive nodes
  
  ## what nodes are active at time at?
  ia <- is.active(nwd, at=at, v=seq_len(network.size(nwd)))
  w <- which(ia)

  ## if n is (the index of) a node in the complete networkDynamic, then cs_ia[n] is 
  ## the index of node n in the collapsed network, assuming n is active in that network
  cs_ia <- cumsum(ia)

  ## take the subset of lasttoggle where both nodes are active at time at
  lasttoggle <- lasttoggle[lasttoggle[,1] %in% w & lasttoggle[,2] %in% w,,drop=FALSE]

  ## recode the active nodes so they are relative to the collapsed network
  lasttoggle[,1:2] <- cs_ia[lasttoggle[,1:2]]
  
  ## ensure correct type
  storage.mode(lasttoggle) <- "integer"    

  nw <- network.collapse(nwd, at = at) #  Convert to a network network.
  nw %n% "time" <- at
  nw %n% "lasttoggle" <- lasttoggle
  nw
}


to.networkDynamic.lasttoggle <- function(nw){
  nwd <- as.networkDynamic(nw)
  if(!is.null(nw %n% "lasttoggle")){
    
    lt.edges <- edgelist_with_lasttoggle(nw)
    
    lt.edges <- lt.edges[lt.edges[,3]>as.integer(-.Machine$integer.max/2),,drop=FALSE] 

    # removed the +1 after lt.edges[,3]
    if(nrow(lt.edges)) nwd <- deactivate.edges(nwd, onset=-Inf, terminus=lt.edges[,3], e=apply(lt.edges[,1:2,drop=FALSE],1,function(e) get.edgeIDs(nwd, e[1], e[2])))
  }
  nwd<-delete.network.attribute(nwd, "time")
  nwd<-delete.network.attribute(nwd, "lasttoggle")
  #class(nwd) <- c("networkDynamic","network")
  #attr(nwd,"end") <- nw %n% "time"
  
  nwd
}

networkDynamic.apply.changes <- function(nwd, changes){
  
  # if there are no changes, just return the existing network
  if(nrow(changes)==0){
    return(nwd)
  }
  storage.mode(changes) <- "integer"
  
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
  
  rv <- cbind(rv, as.integer(-.Machine$integer.max/2)) # default time
  
  if(is(nw, "ergm_state")) {
    lt <- nw$nw0 %n% "lasttoggle"
  } else {
    lt <- nw %n% "lasttoggle"
  }
  if(NROW(lt) == 0) return(rv)
  
  lt <- lt[order(lt[,1], lt[,2]),,drop=FALSE]
  
  # if a non-default time exists, use it instead
  j <- 1
  for(i in seq_len(NROW(rv))) {
    while(j < NROW(lt) && (lt[j,1] < rv[i,1] || (lt[j,1] == rv[i,1] && lt[j,2] < rv[i,2]))) {
      j <- j + 1
    }
    if(rv[i,1] == lt[j,1] && rv[i,2] == lt[j,2]) {
      rv[i,3] <- lt[j,3]
    }
  }
  
  rv
}

unset.offset.call <- function(object){
  if(inherits(object,"call") && object[[1]]=="offset")
    object[[2]]
  else
    object
}

#' @rdname stergm.utils
#' @title An Internal Function for Extracting (Some) Formation and Dissolution Formulas from a Combined Formula
#' 
#' @description This function is used in \code{tergm.EGMME.initialfit} and also when targets or monitoring
#' formulas are specified by characters.  It makes a basic attempt to identify the
#' formation and dissolution formulas within a larger combined formula (which may also
#' include non-separable terms).  Instances of \code{Form} at the top level (which may occur
#' inside \code{offset}) contribute to the formation formula; instances of \code{Persist} and
#' \code{Diss} at the top level (which may also occur inside \code{offset}) contribute to the 
#' dissolution formula.  All other terms are regarded as non-separable; this includes instances 
#' of \code{Form}, \code{Persist}, and \code{Diss} that occur inside other operator terms, 
#' including inside \code{Offset}, and also includes all interactions at the top level (for which
#' the top level term is effectively the interaction operator \code{*} or \code{:}), 
#' whether or not they include \code{Form}, \code{Persist}, and/or \code{Diss}.  
#' The formation and dissolution formulas are obtained by adding 
#' the contributing terms, replacing \code{Form} and \code{Persist} with trivial operators that protect
#' the environments of their formula arguments but have no effect on statistics or coefficient names 
#' (meaning the formulas effectively become cross-sectional), and replacing \code{Diss} by a similar operator
#' that negates statistics.  These are included in the return value as the \code{form} and \code{pers}
#' elements of the list (the "dissolution" formula really being the persistence formula), which also includes 
#' the formula of non-separable terms as \code{nonsep}, and the formula of all terms after replacing 
#' \code{Form}, \code{Persist}, and \code{Diss} as described above as \code{all}.
#' 
#' If usage proves problematic, one may specify the monitoring and/or targets formulas explicitly 
#' (rather than by characters), and one may pass initial coefficient values for the EGMME to avoid
#' running \code{tergm.EGMME.initialfit}.
#'
#' @param formula a \code{formula}.
#'
#' @return A \code{list} containing \code{form}, \code{pers}, \code{nonsep}, and \code{all} formulas as described above.
.extract.fd.formulae <- function(formula) {
  x <- list_rhs.formula(formula)
  signs <- attr(x, "sign")
    
  form <- ~.
  pers <- ~.
  nonsep <- ~.
  all <- ~.

  for(i in seq_along(x)) {
    term <- x[[i]]
    sign <- signs[i]

    if(!is.call(term)) {
      nonsep <- append_rhs.formula(nonsep, structure(list(term), sign=sign))
      all <- append_rhs.formula(all, structure(list(term), sign=sign))
      next
    }

    offset <- identical(quote(offset), term[[1]])

    if(offset) {
      subterm <- term[[2]]
      if(!is.call(subterm)) {
        nonsep <- append_rhs.formula(nonsep, structure(list(term), sign=sign))
        all <- append_rhs.formula(all, structure(list(term), sign=sign))
        next
      }
    } else {
      subterm <- term
    }

    if(identical(quote(Form), subterm[[1]])) {
      if(offset) {
        term[[2]][[1]] <- quote(.P)
      } else {
        term[[1]] <- quote(.P)      
      }
      form <- append_rhs.formula(form, structure(list(term), sign=sign))
      all <- append_rhs.formula(all, structure(list(term), sign=sign))      
    } else if (identical(quote(Persist), subterm[[1]])) {
      if(offset) {
        term[[2]][[1]] <- quote(.P)
      } else {
        term[[1]] <- quote(.P)      
      }
      pers <- append_rhs.formula(pers, structure(list(term), sign=sign))
      all <- append_rhs.formula(all, structure(list(term), sign=sign))      
    } else if (identical(quote(Diss), subterm[[1]])) {
      if(offset) {
        term[[2]][[1]] <- quote(.M)
      } else {
        term[[1]] <- quote(.M)      
      }
      pers <- append_rhs.formula(pers, structure(list(term), sign=sign))
      all <- append_rhs.formula(all, structure(list(term), sign=sign))      
    } else {
      nonsep <- append_rhs.formula(nonsep, structure(list(term), sign=sign))
      all <- append_rhs.formula(all, structure(list(term), sign=sign))      
    }
  }
  
  if(length(form) == 3) {
    form <- form[c(1,3)]
  }

  if(length(pers) == 3) {
    pers <- pers[c(1,3)]
  }

  if(length(nonsep) == 3) {
    nonsep <- nonsep[c(1,3)]
  }

  if(length(all) == 3) {
    all <- all[c(1,3)]
  }
  
  environment(form) <- environment(formula)
  environment(pers) <- environment(formula)
  environment(nonsep) <- environment(formula)
  environment(all) <- environment(formula)
      
  list(form = form, 
       pers = pers, 
       nonsep = nonsep, 
       all = all)
}

# Sum with a Plus sign and label = I
InitErgmTerm..P <- function(nw, arglist, ..., env = baseenv()) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))
  
  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)
  
  if(length(a$formula) == 3) {
    a$formula <- a$formula[c(1,3)]
  }
  a$label <- I

  term <- as.call(c(list(as.name("Sum")),a))
  out <- call.ErgmTerm(term, env = env, nw = nw, ...)
  out$duration <- is.durational(m)
  out
}

# Sum with a Minus sign and label = I
InitErgmTerm..M <- function(nw, arglist, ..., env = baseenv()) {
  a <- check.ErgmTerm(nw, arglist,
                      varnames = c("formula"),
                      vartypes = c("formula"),
                      defaultvalues = list(NULL),
                      required = c(TRUE))

  m <- ergm_model(a$formula, nw, ..., offset.decorate=FALSE)
  
  a$formula[c(2,3)] <- c(-1, ult(a$formula))
  a$label <- I

  term <- as.call(c(list(as.name("Sum")),a))
  out <- call.ErgmTerm(term, env = env, nw = nw, ...)
  out$duration <- is.durational(m)
  out
}
