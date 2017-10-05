#  File R/ergm.godfather.R in package tergm, part of the Statnet suite
#  of packages for network analysis, http://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  http://statnet.org/attribution
#
#  Copyright 2008-2017 Statnet Commons
#######################################################################
#=========================================================================
# This file contains the following 2 functions for computing changestat
# summaries of dynamic networks ??
#   <ergm.godfather>
#   <control.godfather>
#=========================================================================




###########################################################################
# <ergm.godfather>:  make the network a proposal it can't refuse. 
# Each toggle has a timestamp, and this function forces the network to make
# all of the changes at each unique timestamp value (in increasing order)
# keeping track of the change statistics that result. Thus, the final
# product is a matrix of change statistics in which the number of rows is
# determined by the # of unique timestamps and the number of columns is
# determined by the ERGM terms as usual.
#
# --PARAMETERS--
#   formula   : an ergm formula (i.e., nw ~ terms)
#   timestamps: a vector of timestamps for the given 'toggles'
#   toggles   : an edgelist of toggles that corresponds in length to
#               'timestamps'
#   sim       : a stergm sample, as returned by <stergm.getMCMCsample>
#               if 'sim' is not provided, both 'toggles' and
#               'timestamps' should be
#   start     : the start time; this is ignored if 'sim' is provided;
#               default=min(timestamps)
#   end       : the end time; this is ignored if 'sim' is provided;
#               default=max(timestamps)
#   accumulate: whether to proceed to the next timestamp without making
#               the proposed toggles (T or F); FALSE will force all toggles
#               to be realized on the network given in 'formula'
#               ?? if this is TRUE
#   verbose   : whether this and the C function should be verbose (T or F)
#               default=FALSE
#   control   : a list of additional tuning parameters for this function,
#               as returned by <control.godfather>;
#               default=<control.godfather>()
#
# --RETURNED--
#   the dynamic changestats summary as a list of the following:
#    stats     : a matrix, where the i,j entry represents the change in the
#                jth summary statistic between the original network and the
#                network at the ith unique timestamp
#    timestamps: the vector  c(NA, start:end)), where start and end are
#                specified either by attributes of 'sim' or by the 'start'
#                and 'end' inputs or default to the minimum and maximum
#                timestamps
#    newnetwork: the network after all toggles have been made if requested
#                by 'return_new_network' in <control.godfather>;
#                NULL otherwise
#
############################################################################



#' A function to apply a given series of changes to a network.
#' 
#' Gives the network a series of timed proposals it can't refuse. Returns the
#' statistics of the network, and, optionally, the final network.
#' 
#' 
#' @param formula An \code{\link{summary.formula}}-style formula, with
#'   either a \code{\link{network}} or a \code{\link{networkDynamic}}
#'   as the LHS and statistics to be computed on the RHS. If LHS is a
#'   \code{\link{networkDynamic}}, it will be used to derive the
#'   changes to the network whose statistics are wanted. Otherwise,
#'   either \code{changes} or \code{toggles} must be specified, and
#'   the LHS \code{\link{network}} will be used as the starting
#'   network.
#' @param changes A matrix with four columns: time, tail, head, and
#'   new value, describing the changes to be made. Can only be used if
#'   LHS of \code{formula} is not a \code{\link{networkDynamic}}.
#' @param toggles A matrix with three columns: time, tail, and head,
#'   giving the dyads which had changed. Can only be used if LHS of
#'   \code{formula} is not a \code{\link{networkDynamic}}.
#' @param start Time from which to start applying changes. Note that
#'   the first set of changes will take effect at
#'   \code{start+1}. Defaults to the time point 1 before the earliest
#'   change passed.
#' @param end Time from which to finish applying changes. Defaults to
#'   the last time point at which a change occurs.
#' @param end.network Whether to return a network that
#'   results. Defaults to \code{FALSE}.
#' @param stats.start Whether to return the network statistics at
#'   \code{start} (before any changes are applied) as the first row of
#'   the statistics matrix.  Defaults to \code{FALSE}, to produce
#'   output similar to that of
#'   \code{\link[=simulate.stergm]{simulate}} for STERGMs when
#'   \code{output="stats"}, where initial network's statistics are not
#'   returned.
#' @param verbose Whether to print progress messages.
#' @param control A control list generated by
#'   \code{\link{control.tergm.godfather}}.
#' @return If \code{end.network==FALSE} (the default), an
#'   \code{\link{mcmc}} object with the requested network statistics
#'   associed with the network series produced by applying the
#'   specified changes. Its \code{\link{mcmc}} attributes encode the
#'   timing information: so \code{\link{start}(out)} gives the time
#'   point associated with the first row returned, and
#'   \code{\link{end}(out)} out the last. The "thinning interval" is
#'   always 1.
#' 
#' If \code{end.network==TRUE}, return a \code{\link{network}} object with
#' \code{\link{lasttoggle}} "extension", representing the final network, with a
#' matrix of statistics described in the previous paragraph attached to it as
#' an \code{attr}-style attribute \code{"stats"}.
#' @seealso simulate.stergm, simulate.network, simulate.networkDynamic
#' @examples
#' 
#' 
#' g1 <- network.initialize(10, dir=FALSE)
#' g1[1,2] <- 1
#' g1[3,4] <- 1
#' g1 %n% "time" <- 0
#' g1 %n% "lasttoggle" <- -1-rgeom(network.dyadcount(g1),1/4)
#' 
#' dc <- matrix(rnorm(100),10,10); dc <- dc+t(dc)
#' 
#' # Simulate a network, tracking its statistics.
#' simnet <- simulate(g1, formation=~edges, dissolution=~edges, coef.form=-1, coef.diss=1,
#'                    time.slices=50, monitor=~degree(1)+mean.age+degree.mean.age(1)+
#'                                             mean.age(log=TRUE)+degree.mean.age(1,log=TRUE)+
#'                                             degrange(1,3)+mean.age+degrange.mean.age(1,3)+
#'                                             mean.age(log=TRUE)+degrange.mean.age(1,3,log=TRUE)+
#'                                             edge.ages+edgecov(dc)+edgecov.ages(dc),
#'                    output="networkDynamic")
#' 
#' sim.stats <- attr(simnet, "stats")
#' 
#' print(head(sim.stats))
#' sim.stats <- as.matrix(sim.stats)
#' 
#' # Replay the simulation using a networkDynamic, monitoring a potentially different set of
#' # statistics (but same in this case).
#' gf1.stats <- tergm.godfather(simnet~degree(1)+mean.age+degree.mean.age(1)+
#'                                     mean.age(log=TRUE)+degree.mean.age(1,log=TRUE)+
#'                                     degrange(1,3)+mean.age+degrange.mean.age(1,3)+
#'                                     mean.age(log=TRUE)+degrange.mean.age(1,3,log=TRUE)+
#'                                     edge.ages+edgecov(dc)+edgecov.ages(dc),
#'                              start=0, end=50)
#' 
#' print(head(gf1.stats))
#' gf1.stats <- as.matrix(gf1.stats)
#' 
#' # Replay the simulation using the initial network + list of changes.
#' 
#' gf2.stats <- tergm.godfather(g1~degree(1)+mean.age+degree.mean.age(1)+
#'                                 mean.age(log=TRUE)+degree.mean.age(1,log=TRUE)+
#'                                 degrange(1,3)+mean.age+degrange.mean.age(1,3)+
#'                                 mean.age(log=TRUE)+degrange.mean.age(1,3,log=TRUE)+
#'                                 edge.ages+edgecov(dc)+edgecov.ages(dc),
#'                              start=0, end=50, changes=attr(simnet,"changes"))
#' 
#' print(head(gf2.stats))
#' gf2.stats <- as.matrix(gf2.stats)
#' 
#' # We can also compare them to the network statistics summarized.
#' summ.stats <- summary(simnet~degree(1)+mean.age+degree.mean.age(1)+
#'                              mean.age(log=TRUE)+degree.mean.age(1,log=TRUE)+
#'                              degrange(1,3)+mean.age+degrange.mean.age(1,3)+
#'                              mean.age(log=TRUE)+degrange.mean.age(1,3,log=TRUE)+
#'                              edge.ages+edgecov(dc)+edgecov.ages(dc), at=1:50)
#' 
#' print(head(summ.stats))
#' 
#' tol <- sqrt(.Machine$double.eps)
#' # If they aren't all identical, we are in trouble.
#' stopifnot(all.equal(sim.stats,gf1.stats),
#'           all.equal(sim.stats,gf2.stats),
#'           all.equal(sim.stats,summ.stats))
#' 
#' 
#' @export tergm.godfather
tergm.godfather <- function(formula, changes=NULL, toggles=changes[,-4,drop=FALSE],
                           start=NULL, end=NULL,
                           end.network=FALSE,
                           stats.start=FALSE,
                           verbose=FALSE,
                           control=control.tergm.godfather()){
  check.control.class("tergm.godfather", "tergm.godfather")

  nw <- ergm.getnetwork(formula)
  
  formula <- nonsimp.update.formula(formula, nw~., from.new="nw")
  
  if(is.networkDynamic(nw)){
    if(!is.null(toggles)) stop("Network passed already contains change or toggle information.")

    toggles <- do.call(rbind, lapply(nw$mel, function(e) if(length(c(e$atl$active)[is.finite(c(e$atl$active))])) cbind(c(e$atl$active)[is.finite(c(e$atl$active))], e$outl,e$inl) else NULL))
    toggles[,1] <- ceiling(toggles[,1]) # Fractional times take effect at the end of the time step.
    
    net.obs.period<-nw%n%'net.obs.period'
    nwend <- if(!is.null(net.obs.period)) .get.last.obs.time(nw) else NULL
    nwstart<- if(!is.null(net.obs.period)) .get.first.obs.time(nw) else NULL
   
    start <- NVL(start,
                 nwstart,
                 suppressWarnings(min(toggles[,1]))-1
                 )
    if(start==Inf) stop("networkDynamic passed contains no change events or attributes. You must specify start explicitly.")

    end <- NVL(end,
               nwend,
               suppressWarnings(max(toggles[,1]))
               )
    if(end==-Inf) stop("networkDynamic passed contains no change events or attributes. You must specify end explicitly.")

    # The reason why it's > start is that the toggles that took effect
    # at start have already been applied to the network. Conversely,
    # it's <= end because we do "observe" the network at end, so we
    # need to apply the toggles that take effect then.
    toggles <- toggles[toggles[,1]>start & toggles[,1]<=end,,drop=FALSE]

    # This is important, since end is inclusive, but terminus is exclusive.
    if(!all(is.active(nw, onset=start, terminus=end+.Machine$double.eps*end*2, v=seq_len(network.size(nw)), rule="any")
            ==is.active(nw, onset=start, terminus=end+.Machine$double.eps*end*2, v=seq_len(network.size(nw)), rule="all")))
      stop("Network size and/or composition appears to change in the interval between start and end. This is not supported by ergm.godfather() at this time.")

    # Finally, we are ready to extract the network.
  duration.dependent <- if(is.durational(formula)){1} else {0}
  nw <- network.extract.with.lasttoggle(nw, at=start, duration.dependent)

  }else{
    if(is.null(toggles)) stop("Either pass a networkDynamic, or provide change or toggle information.")
      
    start <- NVL(start,
                 attr(toggles, "start"),
                 min(toggles[,1])-1
                 )
    end <- NVL(end,
               attr(toggles, "end"),
               max(toggles[,1])
               )

    # The reason why it's > start is that the toggles that took effect
    # at start have already been applied to the network. Conversely,
    # it's <= end because we do "observe" the network at end, so we
    # need to apply the toggles that take effect then.
    toggles <- toggles[toggles[,1]>start & toggles[,1]<=end,,drop=FALSE]
    
    if(is.null(nw %n% "lasttoggle")) nw %n% "lasttoggle" <- rep(round(-.Machine$integer.max/2), network.dyadcount(nw))
    nw %n% "time" <- start
  }

  if(!is.directed(nw)) toggles[,2:3] <- t(apply(toggles[,2:3,drop=FALSE], 1, sort))
  toggles <- toggles[order(toggles[,1],toggles[,2],toggles[,3]),,drop=FALSE]

  formula <- nonsimp.update.formula(formula, nw~., from.new="nw")
  m <- ergm.getmodel(formula, nw, expanded=TRUE, role="target")
  Clist <- ergm.Cprepare(nw, m)
  m$obs <- summary(m$formula)
  if(end.network){
    maxedges.sd <- sqrt(nrow(toggles)*0.25)*2 # I.e., if each toggle has probability 1/2 of being in a particular direction, this is the s.d. of the number of edges added.
    maxedges <- Clist$nedges + maxedges.sd*control$GF.init.maxedges.mul
  }

  if(verbose) cat("Applying changes...\n")
  repeat{
    z <- .C("godfather_wrapper",
            as.integer(Clist$tails), as.integer(Clist$heads),
            time = if(is.null(Clist$time)) as.integer(0) else as.integer(Clist$time),
            lasttoggle = as.integer(NVL(Clist$lasttoggle,0)),             
            as.integer(Clist$nedges),
            as.integer(Clist$n),
            as.integer(Clist$dir), as.integer(Clist$bipartite),
            as.integer(Clist$nterms), 
            as.character(Clist$fnamestring),
            as.character(Clist$snamestring),
            as.double(Clist$inputs),
            as.integer(nrow(toggles)), as.integer(toggles[,1]),
            as.integer(toggles[,2]), as.integer(toggles[,3]),
            as.integer(start), as.integer(end),
            s = double((1+end-start) * Clist$nstats),
            if(end.network) as.integer(maxedges) else as.integer(0),
            newnwtails = if(end.network) integer(maxedges+1) else integer(0),
            newnwheads = if(end.network) integer(maxedges+1) else integer(0),
            as.integer(verbose),
            status = integer(1), # 0 = OK, MCMCDyn_TOO_MANY_EDGES = 1
            PACKAGE="tergm")

    if(z$status==0) break;
    if(z$status==1){
      maxedges <- 5*maxedges
      if(verbose>0) message("Too many edges encountered in the simulation. Increasing capacity to ", maxedges)
    }
  }

  stats <- matrix(z$s + m$obs, ncol=Clist$nstats, byrow=TRUE)
  colnames(stats) <- m$coef.names
  if(!stats.start) stats <- stats[-1,,drop=FALSE]
  #' @importFrom coda mcmc
  stats <- mcmc(stats, start=if(stats.start) start else start+1)
  
  if(end.network){ 
    if(verbose) cat("Creating new network...\n")
    newnetwork <- newnw.extract(nw,z)
    newnetwork %n% "time" <- z$time
    newnetwork %n% "lasttoggle" <- z$lasttoggle

    attr(newnetwork,"stats")<-stats
    newnetwork
  }else stats
}




####################################################################
# The <control.godfather> function allows for tuning of the
# <ergm.godfather> function
#
# --PARAMETERS--
#   maxedges          : the maximum number of edges to make space
#                       for for the new network; this is ignored
#                       if 5*Clist$nedges is greater; this is also
#                       ignored if 'return_new_network' is FALSE;
#                       default=100000
#
#
# --RETURNED--
#   a list of the above parameters
#
####################################################################
#' Control parameters for [tergm.godfather()].
#'
#' Returns a list of its arguments.
#'
#' @param GF.init.maxedges.mul How much space
#'   is allocated for the edgelist of the final network. It is used
#'   adaptively, so should not be greater than \code{10}.
#' 
#' @export control.tergm.godfather
control.tergm.godfather<-function(GF.init.maxedges.mul=5
              ){
    control<-list()
    for(arg in names(formals(sys.function())))
      control[arg]<-list(get(arg))

    control <- set.control.class("control.tergm.godfather")
    control
  }
