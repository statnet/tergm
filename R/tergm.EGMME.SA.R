#  File R/tergm.EGMME.SA.R in package tergm, part of the Statnet suite of
#  packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################
tergm.EGMME.SA <- function(theta0, nw, model, model.mon,
                            control, proposal, eval.optpars,
                            verbose=FALSE){

  ###### Set the constants and convenience variables. ######
  offsets <- model$etamap$offsettheta # which parameters are offsets?
  p.free <- sum(!model$etamap$offsettheta) # number of free parameters
  p <- length(model$etamap$offsettheta) # total number of parameters
  p.names <- param_names(model)
  
  q <- length(model.mon$etamap$offsettheta) # number of target statistics
  q.names<-param_names(model.mon)

  if(control$SA.plot.progress && !dev.interactive(TRUE)){
    warning("Progress plot requested on a non-interactive graphics device. Ignoring.")
    control$SA.plot.progress <- FALSE # So that we don't print a warning every step.
  }
  
  if(verbose) message("Starting optimization with with coef_0 = ( ", paste(theta0, collapse = " "), " ).")

  ###### Define the optimization run function. ######
  
  do.optimization<-function(states, history, control){
    ind.all <- history$ind.all
    tid.all <- history$tid.all
    oh.all <- history$oh.all
    jitters.all <- history$jitters.all
    
    if(verbose) message("Running stochastic optimization... ", appendLF = FALSE)
    zs <- if(!is.null(ergm.getCluster(control))){
      requireNamespace('parallel')
      # Conveniently, the first argument of tergm.EGMME.SA.Phase2.C
      # is the state of the optimization, so giving clusterApply a
      # list of states will call it for each thread's state.
      if(verbose) {message("Calling tergm.EGMME.SA.Phase2.C:"); message_print(gc())}
      out <- parallel::clusterApply(ergm.getCluster(control), states, tergm.EGMME.SA.Phase2.C, model, model.mon, proposal, control, verbose=verbose)
      if(verbose) message_print(gc())
      out
    }else{
      list(tergm.EGMME.SA.Phase2.C(states[[1]], model, model.mon, proposal, control, verbose=verbose))
    }
    if(verbose) message("Finished. Extracting.")
    for(i in seq_along(states)){

      # Extend the observation index and thread id vectors.
      ind.all <- c(ind.all, sum(tid.all==i) + 1:(control$SA.runlength*control$SA.interval))
      tid.all <- c(tid.all, rep(i, control$SA.runlength*control$SA.interval))
      
      # Extract and store the history of jitters.
      jitters.all <- rbind(jitters.all,zs[[i]]$opt.history[,p+1:p,drop=FALSE])
      colnames(jitters.all) <- p.names
      zs[[i]]$opt.history <- zs[[i]]$opt.history[,-(p+1:p),drop=FALSE]
      
      # Extract and store histroy of trials.
      oh.all <- rbind(oh.all,zs[[i]]$opt.history)
      colnames(oh.all) <- c(p.names,q.names)
    }

    # Limit the size of the full history.
    # TODO: Store a thinned or summarized version of "forgotten" iterations.
    # TODO: Add an option to log it to external storage.
    if(max(ind.all)>control$SA.oh.memory){
      remember.after.ind <- max(ind.all) - control$SA.oh.memory

      # Cut the ind.all last!
      tid.all <-tid.all[ind.all>remember.after.ind]
      oh.all <- oh.all[ind.all>remember.after.ind,,drop=FALSE]
      jitters.all <- jitters.all[ind.all>remember.after.ind,,drop=FALSE]
      ind.all <- ind.all[ind.all>remember.after.ind]

      # Shift the indices so that they would start at 1.
      ind.all <- ind.all - remember.after.ind
    }
    
    history$ind.all <- ind.all
    history$tid.all <- tid.all
    history$jitters.all <- jitters.all
    history$oh.all <- oh.all
   
    min.ind.last <- max(ind.all) - control$SA.runlength*control$SA.interval + 1
    min.ind.keep <- max(ind.all) - max(
                                     max(ind.all)*control$SA.keep.oh,
                                     min(
                                       max(
                                         control$SA.runlength*control$SA.interval*control$SA.keep.min.runs*length(states),
                                         control$SA.keep.min*length(states)
                                         ),
                                       max(ind.all)
                                       )
                                     ) + 1
    
    # Extract and store subhistories of interest.
    
    history$ind <- ind <- ind.all[ind.all>=min.ind.keep]
    history$tid <- tid <-tid.all[ind.all>=min.ind.keep]
    history$ind.last <- ind.last <-ind.all[ind.all>=min.ind.last]
    history$tid.last <- tid.last <- tid.all[ind.all>=min.ind.last]
    history$oh <- oh <- oh.all[ind.all>=min.ind.keep,,drop=FALSE]
    history$oh.last <- oh.last <- oh.all[ind.all>=min.ind.last,,drop=FALSE]
    history$jitters <- jitters <- jitters.all[ind.all>=min.ind.keep,,drop=FALSE]
    history$jitters.last <- jitters.last <- jitters.all[ind.all>=min.ind.last,,drop=FALSE]      

    rm(ind.all, tid.all, jitters.all, oh.all); gc()

    # Plot if requested.
    if(control$SA.plot.progress && dev.interactive(TRUE)){
      requireNamespace('lattice')
      
      get.dev("progress.plot")
      
      thin <- (nrow(oh)-1)%/%(control$SA.max.plot.points/length(states)) + 1
      cols <- floor(sqrt(ncol(oh)))
      layout <- c(cols,ceiling(ncol(oh)/cols))

      #' @importFrom coda mcmc.list
      suppressWarnings(print(lattice::xyplot(window(do.call(mcmc.list,by(as.data.frame(oh),INDICES=list(tid=tid),mcmc,start=min.ind.keep)), thin=thin), panel = function(...) {lattice::panel.xyplot(...);lattice::panel.abline(0, 0)}, as.table = TRUE, layout = layout, xlab=NULL)))
    }
    
    # Extract and return the "states" and the "history".
    list(states = lapply(zs, function(z) list(nw = z$newnetwork,
           nw.diff = z$nw.diff,
           eta = z$eta)
           ),
         history = history)
  }

  # This function essentially runs the chain forward with gradient etc. set to 0.
  do.dummy.run <- function(states, history, control, burnin, steps, eta=NULL){
        control$GainM <- matrix(0, nrow=p, ncol=q)
        control$dejitter <- matrix(0, nrow=p, ncol=p) # Dejitter tries to cancel the effect of jitter on the optimizer.
        control$par.guard <- control$jitter<-rep(0,p)
        control$dev.guard <- rep(0,q)
        control$dev.guard[] <- 1e10 # A huge value
        control$par.guard[!offsets] <- 1e10

        
        control$SA.keep.min.runs <- 0

        control$SA.runlength <- 1
        control$SA.burnin <- burnin
        control$SA.interval <- steps


        
        for(i in seq_along(states)){
          if(!is.null(eta)) states[[i]]$eta <- eta
        }

        do.optimization(states,history,control)
  }
  
  interpolate.par <- function(h.fit, w=diag(1,nrow=ncol(h.fit))){
    x <- t(h.fit[-1,,drop=FALSE])
    y <- -cbind(h.fit[1,])

    c(solve(t(x)%*%w%*%x)%*%t(x)%*%w%*%y)
  }

  #' @importFrom MASS ginv

  V.sandwich <- function(w, G, V.stat=ginv(w)){
    solve(t(G)%*%w%*%G)%*%t(G)%*%w%*%V.stat%*%w%*%G%*%solve(t(G)%*%w%*%G)
  }

  best.states <- function(states, history){
    w <- ginv(cov(history$oh.all[,-(1:p),drop=FALSE]))
    best.i <- which.min(mahalanobis((history$oh.all[-1,-(1:p),drop=FALSE]+history$oh.all[-nrow(history$oh.all),-(1:p),drop=FALSE])/2,0,w,inverted=TRUE))
    best.par <- history$oh.all[best.i,1:p][!offsets]

    lapply(states, function(state){
      if(p.free) state$eta[!model$etamap$offsettheta] <- best.par[seq_len(p.free)]
      state
    })
  }


  
  ##### Construct the initial state. ######

  history <- list()
  history$ind.all <- history$tid.all <- history$oh.all <- history$jitters.all <- state <- NULL

  for(restart in 1:control$SA.restarts){    
    nw.diff <- model.mon$nw.stats - model.mon$target.stats  # nw.diff keeps track of the difference between the current network and the target statistics.
    
    states <- replicate(nthreads(control),
                        {
                          list(nw=nw,
                               eta = ergm.eta(theta0, model$etamap),
                               nw.diff  = nw.diff)
                        },
                        simplify=FALSE
                        )
    
    message('========  Phase 1: Burn in, get initial gradient values, and find a configuration under which all targets vary. ========')
    
    ###### Set up and run the burn-in. ######

    control$changes <- control$collect <- FALSE
    control.phase1<-control
    control.phase1$time.samplesize <- 1
    control.phase1$time.burnin <- control$SA.burnin
    control.phase1$time.interval <- 1
    
    message("Burning in... ", appendLF = FALSE)
    
    zs <- if(!is.null(ergm.getCluster(control))){
      requireNamespace('parallel')
      if(verbose) {message("Calling tergm_MCMC_sample:"); message_print(gc())}
      out <- parallel::clusterApply(ergm.getCluster(control), seq_along(states), function(i) tergm_MCMC_sample(states[[i]]$nw, model, model.mon, proposal, eta=states[[i]]$eta, control=control.phase1, verbose=verbose))
      if(verbose) message_print(gc())
      out
    }else{
      list(tergm_MCMC_sample(states[[1]]$nw, model, model.mon, proposal, eta=states[[1]]$eta, control=control.phase1, verbose=verbose))
    }
    
    message("Done.")
    
    # Update the state with burn-in results.
    for(i in seq_along(zs)){
      states[[i]]$nw <- zs[[i]]$newnetwork
      states[[i]]$nw.diff <- states[[i]]$nw.diff + zs[[i]]$statsmatrix.mon[NROW(zs[[i]]$statsmatrix.mon),]
    }
    
    ###### Gradient estimation and getting to a decent configuration. ######
    
    ## Here, "decent configuration" is defined as one where all target
    ## statistics have some variability --- none are stuck.
    
    # Start with jitter-only.
    
    control$GainM <- matrix(0, nrow=p, ncol=q)
    control$dejitter <- matrix(0, nrow=p, ncol=p) # Dejitter tries to cancel the effect of jitter on the optimizer.
    
    control$par.guard <- control$jitter<-rep(0,p)
    control$dev.guard <- rep(0,q)
    control$jitter[!offsets] <- control$SA.phase1.jitter
    control$dev.guard[] <- 1e10 # A huge value
    control$par.guard[!offsets] <- control$SA.phase1.jitter * 4
    
    ## Adjust the number of time steps between jumps using burn-in.
    edge.ages <- unlist(sapply(states, function(state) (((if(is(state$nw, "ergm_state")) state$nw$nw0 else state$nw) %n% "time") - edgelist_with_lasttoggle(state$nw)[,3] + 1)))
    control$SA.burnin<-control$SA.interval<- round(min(control$SA.max.interval, max(control$SA.min.interval, if(length(edge.ages)>0) control$SA.interval.mul*mean(edge.ages)))/2)
  
  if(is.nan(control$SA.burnin)|is.null(control$SA.burnin)|is.na(control$SA.burnin))
    control$SA.burnin <- control$SA.interval <- 10 # TODO: Kirk : check this
  
    if(verbose>1){
      message("New interval: ", control$SA.interval)
    }  
    
    for(try in 1:control$SA.phase1.tries){
      if(verbose) message('======== Attempt ', try, ' ========') else message('Attempt ', try,' :')
      for(run in 1:control$SA.phase1.minruns){
        tmp <- if(control$SA.restart.on.err) try(do.optimization(states, history, control), silent=!verbose) else do.optimization(states, history, control)
        if(inherits(tmp, "try-error") || all(apply(tmp$history$oh.last[,-(1:p),drop=FALSE],2,var)<sqrt(.Machine$double.eps))){
          message("Something went very wrong. Restarting with smaller gain.")
          control$SA.init.gain <- control$SA.init.gain * control$SA.gain.decay
          do.restart <- TRUE
          break
        }else{
          states <- tmp$states
          history <- tmp$history
          do.restart <- FALSE
        }
      }

      states <- best.states(states, history)

      control$gain <- control$SA.init.gain
      out <- if(control$SA.restart.on.err) try(eval.optpars(states, history, control, TRUE,restart>1,FALSE), silent=!verbose) else eval.optpars(states, history, control, TRUE,restart>1,FALSE)
      if(inherits(out, "try-error") || all(apply(history$oh.last[,-(1:p),drop=FALSE],2,var)<sqrt(.Machine$double.eps))){
        message("Something went very wrong. Restarting with smaller gain.")
        control$SA.init.gain <- control$SA.init.gain * control$SA.gain.decay
        do.restart <- TRUE
        break
      }else do.restart <- FALSE
      control <- out$control

      if(all(!out$ineffectual.pars) && all(!out$bad.fits)){
        message("All parameters have some effect and all statistics are moving. Proceeding to Phase 2.")
        break
      }
      if(try==control$SA.phase1.tries) stop("The optimizer was unable to find a reasonable configuration: one or more statistics are still stuck after multiple tries, and one or more parameters do not appear to have any robust effect.")
    }
    if(do.restart) next
    
    ###### Main optimization run. ######
    
    message('========  Phase 2: Find and refine the estimate. ========')
    
    for(subphase in 1:control$SA.phase2.levels.max){
      if(verbose) message('======== Subphase ', subphase, ' ========') else message('Subphase 2.', subphase, ' ', appendLF=FALSE)
      
      control$gain <- control$SA.init.gain*control$SA.gain.decay^(subphase-1)
      stepdown.count <- control$SA.stepdown.ct
      
      for(regain in 1:control$SA.phase2.repeats){
        tmp <- if(control$SA.restart.on.err) try(do.optimization(states, history, control), silent=!verbose) else do.optimization(states, history, control)
        if(inherits(tmp, "try-error") || all(apply(tmp$history$oh.last[,-(1:p),drop=FALSE],2,var)<sqrt(.Machine$double.eps))){
          message("Something went very wrong. Restarting with smaller gain.")
          control$SA.init.gain <- control$SA.init.gain * control$SA.gain.decay
          do.restart <- TRUE
          break
        }else{
          states <- tmp$states
          history <- tmp$history
          do.restart <- FALSE
        }
        
        if(verbose){
          message("New parameters:")
          for(state in states) message_print(state$eta)
        }

        ## Get updated gain and other values
        out <- if(control$SA.restart.on.err) try(eval.optpars(states, history, control, FALSE,TRUE,TRUE), silent=!verbose) else eval.optpars(states, history, control, FALSE,TRUE,TRUE)
        if(inherits(out, "try-error") || all(apply(history$oh.last[,-(1:p),drop=FALSE],2,var)<sqrt(.Machine$double.eps))){
          message("Something went very wrong. Restarting with smaller gain.")
          control$SA.init.gain <- control$SA.init.gain * control$SA.gain.decay
          do.restart <- TRUE
          break
        }else do.restart <- FALSE
        control <- out$control

        ## Run two tests here:
        ## 1) Are the estimating equations, on average, 0?
        ## 2) Does there appear to be a trend in their sum of squares?

        # Calculate the approximate estimating equation values:
        ys <- history$oh[,-(1:p),drop=FALSE]%*%out$w%*%out$G

        # This test is fast, so no need to thin.
        p.val.1 <- try(approx.hotelling.diff.test(ys)$p.value)
        if(is.na(p.val.1)) p.val.1 <- 0

        # Thin the data to keep from bogging down.
        x <- unique(round(seq(from=1,to=NROW(history$oh),length.out=control$SA.stepdown.maxn)))
        y <- sqrt(mahalanobis(ys,0,ginv(cov(ys)),inverted=TRUE)[x])
        i <- history$ind[x]
        t <- history$tid[x]

        for(thread in unique(t)) i[t==thread] <- rank(i[t==thread])

        #' @importFrom nlme gls corAR1
        fit.2 <- try(summary(gls(y~x,correlation=corAR1(form=~i|t)))$tTable[2,c(1,4)])

         if(!inherits(p.val.1, "try-error") && !inherits(fit.2, "try-error")){
          p.val.2 <- fit.2[2]/2 # We are interested in one-sided decline here.
          est.2 <- fit.2[1]
          if(est.2>0) p.val.2 <- 1-p.val.2 # If it's actually getting worse one-sided p-value is thus.
          
          if(verbose){
            message("Estimating equations = 0 p-value: ", p.val.1, " , trending: ", p.val.2, " .")
          }

          p.vals <- c(p.val.1,p.val.2)

          fisher.pval <- function(p.vals){
            p.vals <- unlist(p.vals)
            df  <- 2*length(p.vals)
            pchisq(-2*sum(log(p.vals)), df, lower.tail=FALSE)
          }
          if(fisher.pval(p.vals)>control$SA.stepdown.p){
            stepdown.count <- stepdown.count - 1
            if(stepdown.count<=0){
              if(verbose) message("Estimating equations do not significantly differ from 0 and neither they nor the parameters exhibit a significant trend. Reducing gain.")
              else message("\\", appendLF = FALSE)
              stepdown.count <- control$SA.stepdown.ct
              if(!verbose) message("")
              break
            }else{
              if(verbose) message("Estimating equations do not significantly differ from 0 and do not exhibit a significant trend.  ", stepdown.count, " / ", control$SA.stepdown.ct, "  to go.")
              else message("\\", appendLF = FALSE)
            }
          }else{
            stepdown.count <- control$SA.stepdown.ct
            if(verbose) message("Estimating equations significantly differ from 0 or exhibit a significant trend. Resetting counter.")
            else message("/", appendLF = FALSE)
          }
        }else{
          if(verbose) message("Problem testing estimating equations. Continuing with current gain.")
          else message("!/", appendLF = FALSE)

        }
        if(do.restart) break
      }

      if(do.restart) break      

      
      #### Decide whether to stop.

      ## Run through minimal number of phases.
      if(subphase<control$SA.phase2.levels.min) next
      
      ## Run three tests here:
      ## 1) Is the stochastic approximation estimate of the GMM sufficiently precise?
      ## 2) Is there strong evidence of nonlinearity?
      ## 3) Does there appear to be a trend in parameter values?

      # Test 1:
      
      par.se <- {
        h <- history$oh[,1:p,drop=FALSE][,!offsets,drop=FALSE]
        apply(h,2,function(x) sd(x)/sqrt(effectiveSize(x)))
      }

      sandwich.se <- sqrt(diag(V.sandwich(out$w,out$G,out$v)))

      if(verbose){
        message("Approximate standard error of the estimate:")
        message_print(sandwich.se)        
        message("Approximate standard error of window means:")
        message_print(par.se)
        message("par. var. / (std. var. + par. var.):")
        message_print(par.se^2/(sandwich.se^2+par.se^2))
      }
      
      # Test 2:

      ys <- history$oh[, -(1:p), drop=FALSE]
      xs <- history$oh[, 1:p, drop=FALSE][,!offsets,drop=FALSE]

      nlin.totest <- rep(c(FALSE,TRUE),c(p.free+1, p.free*2))
      
      nlin.fit <- lm(ys~xs+I(xs^2)+I(xs^3))
      nlin.nfs <- matrix(apply(cbind(resid(nlin.fit)),2,function(x) sum(tapply(x,list(history$tid),length))/sum(tapply(x,list(history$tid),effectiveSize))), nrow=1+p.free*3, ncol=q, byrow=TRUE)
      nlin.coef <- cbind(coef(nlin.fit))
      nlin.vcov <- vcov(nlin.fit)
      
      drop <- apply(is.na(nlin.coef),1,any)
      nlin.coef <- nlin.coef[nlin.totest & !drop,,drop=FALSE]
      nlin.nfs <- nlin.nfs[nlin.totest & !drop,,drop=FALSE]
      nlin.vcov <- t(nlin.vcov[nlin.totest & !drop, nlin.totest & !drop, drop=FALSE]*sqrt(c(nlin.nfs)))*sqrt(c(nlin.nfs))
      chi2 <- mahalanobis(c(nlin.coef),0,ginv(nlin.vcov),inverted=TRUE)

      p.val.1 <- pchisq(chi2, length(nlin.coef), lower.tail=FALSE)
      if(verbose) message("Local nonlinearity p-value: ", p.val.1)

      if(subphase == control$SA.phase2.levels.max){
        if(verbose) message("Maximum number of gain levels exceeded. Stopping.", appendLF = FALSE)
      }else{
        if(all(par.se^2/(sandwich.se^2+par.se^2) > control$SA.phase2.max.mc.se)){
          if(verbose) message("EGMME does not appear to be estimated with sufficient prescision. Continuing.")
          next
        }
        
        if(p.val.1 < control$SA.stop.p){
          if(verbose) message("There is evidence of local nonlinearity. Continuing.")
          next
        }
      }

      ###### If we've gotten this far, proceed to Part 3. ######

      ## Refine the estimate.
      if(inherits(state, "try-error")) stop("Something went wrong too many times. Try better starting values or reducing control$SA.init.gain.") 
      
      eta <- states[[1]]$eta
      
      eta.free <- switch(control$SA.refine,
                         mean = colMeans(history$oh[,1:p,drop=FALSE][,!offsets,drop=FALSE]),
                         linear = interpolate.par(out$oh.fit,out$w),
                         none = if(is.null(ergm.getCluster(control))) eta[!offsets] else stop("No interpolation does not make sense with multithreaded fitting."))
      if(p.free) eta[!model$etamap$offsettheta] <- eta.free[seq_len(p.free)]
      
      if(verbose){
        message("Refining the estimate using the ", control$SA.refine, " method. New estimate:")
        message_print(eta)
      }
      
      ## Sample to estimate standard error and test if we've actually arrived.
      G <- out$G
      w <- out$w
      V.par<-matrix(NA,p,p)
      if(control$SA.se){
        message('========  Phase 3: Simulate from the fit and estimate standard errors. ========')
        
        if(verbose) message("Evaluating target statistics at the estimate.")

        # This is to avoid the latest sample "pushing out" everything else.
        control$SA.keep.min <- round(nrow(history$oh)/length(states))+control$SA.phase3.samplesize.runs*control$SA.runlength*control$SA.interval
        tmp <- do.dummy.run(states, history, control, control$SA.burnin, control$SA.phase3.samplesize.runs*control$SA.runlength*control$SA.interval, eta=eta)
        if(verbose) message("Finished.")

        states <- tmp$states
        history <- tmp$history
        rm(tmp); gc()
        
        min.ind.se <- max(history$ind.all)-control$SA.phase3.samplesize.runs*control$SA.runlength*control$SA.interval + 1
        sm.mon <- history$oh.all[history$ind.all>=min.ind.se,-(1:p),drop=FALSE]
        sm.tid <- history$tid.all[history$ind.all>=min.ind.se]

        V.stat<-cov(sm.mon)
        ee <- sm.mon %*% ginv(V.stat) %*% G
        p.val.1 <- approx.hotelling.diff.test(ee)$p.value

        if(verbose) message("Estimating equation = 0 p-value: ", p.val.1)

        if(subphase == control$SA.phase2.levels.max){
          if(verbose) message("Maximum number of gain levels exceeded. Stopping.", appendLF=FALSE)
        }else{
          if(p.val.1 < control$SA.stop.p){
            if(verbose) message("Simulated values of estimating equations are not centered around 0. Continuing.")
            next
          }else{
            if(verbose) message("Simulated values of estimating equations are centered around 0. Stopping.")          
          }
        }
        
        V.par[!offsets,!offsets]<-V.sandwich(w,G,V.stat)
        sm.mons <- lapply(unique(sm.tid), function(t) sm.mon[sm.tid==t,,drop=FALSE])
        sm.mon <- do.call(mcmc.list, lapply(sm.mons, mcmc, start=control$SA.burnin))
      }else{
        V.par[!offsets,!offsets]<-V.sandwich(w,G,out$v)
        sm.mon <- NULL
      }

      mc.se <- rep(NA, p)
      mc.se[!offsets] <- par.se

      break
    }
    if(!do.restart) break # If we've gotten this far, no restart conditions have been triggered, so we are good.
  }
  
    
  list(newnetwork=if(control$SA.se) zs[[1]]$newnetwork else states[[1]]$nw,
       newnetworks=if(control$SA.se) lapply(zs,"[[","newnetwork") else lapply(states,"[[","nw"),
       init=theta0,
       covar=V.par,
       mc.se = mc.se,
       eta=eta,
       opt.history=history$oh.all,
       sample=sm.mon,
       network=nw)
}

tergm.EGMME.SA.Phase2.C <- function(state, model, model.mon,
                             proposal, control, verbose) {
  on.exit(ergm_Cstate_clear())

  model.comb <- c(model, model.mon)
  proposal$aux.slots <- model.comb$slots.extra.aux$proposal
  
  eta.comb <- c(deInf(state$eta), rep(0,model.mon$etamap$etalength))

  ergmstate <- ergm_state(state$nw, model=model.comb, proposal=proposal)

  maxedges <- max(NVL(control$MCMC.maxedges, Inf), network.edgecount(ergmstate))
  maxchanges <- max(control$MCMC.maxchanges, network.edgecount(ergmstate))
  
  z <- .Call("MCMCDynSArun_wrapper",
             ergmstate,
             as.integer(nparam(model.mon,canonical=TRUE)),
             # Parameter fitting.
             as.double(eta.comb),
             as.double(state$nw.diff),
             as.integer(control$SA.runlength),
             as.double(control$GainM),
             as.double(control$jitter),
             as.double(control$dejitter), # Add a little bit of noise to parameter guesses.
             as.double(control$dev.guard),
             as.double(control$par.guard),
             # MCMC settings.              
             as.integer(control$SA.burnin),
             as.integer(control$SA.interval),
             as.integer(control$MCMC.burnin.min),
             as.integer(control$MCMC.burnin.max),
             as.double(control$MCMC.burnin.pval),
             as.double(control$MCMC.burnin.add),
             as.integer(deInf(maxedges, "maxint")),
             as.integer(maxchanges),
             as.integer(max(verbose-1,0)),
             PACKAGE="tergm")

  if(z$status != 0) stop("DynSA errored with code ", z$status)

  z$state <- update(z$state)

  eta <- z$eta[seq_len(model$etamap$etalength)]
  names(eta) <- param_names(model, canonical = TRUE)

  newnetwork <- as.network(z$state)
  
  list(nw.diff=z$nw.diff,
       newnetwork=newnetwork,
       eta=eta,
       opt.history=matrix(z$opt.history,ncol=2*model.comb$etamap$etalength - model.mon$etamap$etalength,byrow=TRUE))
}
