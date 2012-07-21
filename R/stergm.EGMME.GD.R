stergm.EGMME.GD <- function(theta.form0, theta.diss0, nw, model.form, model.diss, model.mon,
                            control, MHproposal.form, MHproposal.diss,
                            verbose=FALSE){

  eval.optpars <- function(test.G,window,update.jitter){
    for(name in ls(pos=parent.frame())) assign(name, get(name, parent.frame()))
    
    ## Regress statistics on parameters.
    # This uses GLS to account for serial correlation in statistics,
    # since we want p-values. First row is the intercept.

    h <- get(if(window) "oh" else "oh.all")
    control <- get("control")
    state <- get("state")
    x<-h[,1:p,drop=FALSE][,!offsets,drop=FALSE] # #$%^$ gls() doesn't respect I()...
    ys <- h[,-(1:p),drop=FALSE]
    
    x2 <- sweep(x,2,apply(x,2,median),"-")^2
    x3 <- sweep(x,2,apply(x,2,median),"-")^3

    h.fits <-
      if(test.G)
        sapply(1:q,
               function(i){
                 y<-ys[,i]
                 try(gls(y~x+x2+x3, correlation=corARMA(p=2)),silent=TRUE)
               },simplify=FALSE)
      else
        sapply(1:q,
             function(i){
               y<-ys[,i]
               suppressWarnings(try(lmrob(y~x+x2+x3), silent=TRUE))
             },simplify=FALSE)
    
    bad.fits <- sapply(h.fits, inherits, "try-error")
#    bad.fits <-     # Also, ignore fits where the statistics are too concentrated.    
#      bad.fits | (apply(ys,2,function(y){
#        freqs <- table(y)
#        sum(freqs[-which.max(freqs)])
#      })<nrow(h)/2)

    if(all(bad.fits)) stop("The optimization appears to be stuck. Try better starting parameters, lower SA.init.gain, etc.")
    
    ## Grab the coefficients, t-values, and residuals.
    
    h.fit <- h.pvals <- matrix(NA, nrow=p.free+1,ncol=q)
    
    h.pvals[,!bad.fits] <- if(test.G) sapply(h.fits[!bad.fits],function(fit) summary(fit)$tTable[seq_len(p.free+1),4]) else 0
    h.fit[,!bad.fits] <- sapply(h.fits[!bad.fits], coef)[seq_len(p.free+1),]

    h.resid <- matrix(NA, nrow=NROW(ys), ncol=q)
    h.resid[,!bad.fits] <- sapply(h.fits[!bad.fits], resid)
    
    G.signif <- t(h.pvals[-1,,drop=FALSE] < 1-(1-control$SA.phase1.max.p)^(p*q))
    G.signif[is.na(G.signif)] <- FALSE

    ## Compute the variances (robustly) and the statistic weights.
    v <- matrix(NA, q,q)
    v[!bad.fits,!bad.fits] <- covMcd(h.resid[,!bad.fits,drop=FALSE])$cov
    v[is.na(v)] <- 0

    w <- robust.inverse(v)
    
    ## Adjust the number of time steps between jumps.
    edge.ages <- state$nw%n%"time"-ergm.el.lasttoggle(state$nw)[,3]+1
    
    control$SA.interval<- min(control$SA.max.interval, max(control$SA.min.interval, if(length(edge.ages)>0) control$SA.interval.mul*mean(edge.ages)))
    if(verbose>1){
      cat("New interval:",control$SA.interval ,"\n")
    }
    
    ## Detect parameters whose effect we aren't able to reliably detect.
    ineffectual.pars <- !apply(G.signif,2,any)

    if(all(ineffectual.pars)){
      cat("None of the parameters have a detectable effect. Increasing jitter.\n" )
      control$jitter[!offsets] <- control$jitter[!offsets]*2
    }
    else if(any(ineffectual.pars)){
      cat("Parameters",paste.and(p.names[!offsets][ineffectual.pars]),"do not have a detectable effect. Shifting jitter to them.\n" )
      control$jitter[!offsets] <- control$jitter[!offsets] * (ineffectual.pars+1/2) / mean(control$jitter[!offsets] * (ineffectual.pars+1/2))
    }

    ## Evaluate the dstat/dpar gradient matrix.
    G <- t(h.fit[-1,,drop=FALSE])
    G[!G.signif] <- 0
    G[is.na(G)] <- 0

    rownames(w)<-colnames(w)<-rownames(v)<-colnames(v)<-q.names
    
    colnames(G)<-p.names[!offsets]
    rownames(G)<-q.names
    if(verbose>1){
      cat("Most recent parameters:\n")
      cat("Formation:\n")
      print(state$eta.form)
      cat("Dissolution:\n")
      print(state$eta.diss)
      cat("Target differences (most recent):\n")
      print(state$nw.diff)
      cat("Target differences (last run):\n")
      print(colMeans(oh.last[,-(1:p),drop=FALSE]))
      cat("Approximate objective function (most recent):\n")
      print(mahalanobis(oh[nrow(oh),-(1:p),drop=FALSE],0,cov=w,inverted=TRUE))
      cat("Approximate objective function (last run):\n")
      print(mahalanobis(colMeans(oh.last[,-(1:p),drop=FALSE]),0,cov=w,inverted=TRUE))
      cat("Estimated gradient:\n")
      print(G)
      cat("Estimated covariance of statistics:\n")
      print(v)
    }

    # Plot if requested.
    if(control$SA.plot.stats){
      if(!dev.interactive(TRUE)){
        warning("Progress plot requested on a non-interactive graphics device. Ignoring.")
        control$SA.plot.stats <- FALSE # So that we don't print a warning every step.
      }else{
        library(lattice)

        get.dev("gradients")
        G.scl <- sweep(G, 1, apply(G, 1, function(x) sqrt(mean(x^2))), "/")
        G.scl[is.nan(G.scl)] <- 0
        G.scl <- sweep(G.scl, 2, apply(G.scl, 2, function(x) sqrt(mean(x^2))), "/")
        G.scl[is.nan(G.scl)] <- 0
        try(print(.my.levelplot(G.scl,main="Scaled Gradients")),silent=TRUE)

        get.dev("correlations")
        suppressWarnings(try(print(.my.levelplot(cov2cor(v),main="Correlations")),silent=TRUE))

      }
    }

    par.eff <- apply(sweep(G,1,sqrt(diag(v)),"/"),2,function(z)sum(z^2))
    par.eff <- sqrt(par.eff)
    
    
    control$GainM <- matrix(0, nrow=p, ncol=q)
    control$GainM[!offsets,] <- t(sweep(G,2,par.eff,"/")) %*% w * control$gain
    control$GainM[!is.finite(control$GainM)] <- 0

    control$dejitter <- matrix(0, nrow=p, ncol=p)
    control$dejitter[!offsets,!offsets] <- control$GainM[!offsets,,drop=FALSE]%*%G

    rownames(control$GainM) <- rownames(control$dejitter) <- colnames(control$dejitter) <- p.names
    colnames(control$GainM) <- q.names
    
    if(verbose>1){
      cat("New deviation -> coefficient map:\n")
      print(control$GainM)
      cat("New jitter cancelation matrix:\n")
      print(control$dejitter)
    }

    if(update.jitter){
      control$jitter[!offsets] <- apply(oh[,1:p,drop=FALSE][,!offsets,drop=FALSE]-jitters[,!offsets,drop=FALSE],2,sd)*control$SA.phase2.jitter.mul
      names(control$jitter) <- p.names
    }
    
    if(verbose>1){
      cat("New jitter values:\n")
      print(control$jitter)
    }

    control$dev.guard <- apply(h[,-(1:p),drop=FALSE],2,function(x) quantile(abs(x),.9)) * control$SA.guard.mul
    if(verbose>1){
      cat("New deviation guard values:\n")
      print(control$dev.guard)
    }

    control$par.guard <- apply(abs(diff(oh[,1:p,drop=FALSE],lag=control$SA.runlength-1)),2,median) * control$SA.guard.mul
    if(verbose>1){
      cat("New parameter guard values:\n")
      print(control$par.guard)
    }
    
    list(control=control,
         G=G, w=w, v=v, oh.fit=h.fit, ineffectual.pars=ineffectual.pars, bad.fits=bad.fits, state=state)
  }

    stergm.EGMME.SA(theta.form0, theta.diss0, nw, model.form, model.diss, model.mon,
                            control, MHproposal.form, MHproposal.diss, eval.optpars,
                            verbose)
}
