
# evaluate a summary statistic formula at time points defined by at
summary.statistics.networkDynamic <- function(object, at,..., basis=NULL){
  if (missing(at) || !is.numeric(at)){
    stop( "summary.statistics.networkDynamic requires an 'at' parameter specifying the time points at which to summarize the formula")
  }
  if(!is.null(basis)) object <- ergm.update.formula(object, basis~.)
  duration.dependent <- is.durational(object)
  t(rbind(sapply(at,
              function(t){
                nw <- network.extract.with.lasttoggle(ergm.getnetwork(object), t, duration.dependent)
                f <- ergm.update.formula(object, nw~., from.new="nw")
                # make sure this is dispatched to the .network  and not .networkDynamic version
                # of summary statistics to avoid recurisve calls
                summary.statistics.network(f,...)
              }
          )
      )
  )
}
