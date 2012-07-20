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
# negative, white for 0, and ensures that the scales on both sides are
# the same.

.my.levelplot <- function(m,levels=80,...){
  bound <- max(na.omit(c(abs(m))))

  levelplot(m, at=unique(c(seq(-bound,0,length.out=levels/2+1),
                 seq(0,+bound,length.out=levels/2+1))),
            col.regions=c(hsv(h=0,s=seq(1,2/levels,length.out=levels/2),v=1),
              hsv(h=.5,s=seq(2/levels,1,length.out=levels/2),v=1)),
            ...)
}
