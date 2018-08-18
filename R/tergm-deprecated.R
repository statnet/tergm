#' @name tergm-deprecated
#' @rdname tergm-deprecated
#' @title Functions that will no longer be supported in future releases of the package
#' @description Functions that have been superceed, were never documented, or will be removed from the package for other reasons
#' @keywords misc internal
NULL

.dep_method <- local({
  warned <- c()
  function(generic, class){
    fullname <- paste(generic,class,sep=".")
    if(! fullname%in%warned){
      me <- sys.call(-1)[[1]]
      if(length(me)>1 && me[[1]]=="::") me <- me[[3]]
      parent <- sys.call(-2)[[1]]
      if(length(parent)>1 && parent[[1]]=="::") parent <- parent[[3]]
      if(me==fullname && NVL(parent,"")!=generic){
        do.call(".Deprecated", list(msg=paste0("You appear to be calling ", fullname,"() directly. ", fullname,"() is a method, and will not be exported in a future version of ", sQuote("ergm"),". Use ", generic, "() instead, or getS3method() if absolutely necessary."), old=fullname))
        warned <<- c(warned, fullname)
      }
    }
  }
})

# Only evaluate deprecation warning once per function.
.dep_once <- local({
  warned <- c()
  function(...){
    me <- sys.call(-1)
    myname <- as.character(me[[1]])
    if(length(myname)>1 && myname[[1]]=="::") myname <- myname[[3]]
    if(! myname%in%warned){
      do.call(".Deprecated", modifyList(list(old=myname),list(...)))
      warned <<- c(warned, myname)
    }
  }
})
