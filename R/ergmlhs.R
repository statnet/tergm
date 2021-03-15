.combine_ergmlhs <- function(nwl, ignore.settings=c()){
  ergml <- nwl %>% map(`%n%`, "ergm") %>% map(NVL, list())
  l <- ergml[[1]]

  for(i in seq_along(ergml)[-1]){
    l2 <- ergml[[i]]
    all_names <- setdiff(union(names(l),names(l2)), ignore.settings)
    for(name in all_names){
      if(!identical(l[[name]], l2[[name]])){
        message("Setting ", sQuote(name), " differs between network ", i, " and some prior network and will overwrite it. Is probably benign (e.g., environment of a formula that does not reference its environment), or the networks may have inconsistent ERGM settings.")
        l[[name]] <- l2[[name]]
      }
    }
  }
  l
}
