mkStartupMessage <- function(pkgname){
  INST_MAP <- list(washington.edu="University of Washington",
                   uw.edu="University of Washington",
                   psu.edu="Penn Statue University",
                   uci.edu="University of California -- Irvine",
                   ucla.edu="University of California -- Los Angeles",
                   nyu.edu="New York University") 

  desc <- packageDescription(pkgname)
  pns <- eval(parse(text=desc$`Authors@R`))
  pnnames <- format(pns, include=c("given","family"))
  pninsts <- sapply(pns, function(pn) NVL(INST_MAP[[gsub(".*?([^.@]+\\.[^.]{2,4})$","\\1",pn$email)]],""))

  authors <- sapply(pns, function(pn) "aut" %in% pn$role)

  pnlines <- ifelse(pninsts=="", pnnames, paste(pnnames,pninsts, sep=", "))
  
  copylist <- paste("Copyright (c) ",substr(desc$Date,1,4),", ",sep="")
  copylist <- paste(copylist, pnlines[authors][1],"\n",
                    paste(
                      paste(rep(" ",nchar(copylist)),collapse=""),
                      c(pnlines[authors][-1],if(sum(!authors)) "with contributions from",pnlines[!authors]),sep="",collapse="\n"),
                    sep="") 
     paste("\n",desc$Package,": version ", desc$Version, ', created on ', desc$Date, '\n',copylist,"\n",
          'Based on "statnet" project software (statnet.org).\n',
          'For license and citation information see statnet.org/attribution\n',
          'or type citation("',desc$Package,'").\n', sep="")
}

.onAttach <- function(lib, pkg){
  packageStartupMessage(mkStartupMessage("tergm"))
  
  .RegisterMHPs()
  .RegisterConstraintImplications()
}

.RegisterMHPs <- function(){
  ergm.MHP.table("c", "Bernoulli", "atleast",  0, "random", "formationMLE")
  ergm.MHP.table("c", "Bernoulli", "atleast+bd",  0, "random", "formationMLE")
  ergm.MHP.table("c", "Bernoulli", "atleast",  1, "TNT", "formationMLETNT")
  ergm.MHP.table("c", "Bernoulli", "atleast+bd",  1, "TNT", "formationMLETNT")
  ergm.MHP.table("c", "Bernoulli", "atmost",  0, "random", "dissolutionMLE")
  ergm.MHP.table("c", "Bernoulli", "atmost+bd",  0, "random", "dissolutionMLE")
  ergm.MHP.table("c", "Bernoulli", "atleast+observed",  0, "random", "formationNonObservedMLE")
  ergm.MHP.table("c", "Bernoulli", "atmost+observed",  0, "random", "dissolutionNonObservedMLE")
  ergm.MHP.table("f", "Bernoulli", "",  0, "random", "formation")
  ergm.MHP.table("f", "Bernoulli", "bd",  0, "random", "formation")
  ergm.MHP.table("f", "Bernoulli", "",  1, "TNT", "formationTNT")
  ergm.MHP.table("f", "Bernoulli", "bd",  1, "TNT", "formationTNT")
  ergm.MHP.table("d", "Bernoulli", "",  0, "random", "dissolution")
  ergm.MHP.table("d", "Bernoulli", "bd",  0, "random", "dissolution")
}

.RegisterConstraintImplications <- function(){
  ergm.ConstraintImplications("atleast", c())
  ergm.ConstraintImplications("atmost", c())
}
