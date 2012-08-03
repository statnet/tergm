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
