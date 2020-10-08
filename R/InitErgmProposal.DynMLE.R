InitErgmProposal.staticDiscordTNT <- function(arguments, nw, model) {
  dissolvable <- as.rlebdm(arguments$constraints$discordTNT$nw)
  formable <- !dissolvable
  dissolvable <- ergm_dyadgen_select(arguments, nw, NULL, dissolvable)
  formable <- ergm_dyadgen_select(arguments, nw, NULL, formable)
  proposal <- list(name = "staticDiscordTNT", pakcage="ergm", formable=formable, dissolvable=dissolvable, inputs=0.5)
  proposal
}
