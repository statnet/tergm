InitErgmProposal.staticDiscordTNT <- function(arguments, nw, model) {
  proposal <- list(name = "staticDiscordTNT", pakcage="ergm", iinputs=to_ergm_Cdouble(arguments$constraints$discordTNT$nw))
  proposal
}
