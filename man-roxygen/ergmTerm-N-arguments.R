#  File man-roxygen/ergmTerm-N-arguments.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free, open
#  source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2025 Statnet Commons
################################################################################
#' @param lm,subset,weights,contrasts,offset,label **[NetSeries()] LHS only** arguments to specify time-varying parameters. See [`N()`][ergm.multi::N-ergmTerm] term operator in the \pkg{ergm.multi} for details. `lm` formula may reference `.Time` for the network's time index, `.TimeID` for the its index in the network series (where the initial network is 1 and the first modelled network is 2), and `.TimeDelta` for the time elapsed between the network and the immediately previous network in the series.
