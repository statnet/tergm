#  File man-roxygen/ergmTerm-emptyval.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2022 Statnet Commons
################################################################################
#' @param emptyval can be used to specify the value returned if the network does not have any actors
#'   with degree in the specified range. This is, technically, an arbitrary value, but it should
#'   not have a substantial effect unless a non-negligible fraction of
#'   networks at the parameter configuration of interest has no actors
#'   with specified degree.
