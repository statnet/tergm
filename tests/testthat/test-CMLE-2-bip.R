#  File tests/dynamic_MLE_2.R in package tergm, part of the
#  Statnet suite of packages for network analysis, https://statnet.org .
#
#  This software is distributed under the GPL-3 license.  It is free,
#  open source, and has the attribution requirements (GPL Section 7) at
#  https://statnet.org/attribution .
#
#  Copyright 2008-2022 Statnet Commons
################################################################################

attach(CMLE.tools)
do.run_2(FALSE, m, prop.weights=c("default","random"))
detach(CMLE.tools)
