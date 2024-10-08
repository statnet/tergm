/*  File src/kvint.h in package tergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2024 Statnet Commons
 */
#ifndef INCLUDE_SRC_KVINT_H_
#define INCLUDE_SRC_KVINT_H_

#include <Rinternals.h>
#include "ergm_kvec.h"

typedef kvec_t(int) kvint;
SEXP kvint_to_SEXP(kvint v);

#endif  // INCLUDE_SRC_KVINT_H_
