/*  File src/kvint.c in package tergm, part of the
 *  Statnet suite of packages for network analysis, https://statnet.org .
 *
 *  This software is distributed under the GPL-3 license.  It is free,
 *  open source, and has the attribution requirements (GPL Section 7) at
 *  https://statnet.org/attribution .
 *
 *  Copyright 2008-2024 Statnet Commons
 */
#include "kvint.h"

SEXP kvint_to_SEXP(kvint v) {
    int size = kv_size(v);
    SEXP out = PROTECT(allocVector(INTSXP, size));
    int *out_ptr = INTEGER(out);
    for (int i = 0; i < size; i++) {
        out_ptr[i] = kv_A(v, i);
    }
    UNPROTECT(1);
    return out;
}
