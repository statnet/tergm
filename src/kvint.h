#ifndef INCLUDE_SRC_KVINT_H_
#define INCLUDE_SRC_KVINT_H_

#include <Rinternals.h>
#include "ergm_kvec.h"

typedef kvec_t(int) kvint;
SEXP kvint_to_SEXP(kvint v);

#endif  // INCLUDE_SRC_KVINT_H_
