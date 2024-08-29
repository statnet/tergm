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
