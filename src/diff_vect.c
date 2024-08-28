#include <stdlib.h>
#include "diff_vect.h"

DiffVec diff_new(int capacity) {
    DiffVec diff = {
        .size = 0,
        .capacity = capacity,
        .content = malloc(capacity * sizeof(int)),
    };
    return diff;
}

void diff_clear(DiffVec* diff) {
    free(diff->content);
}

void diff_append(DiffVec* diff, int value) {
    if (diff->size == diff->capacity) {
        diff->capacity *= 2;
        diff->content = realloc(diff->content, diff->capacity * sizeof(int));
    }
    diff->content[diff->size++] = value;
}

SEXP diff_to_SEXP(DiffVec* diff) {
    SEXP out = PROTECT(allocVector(INTSXP, diff->size));
    int *out_ptr = INTEGER(out);
    for (int i = 0; i < diff->size; i++) {
        out_ptr[i] = diff->content[i];
    }
    UNPROTECT(1);
    return out;
}
