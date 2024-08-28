#ifndef INCLUDE_SRC_DIFF_VECT_H_
#define INCLUDE_SRC_DIFF_VECT_H_

#include <Rinternals.h>

typedef struct {
    int size;
    int capacity;
    int *content;
} DiffVec;

DiffVec diff_new(int capacity);
void diff_clear(DiffVec* diff);
SEXP diff_to_SEXP(DiffVec* diff);
void diff_append(DiffVec* diff, int value);

#endif  // INCLUDE_SRC_DIFF_VECT_H_
