#ifndef _UTIL_H
#define _UTIL_H

// This file provides definitions for various PSO utility functions.

#include <stdbool.h>

#include "transform.h"

// This function randomly permutes a list of integers.

void util_list_shuffle(RNG_STATE_T state, uint64_t *list, size_t len);

/*
   This function reads in a *d*-vector v from *in* and writes Mv + b to *out*,
   where M is a diagonal matrix formed from *m* and b is a vector formed from
   *b*.
*/

void util_list_map(double *in, double *out, double *m, double *b, size_t d);

/*
   This function returns the Euclidean distance between the *d*-vectors *v* and
   *w*.
*/

double util_list_dist(double *v, double *w, size_t d);

/*
   This function writes out coordinates in the hypercube [0, 1]^*d* according
   to Latin Hypercube Sampling. The space is divided into a grid of cubes with
   side length 1 / *n*, and no position will share the same grid coordinate
   with any other position. The ith position will start at an offset of
   *d* x *i* in *coords*, where i ranges from 0 to *n* - 1. It will return
   false if either integer parameter is 0 or on a memory allocation error.
*/

bool util_array_lhs(RNG_STATE_T state, double *coords, uint64_t n, size_t d);

#endif
