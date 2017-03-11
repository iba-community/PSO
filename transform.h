#ifndef _TRANSFORM_H
#define _TRANSFORM_H

/*
   This is the interface for the supplied functions which transform a random
   bitstream into values from various distributions. Each function takes an RNG
   state as an input, possibly followed by distribution parameters.
*/

#include <stdlib.h>

#include "rng.h"

/*
   This constant defines the number of tries to use with rejection sampling
   before the function "gives up." If p = 1/2, the probability of reaching this
   condition on truly random data is 2^(-128). As Bruce Schneier puts it, this
   is significantly smaller than the chance of being killed by a meteor while
   reading this comment. Still alive? Good. Its main use is to prevent
   functions from hanging indefinitely in the case of RNG failure.
*/

#ifndef TRANSFORM_MAX_TRIES
#define TRANSFORM_MAX_TRIES 128
#endif

/*
   This constant gives the maximum number of dimensions for hypherspherical
   sampling. It saves us from having to dynamically allocate memory. This also
   places a bound on the maxiumum dimensionality of PSO.
*/

#ifndef TRANSFORM_MAX_DIM
#define TRANSFORM_MAX_DIM 50
#endif

/*
   This function returns a non-negative integer value uniformly drawn from
   [*l, *u*]. If *u* <= *l*, it will simply return *l*. It uses an integer
   base-2 logarithm to determine the appropriate interval for rejection
   sampling.
*/

uint64_t transform_integer(RNG_STATE_T state, uint64_t l, uint64_t u);

/*
   This function returns a floating point value which comes from a "good
   enough" approximately uniform distribution on [*l*, *u*). Since the
   probability of getting exactly one of the bounds is miniscule, the topology
   of the interval is largely irrelevant.
*/

double transform_real(RNG_STATE_T state, double l, double u);

/*
   This function returns a number drawn from a normal distribution with mean
   *mu* and standard deviation *sigma*.
*/

double transform_normal(RNG_STATE_T state, double mu, double sigma);

/*
   This function uniformly generates a position vector whose distance from *c*
   is less than or equal to *r*. The number of dimensions is given by *d*,
   which is clipped if it exceeds TRANSFORM_MAX_DIM. The result is written over
   *c*.
*/

void transform_hypersphere(RNG_STATE_T state, double r, double *c, size_t d);

#endif
