#ifndef _RNG_H
#define _RNG_H

/*
   This is the interface that a random number generator supplying PSO must
   abide by. Disregarding the details of initialization, its function is to
   supply entropy to the algorithm in 64-bit blocks. The functions outlined in
   *transform.h* take an instance of such a generator as an argument and return
   values drawn from the distributions needed by PSO.
*/

#include <stdint.h>

typedef void *RNG_STATE_T;
typedef void *RNG_SEED_T;

/*
   This function should return a human-readable name identifying the RNG. It
   doesn't have to be unique and shouldn't be used for detecting its type in
   code.
*/

char *rng_name(void);

/*
   This function should return an identifier of the implementation that's
   globally unique. It doesn't have to be human-readable and should be used in
   code to determine how the seed should be created. We recommend a 128-bit
   hexadecimal value generated from a TRNG.
*/

char *rng_uid(void);

/*
   This function allocates a new instance of the opaque state structure needed
   by the RNG. It should return NULL on failure. The state is then passed to
   (and updated by) every subsequent call to the RNG. It should always be freed
   by the corresponding *rng_free_state()* function.
*/

RNG_STATE_T rng_allocate_state(void);

/*
   This functions frees any resources used by a previously allocated state. It
   should be harmless when called on NULL, but no guarantees are made
   regarding previosly freed or otherwise invalid states. Therefore, the
   pointer passed to it should be subsequently set to NULL to avoid common
   errors.
*/

void rng_free_state(RNG_STATE_T state);

/*
   This function should prepare the state for use by the RNG, taking into
   account the seed data provided. The internal structure of *state* never
   needs to be known to the user of these functions, but they will have to use
   the proper *seed* format as specified in the documentation of the particular
   RNG. It is possible that *seed* will be NULL for some generators.
*/

void rng_initialize_state(RNG_STATE_T state, RNG_SEED_T seed);

/*
   This function deterministically derives a seed from a human-readable,
   null-terminated value and writes it to *seed*. It leaves *seed* unchanged if
   *phrase* isn't properly formatted or seed derivation isn't implemented for
   this generator (or if either argument is NULL). However, all these conditions
   can be detected in advance and should not occur in production code.
*/

void rng_derive_seed(RNG_SEED_T seed, char *phrase);

/*
   This function returns a 64-bit value (which should ideally be drawn from a
   uniform distribution) and updates the internal state.
*/

uint64_t rng_next_block(RNG_STATE_T state);

#endif
