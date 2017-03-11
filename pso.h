#ifndef _PSO_H
#define _PSO_H

// This file provides definitions for the main PSO algorithm.

#include <stdbool.h>

#include "transform.h"

/*
   This constant defines the maximum value of K allowed. It should never need
   to be anywhere as high in practice.
*/

#ifndef PSO_MAX_NEIGHBORS
#define PSO_MAX_NEIGHBORS 25
#endif

// This constant defines the maximum PSO swarm size.

#ifndef PSO_MAX_SWARM_SIZE
#define PSO_MAX_SWARM_SIZE 50
#endif

// This definition is for the fitness function that will be supplied to PSO.

typedef double (*PSO_FITNESS_T)(double *pos);

typedef struct
{
    double pos[TRANSFORM_MAX_DIM];

    double fitness;
} PSO_RESULTS_T;

typedef struct
{
    double x[TRANSFORM_MAX_DIM];

    double tmp[TRANSFORM_MAX_DIM];

    double v[TRANSFORM_MAX_DIM];

    double p[TRANSFORM_MAX_DIM];

    double q;

    uint64_t N[PSO_MAX_NEIGHBORS + 1];

    double l[TRANSFORM_MAX_DIM];

    double m;
} PSO_PARTICLE_T;

typedef struct
{
    RNG_STATE_T state;

    pthread_mutex_t mutex;

    PSO_FITNESS_T fitness;

    size_t dim;

    size_t size;

    size_t max_evals;

    size_t k;

    double best_fitness;

    double omega;

    double c;

    PSO_PARTICLE_T particles[PSO_MAX_SWARM_SIZE];

    uint64_t indices[PSO_MAX_SWARM_SIZE];

    double lower[TRANSFORM_MAX_DIM];

    double coefs[TRANSFORM_MAX_DIM];

    double best_pos[TRANSFORM_MAX_DIM];
} PSO_SWARM_T;

/*
   This function initializes the structure that keeps track of the PSO swarm.
   It should be provided a pointer to the swarm, the fitness function, the
   lower and upper parameter bounds, the dimension of the search space, the
   swarm size, and a phrase to initialize the RNG (which can be NULL if not
   applicable). It will return true on success and false on invalid parameters
   or a memory allocation error.
*/

bool pso_initialize(
        PSO_SWARM_T *swarm,
        PSO_FITNESS_T fitness,
        double c,
        double omega,
        double *lower,
        double *upper,
        size_t dim,
        size_t size,
        size_t max_evals,
        size_t k,
        char *phrase
        );

/*
   This function computes the fitness of a given position within the hypercube
   by applying the appropriate affine transform before sending the coordinates
   to the fitness function. It uses *tmp* for temporary storage, and both
   arrays should be large enough to hold *swarm*->dim doubles.
*/

double pso_compute_fitness(PSO_SWARM_T *swarm, double *pos, double *tmp);

/*
   This function shuffles the list of particles in the swarm. It should be
   called each iteration before any other computations are performed.
*/

void pso_shuffle(PSO_SWARM_T *swarm);

/*
   This function performs the intermediate evaluation steps on particles with
   indices in the interval [*begin*, *end*]. If the fitness function is
   thread-safe and this function is called on disjoint intervals, then it is
   thread-safe too. It is where the bulk of the work takes place.
*/

void pso_evaluate_interval(PSO_SWARM_T *swarm, size_t begin, size_t end);

/*
   This function shoud be called after each interval in a partition of the
   swarm has been evaluated. It returns true if the swarm is ready for another
   iteration and false if the computation has terminated.
*/

bool pso_finalize(PSO_SWARM_T *swarm);

/*
   This function writes the current best position and corresponding fitness
   function value to the *results* structure.
*/

void pso_write_optimum(PSO_SWARM_T *swarm, PSO_RESULTS_T *results);

/*
   This function frees all the memory held by an initialized swarm. Note that
   it does not free the *swarm* structure itself, as it is not necessarily
   dynamically allocated.
*/

void pso_free(PSO_SWARM_T *swarm);

#endif
