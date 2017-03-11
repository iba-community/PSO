#define _GNU_SOURCE

#include <string.h>

#ifndef EXCLUDE_LINUX
#include <linux/random.h>
#include <sys/syscall.h>
#include <unistd.h>
#endif

#include "pso.h"
#include "rng.h"
#include "util.h"

#define TEST_RANDOM_XORSHIFT_UID    "1f3a3ccab4d1cc0447e2f8c07f35cce7"
#define TEST_RANDOM_URANDOM_UID     "58fd3704b7de783c46c1e9d11f8fe3e2"

static RNG_STATE_T initialize_rng(char *phrase)
{
    RNG_STATE_T state = rng_allocate_state();

    if (!state)
        return NULL;

    char *uid = rng_uid();

    if (strcmp(uid, TEST_RANDOM_XORSHIFT_UID) == 0)
    {
        uint64_t seed[2];

        if (phrase)
            rng_derive_seed(seed, phrase);
        else
        {
#ifndef EXCLUDE_LINUX
            syscall(SYS_getrandom, seed, 16, 0);
#else
            rng_free_state(state);

            return NULL;
#endif
        }


        rng_initialize_state(state, seed);
    }
#ifndef EXCLUDE_LINUX
    else if (strcmp(uid, TEST_RANDOM_URANDOM_UID) == 0)
        rng_initialize_state(state, NULL);
#endif
    else
    {
        rng_free_state(state);

        return NULL;
    }

    return state;
}

static void generate_topology(PSO_SWARM_T *swarm)
{
    for (size_t i = 0; i < swarm->size; ++i)
    {
        PSO_PARTICLE_T *particle = swarm->particles + i;

        particle->N[swarm->k] = i;

        for (size_t j = 0; j < swarm->k; ++j)
            particle->N[j] = transform_integer(
                    swarm->state,
                    0,
                    swarm->size - 1
                    );
    }
}

static void broadcast(PSO_SWARM_T *swarm, size_t index)
{
    PSO_PARTICLE_T *particle = swarm->particles + index;

    for (size_t i = 0; i <= swarm->k; ++i)
    {
        PSO_PARTICLE_T *peer = swarm->particles + particle->N[i];

        peer->m = particle->q;

        memcpy(peer->l, particle->p, swarm->dim * sizeof(double));
    }
}

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
        )
{
    // Zero values check.
    if (!(swarm && fitness && lower && upper && dim && size && max_evals && k))
        goto pso_initialize_error_1;

    // Upper bounds check.
    if (
            dim > TRANSFORM_MAX_DIM ||
            size > PSO_MAX_SWARM_SIZE ||
            k > PSO_MAX_NEIGHBORS
       )
        goto pso_initialize_error_1;

    RNG_STATE_T state = initialize_rng(phrase);

    if (!state)
        goto pso_initialize_error_1;

    size_t len = dim * sizeof(double);

    double *coords = malloc(size * len);

    if (!coords)
        goto pso_initialize_error_2;

    // Initialize coordinates with Latin Hypercube Sampling.
    if (!util_array_lhs(state, coords, size, dim))
        goto pso_initialize_error_3;

    // Initialize constants.
    swarm->fitness = fitness;
    swarm->c = c;
    swarm->omega = omega;
    swarm->dim = dim;
    swarm->size = size;
    swarm->max_evals = max_evals;
    swarm->k = k;

    swarm->state = state;

    // Initialize affine transform parameters
    memcpy(swarm->lower, lower, len);

    for (size_t i = 0; i < dim; ++i)
        swarm->coefs[i] = upper[i] - lower[i];

    // Do per-particle initialization.
    generate_topology(swarm);

    for (size_t i = 0; i < size; ++i)
    {
        swarm->indices[i] = i;

        PSO_PARTICLE_T *particle = swarm->particles + i;

        // Copy coordinates to particle's storage.
        double *pos = coords + i * dim;

        memcpy(particle->x, pos, len);
        memcpy(particle->p, pos, len);
        memcpy(particle->l, pos, len);

        // Evaluate fitness.
        particle->q = pso_compute_fitness(swarm, particle->x, particle->tmp);
        particle->m = particle->q;

        if (i == 0 || particle->q < swarm->best_fitness)
        {
            swarm->best_fitness = particle->q;

            memcpy(swarm->best_pos, pos, len);
        }

        // Initialize velocity.
        for (size_t j = 0; j < dim; ++j)
            particle->v[j] = transform_real(
                    state,
                    -particle->x[j],
                    1 - particle->x[j]
                    );
    }

    free(coords);

    for (size_t i = 0; i < size; ++i)
    {
        PSO_PARTICLE_T *particle = swarm->particles + i;

        size_t j = 0;

        for (; j < k; ++j)
            if (particle->q >= (swarm->particles + particle->N[j])->q)
                break;

        if (j == k)
            broadcast(swarm, i);
    }

    return true;

pso_initialize_error_3:
    free(coords);
pso_initialize_error_2:
    rng_free_state(state);
pso_initialize_error_1:
    return false;
}

double pso_compute_fitness(PSO_SWARM_T *swarm, double *pos, double *tmp)
{
    util_list_map(pos, tmp, swarm->coefs, swarm->lower, swarm->dim);

    return swarm->fitness(tmp);
}

void pso_shuffle(PSO_SWARM_T *swarm)
{
    util_list_shuffle(swarm->state, swarm->indices, swarm->size);
}

void pso_evaluate_interval(PSO_SWARM_T *swarm, size_t begin, size_t end)
{
    for (size_t i = begin; i <= end; ++i)
    {
        PSO_PARTICLE_T *particle = swarm->particles + swarm->indices[i];

        if (util_list_dist(particle->p, particle->l, swarm->dim) == 0)
        {
            for (size_t j = 0; j < swarm->dim; ++j)
                particle->tmp[j] =
                    particle->x[j] +
                    swarm->c / 2 * (
                            particle->p[j] -
                            particle->x[j]
                            );
        }
        else
        {
            for (size_t j = 0; j < swarm->dim; ++j)
                particle->tmp[j] =
                    particle->x[j] +
                    swarm->c / 3 * (
                            particle->p[j] +
                            particle->l[j] -
                            2 * particle->x[j]
                            );

        }

        transform_hypersphere(
                swarm->state,
                util_list_dist(particle->x, particle->tmp, swarm->dim),
                particle->tmp,
                swarm->dim
                );

        for (size_t j = 0; j < swarm->dim; ++j)
        {
            particle->v[j] =
                swarm->omega * particle->v[j] +
                particle->tmp[j] -
                particle->x[j];

            particle->x[j] += particle->v[j];

            if (particle->x[j] < 0)
            {
                particle->x[j] = 0;
                particle->v[j] *= -0.5;
            }
            else if (particle->x[j] > 1)
            {
                particle->x[j] = 1;
                particle->v[j] *= -0.5;
            }
        }

        double fitness = pso_compute_fitness(
                swarm,
                particle->x,
                particle->tmp
                );

        if (fitness < particle->q)
        {
            memcpy(particle->p, particle->x, swarm->dim * sizeof(double));

            particle->q = fitness;
        }
    }
}

bool pso_finalize(PSO_SWARM_T *swarm)
{
    double old_fitness = swarm->best_fitness;

    for (size_t i = 0; i < swarm->size; ++i)
    {
        PSO_PARTICLE_T *particle = swarm->particles + swarm->indices[i];

        if (particle->q < swarm->best_fitness)
        {
            swarm->best_fitness = particle->q;

            memcpy(swarm->best_pos, particle->p, swarm->dim * sizeof(double));
        }

        if (particle->q < particle->m)
            broadcast(swarm, swarm->indices[i]);
    }

    if (swarm->best_fitness == old_fitness)
        generate_topology(swarm);

    if (swarm->max_evals >= swarm->size)
    {
        swarm->max_evals -= swarm->size;

        return true;
    }
    else
        return false;
}

void pso_write_optimum(PSO_SWARM_T *swarm, PSO_RESULTS_T *results)
{
    util_list_map(
            swarm->best_pos,
            results->pos,
            swarm->coefs,
            swarm->lower,
            swarm->dim
            );

    results->fitness = swarm->best_fitness;
}

void pso_free(PSO_SWARM_T *swarm)
{
    rng_free_state(swarm->state);
}
