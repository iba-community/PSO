#include <math.h>

#include <pthread.h>

#include "transform.h"

uint64_t transform_integer(RNG_STATE_T state, uint64_t l, uint64_t u)
{
    if (u <= l)
        return l;

    uint64_t diff = u - l;

    // The smallest integer such that 2^*pos* > diff.
    int pos = 8 * sizeof(long long) - __builtin_clzll(diff);

    uint64_t mask = (pos == 64) ? UINT64_MAX : (1 << pos) - 1;

    unsigned num_blocks = 64 / pos;

    unsigned rem = 64 % num_blocks;

    for (unsigned i = 0; i < TRANSFORM_MAX_TRIES; ++i)
    {
        uint64_t candidate = rng_next_block(state) >> rem;

        for (unsigned j = 0; j < num_blocks; ++j)
        {
            uint64_t res = candidate & mask;

            if (res <= diff)
                return l + res;

            candidate >>= pos;
        }
    }

    return (l + u) / 2;
}

// 2^53
#define TRANSFORM_DENOMINATOR 0x20000000000000

double transform_real(RNG_STATE_T state, double l, double u)
{
    double raw = (double)(rng_next_block(state) >> (64 - 53));

    return fma(u - l, raw / TRANSFORM_DENOMINATOR, l);
}

double transform_normal(RNG_STATE_T state, double mu, double sigma)
{
    for (unsigned i = 0 ; i < TRANSFORM_MAX_TRIES; ++i)
    {
        double u = transform_real(state, -1, 1);
        double v = transform_real(state, -1, 1);

        double s = u * u + v * v;

        if (s > 0 && s < 1)
            return fma(u * sqrt(-2 * log(s) / s), sigma, mu);
    }

    return mu;
}

/*
   Since this is the only function which is supposed to be called in a
   multithreaded environment, it includes a static mutex to make sure the RNG
   remains in a consistent state.
*/

void transform_hypersphere(RNG_STATE_T state, double r, double *c, size_t d)
{
    static pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;

    double scratch[TRANSFORM_MAX_DIM];

    d = (d > TRANSFORM_MAX_DIM) ? TRANSFORM_MAX_DIM : d;

    double factor = 0;

    pthread_mutex_lock(&mutex);

    double scale = r * pow(transform_real(state, 0, 1), (double)1 / d);

    for (size_t i = 0; i < d; ++i)
    {
        double val = transform_normal(state, 0, 1);

        scratch[i] = val;

        factor += val * val;
    }

    pthread_mutex_unlock(&mutex);

    factor = scale / sqrt(factor);

    for (size_t i = 0; i < d; ++i)
        c[i] = fma(scratch[i], factor, c[i]);
}
