#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>

#include <linux/random.h>
#include <sys/syscall.h>
#include <unistd.h>

#include "rng.h"

typedef struct
{
    uint64_t buf[32];

    unsigned ctr;
} URANDOM_STATE_T;

char *rng_name(void)
{
    return "/dev/urandom";
}

char *rng_uid(void)
{
    return "58fd3704b7de783c46c1e9d11f8fe3e2";
}

RNG_STATE_T rng_allocate_state(void)
{
    return malloc(sizeof(URANDOM_STATE_T));
}

void rng_free_state(RNG_STATE_T state)
{
    free(state);
}

void rng_initialize_state(RNG_STATE_T state, RNG_SEED_T seed)
{
    ((URANDOM_STATE_T *)state)->ctr = 32;
}

void rng_derive_seed(RNG_SEED_T seed, char *phrase)
{
    // Do nothing.
}

uint64_t rng_next_block(RNG_STATE_T state)
{
    URANDOM_STATE_T *s = (URANDOM_STATE_T *)state;

    if (s->ctr == 32)
    {
        syscall(SYS_getrandom, s->buf, 256, 0);

        s->ctr = 0;
    }

    return s->buf[s->ctr++];
}
