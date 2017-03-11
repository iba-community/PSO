// This file gives an implementation of the xorshift128+ generator.

#include <stdlib.h>
#include <string.h>

#include "rng.h"

char *rng_name(void)
{
    return "xorshift128+";
}

char *rng_uid(void)
{
    return "1f3a3ccab4d1cc0447e2f8c07f35cce7";
}

RNG_STATE_T rng_allocate_state(void)
{
    return malloc(16);
}

void rng_free_state(RNG_STATE_T state)
{
    free(state);
}

// The seed value (as pointed to by *seed*) must be 128 bits.
void rng_initialize_state(RNG_STATE_T state, RNG_SEED_T seed)
{
    uint64_t *s = (uint64_t *)seed;

    // We can't have a seed of 0, so use some whimsical values instead.
    if (!(s[0] | s[1]))
    {
        s[0] = 0xbaddbeefbaddcafe;
        s[1] = 0xcafed00d8badf00d;
    }

    memcpy(state, seed, 16);

    // Spin the generator a bit to mitigate "bad" seeds.
    for (unsigned i = 0; i < 128; ++i)
        rng_next_block(state);
}

/* 128-bit multiplication emulated with 4 64-bit values, each representing a
   32-bit value plus carry. The high 128 bits of the product are ignored.
*/

#define XORSHIFT_MASK 0x00000000ffffffff;

static void mul128(uint64_t *a, uint64_t *b)
{
    uint64_t c[4];

    // Compute word 0 and carry once.
    c[0] = a[0] * b[0];
    c[1] = c[0] >> 32;
    c[0] &= XORSHIFT_MASK;

    // Compute word 1 and carry twice.
    c[1] += a[0] * b[1];
    c[2] = c[1] >> 32;
    c[1] &= XORSHIFT_MASK;

    c[1] += a[1] * b[0];
    c[2] += c[1] >> 32;
    c[1] &= XORSHIFT_MASK;

    // Compute word 2 and carry thrice.
    c[2] += a[0] * b[2];
    c[3] = c[2] >> 32;
    c[2] &= XORSHIFT_MASK;

    c[2] += a[1] * b[1];
    c[3] += c[2] >> 32;
    c[2] &= XORSHIFT_MASK;

    c[2] += a[2] * b[0];
    c[3] += c[2] >> 32;
    c[2] &= XORSHIFT_MASK;

    // Compute word 3 and discard carries.
    c[3] += a[0] * b[3] + a[1] * b[2] + a[2] * b[1] + a[3] * b[0];
    c[3] &= XORSHIFT_MASK;

    // Copy result back to *a*.
    memcpy(a, c, 32);
}

// Implements the FNV-1a 128-bit hash.
void rng_derive_seed(RNG_SEED_T seed, char *phrase)
{
    if (!seed || !phrase)
        return;

    uint64_t hash[] = { 0x6295c58d, 0x62b82175, 0x07bb0142, 0x6c62272e };
    uint64_t fnv_prime[] = { 0x0000013b, 0, 0x01000000, 0 };

    for (; *phrase; ++phrase)
    {
        hash[0] ^= *phrase;

        mul128(hash, fnv_prime);
    }

    uint64_t *s = (uint64_t *)seed;

    s[0] = hash[1] << 32 | hash[0];
    s[1] = hash[3] << 32 | hash[2];
}

uint64_t rng_next_block(RNG_STATE_T state)
{
    uint64_t *s = (uint64_t *)state;

    uint64_t x = s[0];
    uint64_t y = s[1];

    s[0] = y;

    x ^= (x << 23);

    return (s[1] = x ^ y ^ (x >> 18) ^ (y >> 5)) + y;
}
