#include <math.h>
#include <string.h>

#include "util.h"

void util_list_shuffle(RNG_STATE_T state, uint64_t *list, size_t len)
{
    for (size_t i = len - 1; i > 0; --i)
    {
        uint64_t j = transform_integer(state, 0, i);

        uint64_t t = list[i];
        list[i] = list[j];
        list[j] = t;
    }
}

void util_list_map(double *in, double *out, double *m, double *b, size_t d)
{
    for (size_t i = 0; i < d; ++i)
        out[i] = fma(in[i], m[i], b[i]);
}

double util_list_dist(double *v, double *w, size_t d)
{
    double dist = 0;

    for (size_t i = 0; i < d; ++i)
    {
        double diff = v[i] - w[i];

        dist += diff * diff;
    }

    return sqrt(dist);
}

bool util_array_lhs(RNG_STATE_T state, double *coords, uint64_t n, size_t d)
{
    if (n == 0 || d == 0)
        return false;

    size_t num_entries = n * d;

    uint64_t *array = malloc(num_entries * 8);

    if (!array)
        return false;

    double eps = (double)1 / n;

    /*
       Initialize row-major array where the columns are grid coordinates for
       each particle. They start out as <0, 0, ...>, <1, 1, ...>, etc.
    */

    for (size_t i = 0; i < n; ++i)
        array[i] = i;

    for (size_t i = n; i < num_entries; i += n)
        memcpy(array + i, array, n * 8);

    /*
       Shuffle the rows of the array. This places the particles randomly in the
       grid while ensuring that no particle will share a "hyper-row" with
       another particle.
    */

    for (size_t i = 0; i < num_entries; i += n)
        util_list_shuffle(state, array + i, n);

    /*
       Convert columns into real coordinates. We want each particle to be
       randomly placed within its grid cube.
    */

    for (size_t i = 0; i < n; ++i)
    {
        for (size_t j = 0; j < d; ++j)
            coords[i * d + j] = fma(
                    array[j * n + i],
                    eps,
                    transform_real(state, 0, eps)
                    );
    }

    free(array);

    return true;
}
