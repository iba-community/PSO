#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <pthread.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>

#include "pso.h"

#include "nord.dat"

#ifndef NT
#define NT 1
#endif

char *names[] = 
{
    "B_0",
    "h",
    "beta_B",
    "beta_I",
    "eta",
    "gamma",
    "delta",
    "omega"
};

typedef struct
{
    PSO_SWARM_T *swarm;
    size_t begin;
    size_t end;
} JOB_T;

static int sirb(
        double time,
        const double *depvars,
        double *dydt,
        void *raw_params
        )
{
    double *params = (double *)raw_params;

    // Give everything names for instructive purposes.

    double S = depvars[0];
    double I = depvars[1];
    double R = depvars[2];
    double B = depvars[3];

    double N = S + I + R;
    double b = params[0];
    double d = params[1];
    double kappa = params[2];
    double beta_B = params[3];
    double beta_I = params[4];
    double eta = params[5];
    double gamma = params[6];
    double delta = params[7];
    double omega = params[8];

    dydt[0] = b*N - d*S - beta_B*(B*S/(kappa + B)) - beta_I*(S*I/N) + omega*R;
    dydt[1] = -d*I + beta_B*(B*S/(kappa + B)) + beta_I*(S*I/N) - gamma*I;
    dydt[2] = -d*R + gamma*I - omega*R;
    dydt[3] = eta*I - delta*B;

    return GSL_SUCCESS;
}

static bool solve(
        double *params,
        double *initial,
        double *timeline,
        size_t timeline_len,
        double *output
        )
{
    gsl_odeiv2_system system =
    {
        .function = sirb,
        .jacobian = NULL,
        .dimension = 4,
        .params = params
    };

    gsl_odeiv2_step *step = gsl_odeiv2_step_alloc(gsl_odeiv2_step_rkf45, 4);

    gsl_odeiv2_control *control = gsl_odeiv2_control_y_new(1e-6, 1e-3);

    gsl_odeiv2_evolve *evolve = gsl_odeiv2_evolve_alloc(4);

    double t = 0;
    double h = 1e-6;

    bool success = false;

    for (size_t i = 0; i < timeline_len; ++i)
    {
        double t1 = timeline[i];

        while (t < t1)
        {
            if (gsl_odeiv2_evolve_apply(
                    evolve,
                    control,
                    step,
                    &system,
                    &t,
                    t1,
                    &h,
                    initial
                    ) != GSL_SUCCESS)
                goto solve_finish;
        }

        memcpy(output + 4 * i, initial, 4 * sizeof(double));
    }

    success = true;

solve_finish:
    gsl_odeiv2_evolve_free(evolve);
    gsl_odeiv2_control_free(control);
    gsl_odeiv2_step_free(step);

    return success;
}

static double fitness(double *pos)
{
    double initial[] = { S_INIT, I_INIT, 0, 0 };
    double params[9] = { 0.000072, 0.000044, 1e6 };
    double output[4 * TIMES_LEN];

    initial[1] /= pos[1];
    initial[3] = pos[0] * 1e6;

    memcpy(params + 3, pos + 2, 6 * sizeof(double));

    if (!solve(params, initial, times, TIMES_LEN, output))
    {
        fputs("Solver error!\n", stderr);
        exit(EXIT_FAILURE);
    }

    double mad = 0;

    for (size_t i = 0; i < TIMES_LEN; ++i)
        mad += fabs(output[4 * i + 1] * pos[1] - ivals[i]);

    return mad / TIMES_LEN;
}

static void partition(PSO_SWARM_T *swarm, JOB_T *jobs)
{
    size_t per_thread = swarm->size / NT;

    size_t remainder = swarm->size % NT;

    for (size_t i = 0; i < remainder; ++i)
    {
        jobs[i].swarm = swarm;
        jobs[i].begin = i * (per_thread + 1);
        jobs[i].end = i * (per_thread + 1) + per_thread;
    }

    for (size_t i = remainder; i < NT; ++i)
    {
        jobs[i].swarm = swarm;
        jobs[i].begin = remainder + i * per_thread;
        jobs[i].end = remainder + (i + 1) * per_thread - 1;
    }
}

void *task(void *data)
{
    JOB_T *job = (JOB_T *)data;

    pso_evaluate_interval(job->swarm, job->begin, job->end);

    return NULL;
}

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        fprintf(stderr, "Usage: %s \"Seed phrase\"\n", argv[0]);

        return EXIT_FAILURE;
    }

    PSO_SWARM_T swarm;
    PSO_RESULTS_T results;

    double lower[] =
    {
        0,
        0.001,
        0.001,
        0.001,
        0.001,
        1.0 / 14,
        1.0 / 40,
        0.0001
    };

    double upper[] =
    {
        1,
        1,
        1,
        0.5,
        1,
        1.0 / 2,
        1.0 / 3,
        1.0 / 360
    };

    size_t max_evals = 2000000;

    if (!pso_initialize(
                &swarm,
                fitness,
                1.193,
                0.721,
                lower,
                upper,
                8,
                40,
                max_evals,
                3,
                argv[1]
                ))
    {
        fputs("Failed to initialize swarm!\n", stderr);

        return EXIT_FAILURE;
    }

    JOB_T jobs[NT];

    pthread_t threads[NT];

    partition(&swarm, jobs);

    do
    {
        printf(
                "\rProgress: %.0f%%",
                100 * (1 - (double)swarm.max_evals / max_evals)
              );

        fflush(stdout);

        pso_shuffle(&swarm);

        for (size_t i = 0; i < NT; ++i)
            if (pthread_create(threads + i, NULL, task, jobs + i) != 0)
            {
                fputs("Thread creation error!\n", stderr);

                return EXIT_FAILURE;
            }

        for (size_t i = 0; i < NT; ++i)
            if (pthread_join(threads[i], NULL) != 0)
            {
                fputs("Thread join error!\n", stderr);

                return EXIT_FAILURE;
            }

    } while (pso_finalize(&swarm));

    pso_write_optimum(&swarm, &results);

    printf("\nFitness: %.2f\n", results.fitness);

    for (unsigned i = 0; i < 8; ++i)
        printf("%s:\t%.6e\n", names[i], results.pos[i]);

    return EXIT_SUCCESS;
}
