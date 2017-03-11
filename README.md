# Particle Swarm Optimization

This repo contains code for PSO and its application to a modeling problem. Once you've created `libpso.a`, you can use it in your own projects. If you cannot count on dynamic linkage to a newish version of GSL, you can statically link it for a more robust binary. Also provided are two RNG modules (one of which must also be compiled in). The header files contain lots of useful documentation. The file `model.c` is not intended to be reused, but rather it serves to demonstrate the PSO library in action.

To compile:

    gcc -std=c99 -O2 -c {pso,transform,util}.c
Add the flag `-DEXCLUDE_LINUX` to remove dependence on the `getrandom()` syscall.

    ar rcs libpso.a *.o
    gcc -std=c99 -O2 -DNT=<number of cores> -c {model,xorshift}.c
Reproducibility from a deterministic generator is not guaranteed for `NT > 1`. You'll have to tweak `model.c` to make `urandom` work and `pso.c` if you want a custom RNG instead.

    gcc -L. -o model {model,xorshift}.o -l{gsl,gslcblas,pso,m} -pthread
