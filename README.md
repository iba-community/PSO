To compile:

    gcc -std=c99 -O2 -c {pso,transform,util}.c
Add the flag `-DEXCLUDE_LINUX` to remove dependence on the `getrandom()` syscall.

    ar rcs libpso.a *.o
    gcc -std=c99 -O2 -DNT=<number of cores> -c {model,xorshift}.c
Reproducibility from a deterministic generator is not guaranteed for `NT > 1`. You'll have to tweak `model.c` to make `urandom` work and `pso.c` if you want a custom RNG instead.

    gcc -L. -o model {model,xorshift}.o -l{gsl,gslcblas,pso,m} -pthread
