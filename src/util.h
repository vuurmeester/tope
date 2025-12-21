#pragma once

/* alloca: */
#ifdef _WIN32
#include <malloc.h>
#else
#include <alloca.h>
#endif

/** Swap memory content. */
void memswp(void* ptr1, void* ptr2, int numbytes);

/** Reset random number generator. */
void random_reset(void);

/** Random number between 0 and 1. */
double random_getdouble(void);

/** Random number in range [lo, hi>. */
int random_getint(int lo, int hi);

/** Time difference [s] since last call. */
double tope_gettime(void);
