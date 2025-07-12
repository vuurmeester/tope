#pragma once

/* alloca: */
#ifdef _WIN32
#include <malloc.h>
#else
#include <alloca.h>
#endif

/** Swap memory content. */
void memswp(void* ptr1, void* ptr2, int numbytes);

/** Random number between 0 and 1. */
double random_getdouble();

/** Time difference [s] since last call. */
double clock_gettimediff();
