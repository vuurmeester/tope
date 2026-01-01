#ifdef _WIN32
#include <windows.h>
#else
#define _POSIX_C_SOURCE 199309L
#include <sys/time.h>
#include <time.h>
#endif

#include <assert.h>
#include <limits.h>
#include <string.h>

#include "util.h"

static unsigned s_seed = 0;
static double const s_nrmfac = 1.0 / ((double)UINT_MAX + 1.0);



void memswp(void* ptr1, void* ptr2, int numbytes)
{
  assert(
      (char*)ptr1 + numbytes <= (char*)ptr2 ||
      (char*)ptr2 + numbytes <= (char*)ptr1
  );  /* non-overlapping */
  void* tmp = alloca(numbytes);
  memcpy(tmp, ptr1, numbytes);
  memcpy(ptr1, ptr2, numbytes);
  memcpy(ptr2, tmp, numbytes);
}



void random_reset(void)
{
  s_seed = 0;
}



double random_getdouble(void)
{
  /* Returns double x, where 0.0 <= x < 1.0: */
  s_seed = 1664525U * s_seed + 1013904223U;
  return s_nrmfac * (double)s_seed;
}



int random_getint(int lo, int hi)
{
  assert(lo < hi);
  /* Returns int n, where lo <= n < hi: */
  s_seed = 1664525U * s_seed + 1013904223U;
  int n = s_seed % (hi - lo) + lo;
  return n;
}



#ifdef _WIN32
static double const maxdouble = (double)UINT_MAX + 1.0;
static double large2double(LARGE_INTEGER const* i)
{
  return (double)i->LowPart + maxdouble * (double)i->HighPart;
}
#else
static double const maxdouble = 1.0e6;
#endif



static void getcounts(unsigned* hi, unsigned* lo)
{
#ifdef _WIN32
  LARGE_INTEGER temp;
  QueryPerformanceCounter(&temp);
  *hi = temp.HighPart;
  *lo = temp.LowPart;
#else
  struct timeval tv;
  gettimeofday(&tv, 0);
  *hi = tv.tv_sec;
  *lo = tv.tv_usec;
#endif
}



static double getperiod(void)
{
#ifdef _WIN32
  LARGE_INTEGER temp;
  QueryPerformanceFrequency(&temp);
  double period = 1.0 / large2double(&temp);
  return period;
#else
  return 1.0e-6;
#endif
}



double tope_gettime(void)
{
  unsigned hi;
  unsigned lo;
  double period = getperiod();
  getcounts(&hi, &lo);
  return (maxdouble * (double)hi + (double)lo) * period;
}
