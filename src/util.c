#include <limits.h>
#include <string.h>
#ifdef _WIN32
#include <windows.h>
#else
#define _POSIX_C_SOURCE 199309L
#include <sys/time.h>
#include <time.h>
#endif

#include "util.h"

static unsigned s_seed = 0;
static double const s_nrmfac = 1.0 / ((double)UINT_MAX + 1.0);
static double s_time = 0.0;



void memswp(void* ptr1, void* ptr2, int numbytes)
{
  void* tmp = alloca(numbytes);
  memcpy(tmp, ptr1, numbytes);
  memcpy(ptr1, ptr2, numbytes);
  memcpy(ptr2, tmp, numbytes);
}



unsigned random_getuint()
{
  s_seed = 1664525U * s_seed + 1013904223U;
  return s_seed;
}



double random_getdouble()
{
  /* Returns double x, where 0.0 <= x < 1.0: */
  return s_nrmfac * (double)random_getuint();
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



double polytoop_clock_getperiod()
{
#ifdef _WIN32
  LARGE_INTEGER temp;
  double period;
  QueryPerformanceFrequency(&temp);
  period = 1.0 / large2double(&temp);
  return period;
#else
  return 1.0e-6;
#endif
}



double polytoop_clock_gettime()
{
  unsigned hi;
  unsigned lo;
  double period = polytoop_clock_getperiod();
  getcounts(&hi, &lo);
  return (maxdouble * (double)hi + (double)lo) * period;
}



double polytoop_clock_sleep(double sleeptime)
{
  /* Initialize sleep parameters: */
  double starttime = polytoop_clock_gettime();
  double endtime = starttime + sleeptime;
  double timenow = starttime;

  /* Short sleeps: */
  while (timenow < endtime) {
#ifdef _WIN32
    Sleep(0);
#else
    struct timespec ts;
    ts.tv_sec = 0;
    ts.tv_nsec = 0;
    nanosleep(&ts, NULL);
#endif
    timenow = polytoop_clock_gettime();
  }

  /* Return actual time slept: */
  return timenow - starttime;
}



unsigned clock_getticks()
{
  unsigned hi;
  unsigned lo;

  getcounts(&hi, &lo);
  return lo;
}



double polytoop_clock_gettimediff()
{
  double newTime = polytoop_clock_gettime();
  double diff = newTime - s_time;
  s_time = newTime;
  return diff;
}
