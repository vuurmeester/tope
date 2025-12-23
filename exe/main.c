#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include "../src/math.h"
#include "../src/util.h"

#include <tope.h>

static int linenr = 0;

static void error(char const* str, ...)
{
  va_list args;
  va_start(args, str);
  char fmt[256];
  sprintf(fmt, "error on line %d: ", linenr);
  strcat(fmt, str);
  strcat(fmt, "\n");
  vprintf(fmt, args);
  va_end(args);
  exit(1);
}



static char* read_line()
{
  #define BUFSIZE 1024
  static char line[BUFSIZE];
  ++linenr;
  if (fgets(line, BUFSIZE, stdin) == NULL) {
    error("could not read line");
  }
  return line;
}



int main(void)
{
  char* start = NULL;
  char* next = NULL;

  /* Read dimension: */
  start = read_line();
  int d = strtol(start, &next, 10);
  if (next == start) {
    error("could not read dimension");
  }
  if (d < 2) {
    error("dimension must be greater than 1");
  }
  printf("dimension %d\n", d);

  /* Read point count: */
  start = read_line();
  int npoints = strtol(start, &next, 10);
  if (next == start) {
    error("could not read point count");
  }
  printf("npoints %d\n", npoints);

  double* points = malloc(d * npoints * sizeof(double));
  for (int i = 0; i < npoints; ++i) {
    /* Read vertex: */
    start = read_line();
    for (int j = 0; j < d; ++j) {
      /* Read point component: */
      points[i * d + j] = strtod(start, &next);
      if (next == start) {
        error("could not read point %d component %d", i, j);
      }
      start = next;
    }
  }

  if (npoints < 20 && d < 4) {
    mat_print(npoints, d, points);
  }
  
  /* Add points to tope: */
  double starttime = tope_gettime();
  Tope* tope = tope_frompoints(npoints, d, points);
  double endtime = tope_gettime();
  int nfacets = tope_getnumfacets(tope);
  int nridges = tope_getnumridges(tope);
  int nverts = tope_getnumvertices(tope);
  printf("facets      = %d\n"  , nfacets);
  printf("ridges      = %d\n"  , nridges);
  printf("verts       = %d\n"  , nverts);
  printf("time        = %g\n", endtime - starttime);
  printf("memory      = %llu\n\n", tope_bytes_used(tope));

  tope_delete(tope);
  free(points);

  return 0;
}
