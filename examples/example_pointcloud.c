#include <stdio.h>
#include <stdlib.h>

#include "../src/math.h"
#include "../src/util.h"

#include <tope.h>



double benchmark(int ntests, int npoints, int ndims, bool cospherical)
{
  double* points = malloc(npoints * ndims * sizeof(double));

  random_reset();

  printf("***** tope test *****\n");
  printf("ntests      = %d\n", ntests);
  printf("npoints     = %d\n", npoints);
  printf("ndims       = %d\n", ndims);
  printf("cospherical = %s\n", cospherical ? "true" : "false");

  /* Random points: */
  for (int ipoint = 0; ipoint < npoints; ++ipoint) {
    for (int idim = 0; idim < ndims; ++idim) {
      points[ipoint * ndims + idim] = random_getdouble() - 0.5;
    }
    if (cospherical) {
      vec_normalize(ndims, &points[ipoint * ndims]);
    }
  }

  double start = tope_gettime();
  int nfacets = 0;
  int nridges = 0;
  int nverts = 0;
  for (int itest = 0; itest < ntests; ++itest) {
    /* Add points to tope: */
    Tope* tope = tope_frompoints(npoints, ndims, points);

    /* Accumulate total number of facets and vertices created: */
    nfacets += tope_getnumfacets(tope);
    nridges += tope_getnumridges(tope);
    nverts += tope_getnumvertices(tope);

    /* Clean up: */
    tope_delete(tope);
  }
  double dt = tope_gettime() - start;
  printf("facets      = %d\n", nfacets);
  printf("ridges      = %d\n", nridges);
  printf("verts       = %d\n", nverts);
  printf("time        = %g\n\n", dt);

  free(points);

  return dt;
}



int main(void)
{
  double t = 0.0;

  /* 4D non-cospherical: */
  t += benchmark(200, 32, 4, false);
  t += benchmark(200, 64, 4, false);

  /* 3D non-cospherical: */
  t += benchmark(200, 64, 3, false);
  t += benchmark(200, 128, 3, false);

  /* 3D cospherical: */
  t += benchmark(200, 64, 3, true);
  t += benchmark(200, 128, 3, true);

  /* Many small hulls: */
  t += benchmark(2000, 16, 3, true);

  /* One large cospherical tope: */
  t += benchmark(1, 20000, 3, true);

  /* Large 4D tope: */
  t += benchmark(1, 200000, 4, false);

  /* Large 5D tope: */
  t += benchmark(1, 50000, 5, false);

  printf("t_total = %g\n", t);

  return 0;
}
