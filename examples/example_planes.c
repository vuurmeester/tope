#include <stdio.h>
#include <stdlib.h>

#include "../src/math.h"
#include "../src/util.h"

#include <tope.h>



void benchmark(int ntests, int nplanes, int ndims)
{
  int itest;
  int idim;
  int iplane;
  int nfacets = 0;
  int nverts = 0;
  int nridges = 0;
  double start;
  double dt;
  double* center;
  double* normals;
  double* dists;
  Tope* tope;

  start = tope_gettime();
  random_reset();
  for (itest = 0; itest < ntests; ++itest) {
    /* Random center: */
    center = malloc(ndims * sizeof(double));
    for (idim = 0; idim < ndims; ++idim) {
      center[idim] = 10.0 * (random_getdouble() - 0.5);
    }

    /* Random planes: */
    normals = malloc(nplanes * ndims * sizeof(double));
    dists = malloc(nplanes * sizeof(double));
    for (iplane = 0; iplane < nplanes; ++iplane) {
      for (idim = 0; idim < ndims; ++idim) {
        normals[iplane * ndims + idim] = random_getdouble() - 0.5;
      }
      vec_normalize(ndims, &normals[iplane * ndims]);
      dists[iplane] = 0.1 + random_getdouble() +
                      vec_dot(ndims, center, &normals[iplane * ndims]);
    }

    /* Create tope object: */
    tope = tope_fromplanes(nplanes, ndims, normals, dists, center);

    /* Accumulate total number of vertices created: */
    if (tope) {
      nfacets += tope_getnumfacets(tope);
      nverts += tope_getnumvertices(tope);
      nridges += tope_getnumridges(tope);
      tope_delete(tope);
    }

    free(dists);
    free(normals);
    free(center);
  }
  dt = tope_gettime() - start;
  printf("ntests      = %d\n", ntests);
  printf("nplanes     = %d\n", nplanes);
  printf("ndims       = %d\n", ndims);
  printf("facets      = %d\n", nfacets);
  printf("ridges      = %d\n", nridges);
  printf("verts       = %d\n", nverts);
  printf("time        = %g\n\n", dt);
}



int main()
{
  benchmark(500, 32, 4);
  benchmark(500, 64, 4);
  benchmark(500, 64, 3);
  benchmark(500, 128, 3);

  return 0;
}
