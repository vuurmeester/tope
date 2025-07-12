#include <stdlib.h>
#include <stdio.h>

#include "../src/util.h"
#include "../src/math.h"

#include <polytoop.h>



void benchmark(int ntests, int nplanes, int ndims)
{
  int nfacets;
  int itest;
  int idim;
  int iplane;
  double dt;
  double* normals;
  double* dists;
  double* center;
  Polytoop* polytoop;

  clock_gettimediff();
  nfacets = 0;
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
      dists[iplane] = 0.1 + random_getdouble() + vec_dot(ndims, center, &normals[iplane * ndims]);
    }

    /* Create polytoop object: */
    polytoop = polytoop_new();

    /* Add planes to polytoop: */
    polytoop_fromplanes(polytoop, nplanes, ndims, normals, dists);

    /* Accumulate total number of vertices created: */
    nfacets += polytoop_getnumfacets(polytoop);

    /* Clean up: */
    polytoop_delete(polytoop);
    free(dists);
    free(normals);
    free(center);
  }
  dt = clock_gettimediff();
  printf("ntests      = %d\n", ntests);
  printf("nplanes     = %d\n", nplanes);
  printf("ndims       = %d\n", ndims);
  printf("facets      = %d\n", nfacets);
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
