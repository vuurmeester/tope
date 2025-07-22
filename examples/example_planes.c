#include <stdio.h>
#include <stdlib.h>

#include "../src/math.h"
#include "../src/util.h"

#include <polytoop.h>



void benchmark(int ntests, int nplanes, int ndims)
{
  random_reset();

  double start = polytoop_gettime();
  int nfacets = 0;
  int nverts = 0;
  int nridges = 0;
  for (int itest = 0; itest < ntests; ++itest) {
    /* Random center: */
    double* center = malloc(ndims * sizeof(double));
    for (int idim = 0; idim < ndims; ++idim) {
      center[idim] = 10.0 * (random_getdouble() - 0.5);
    }

    /* Random planes: */
    double* normals = malloc(nplanes * ndims * sizeof(double));
    double* dists = malloc(nplanes * sizeof(double));
    for (int iplane = 0; iplane < nplanes; ++iplane) {
      for (int idim = 0; idim < ndims; ++idim) {
        normals[iplane * ndims + idim] = random_getdouble() - 0.5;
      }
      vec_normalize(ndims, &normals[iplane * ndims]);
      dists[iplane] = 0.1 + random_getdouble() +
                      vec_dot(ndims, center, &normals[iplane * ndims]);
    }

    /* Create polytoop object: */
    Polytoop* polytoop = polytoop_fromplanes(nplanes, ndims, normals, dists);

    /* Accumulate total number of vertices created: */
    if (polytoop) {
      nfacets += polytoop_getnumfacets(polytoop);
      nverts += polytoop_getnumvertices(polytoop);
      nridges += polytoop_getnumridges(polytoop);
      polytoop_delete(polytoop);
    }

    free(dists);
    free(normals);
    free(center);
  }
  double dt = polytoop_gettime() - start;
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
