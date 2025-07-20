#include <stdio.h>
#include <stdlib.h>

/*#include <qhull/libqhull.h>
#include <qhull/io.h>*/

#include <stdio.h>
#include <stdlib.h>

#include "../src/math.h"
#include "../src/util.h"

#include <polytoop.h>



double benchmark(int ntests, int npoints, int ndims, int cospherical)
{
  int nfacets;
  int nridges;
  int nverts;
  int itest;
  int idim;
  int ipoint;
  double dt;
  double* points;
  Polytoop* polytoop;

  random_reset();

  points = malloc(npoints * ndims * sizeof(double));

  printf("***** polytoop test *****\n");
  printf("ntests      = %d\n", ntests);
  printf("npoints     = %d\n", npoints);
  printf("ndims       = %d\n", ndims);
  printf("cospherical = %s\n", cospherical ? "true" : "false");

  double start = clock_gettime();
  nfacets = 0;
  nridges = 0;
  nverts = 0;
  for (itest = 0; itest < ntests; ++itest) {
    /* Random points: */
    for (ipoint = 0; ipoint < npoints; ++ipoint) {
      for (idim = 0; idim < ndims; ++idim) {
        points[ipoint * ndims + idim] = random_getdouble() - 0.5;
      }
      if (cospherical) {
        vec_normalize(ndims, &points[ipoint * ndims]);
      }
    }

    /* Create polytoop object: */
    polytoop = polytoop_new();

    /* Add points to polytoop: */
    polytoop_frompoints(polytoop, npoints, ndims, points);

    /* Accumulate total number of facets and vertices created: */
    nfacets += polytoop_getnumfacets(polytoop);
    nridges += polytoop_getnumridges(polytoop);
    nverts += polytoop_getnumvertices(polytoop);

    /* Clean up: */
    polytoop_delete(polytoop);
  }
  dt = clock_gettime() - start;
  printf("facets      = %d\n", nfacets);
  printf("ridges      = %d\n", nridges);
  printf("verts       = %d\n", nverts);
  printf("time        = %g\n\n", dt);

  free(points);

  return dt;
}



void benchqhull(int ntests, int npoints, int ndims, int cospherical)
{
  int nfacets;
  int nverts;
  int itest;
  int idim;
  int ipoint;
  double dt;
  double* points;

  points = malloc(npoints * ndims * sizeof(double));

  printf("***** qhull test *****\n");
  printf("ntests      = %d\n", ntests);
  printf("npoints     = %d\n", npoints);
  printf("ndims       = %d\n", ndims);
  printf("cospherical = %s\n", cospherical ? "true" : "false");

  double start = clock_gettime();
  nfacets = 0;
  nverts = 0;
  for (itest = 0; itest < ntests; ++itest) {
    /* Random points: */
    for (ipoint = 0; ipoint < npoints; ++ipoint) {
      for (idim = 0; idim < ndims; ++idim) {
        points[ipoint * ndims + idim] = random_getdouble() - 0.5;
      }
      if (cospherical) {
        vec_normalize(ndims, &points[ipoint * ndims]);
      }
    }

    /* Call qhull: */
    /*qh_new_qhull(ndims, npoints, points, 0, "qhull Qt", 0, stderr);*/

    /* Accumulate total number of facets and vertices created: */
    /*nfacets += qh_qh.num_facets;
    nverts += qh_qh.num_vertices;*/

    /* Clean up: */
    /*qh_freeqhull(1);
    {
      int curlong;
      int totlong;
      qh_memfreeshort(&curlong, &totlong);
    }*/
  }
  dt = clock_gettime() - start;
  printf("facets      = %d\n", nfacets);
  printf("verts       = %d\n", nverts);
  printf("time        = %g\n\n", dt);

  free(points);
}



int main()
{
  // #ifndef NDEBUG
  //   benchmark(1, 64, 8, 0);
  //   benchqhull(1, 64, 8, 0);
  //   return 0;
  // #endif
  double t = 0.0;
  /* 4D non-cospherical: */
  t += benchmark(200, 32, 4, 0);
  t += benchmark(200, 64, 4, 0);

  /* 3D non-cospherical: */
  t += benchmark(200, 64, 3, 0);
  t += benchmark(200, 128, 3, 0);

  /* 3D cospherical: */
  t += benchmark(200, 64, 3, 1);
  t += benchmark(200, 128, 3, 1);

  /* Many small hulls: */
  t += benchmark(2000, 16, 3, 1);

  /* One large polytoop: */
  t += benchmark(1, 20000, 3, 1);

  /* Large 4D polytoop: */
  t += benchmark(1, 200000, 4, 0);

  /* Large 5D polytoop: */
  t += benchmark(1, 50000, 5, 0);

  printf("t_qhull = %g\n", t);

  return 0;
}
