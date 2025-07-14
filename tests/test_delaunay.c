#include <stdlib.h>
#include <stdio.h>

#include "../src/util.h"
#include "../src/math.h"

#include <polytoop.h>

#define DIM 2
#define NPOINTS 5



int main()
{
  int m;
  int ifacet;
  int ivertex;
  int indices[DIM + 1];
  double points[NPOINTS * DIM];
  double weights[DIM + 1];
  double xi[DIM];
  double* positions;
  Polytoop* polytoop;
  polytoop_Facet* facet;

  /* Define a square. */
  m = 0;

  points[m * DIM + 0] = -2.0;
  points[m * DIM + 1] = -2.0;
  ++m;

  points[m * DIM + 0] = -2.0;
  points[m * DIM + 1] = 2.0;
  ++m;

  points[m * DIM + 0] = 2.0;
  points[m * DIM + 1] = -2.0;
  ++m;

  points[m * DIM + 0] = 2.0;
  points[m * DIM + 1] = 2.0;
  ++m;

  /* Internal point: */
  points[m * DIM + 0] = 0.0;
  points[m * DIM + 1] = 0.0;
  ++m;

  /* Create polytoop object: */
  polytoop = polytoop_new();

  /* Add points to polytoop: */
  polytoop_delaunay(polytoop, m, 2, points);

  if (polytoop_getnumfacets(polytoop) != 4) {
    return -1;
  }

  /* Interpolation: */
  xi[0] = 1;
  xi[1] = 4;
  polytoop_interpolate(polytoop, xi, indices, weights);

  /* Sum of weights should be 1: */
  if (fabs(vec_sum(DIM + 1, weights) - 1.0) > 1.0e-6) {
    return -1;
  }

  /* Print the polytoop: */
  for (facet = polytoop_firstfacet(polytoop), ifacet = 0; facet != NULL; facet = polytoop_facet_nextfacet(facet), ++ifacet) {
    printf("Facet %d:\n", ifacet + 1);
    positions = malloc(DIM * polytoop_facet_getnumvertices(facet) * sizeof(double));
    for (ivertex = 0; ivertex < polytoop_facet_getnumvertices(facet); ++ivertex) {
      polytoop_vertex_getposition(polytoop_facet_getvertex(facet, ivertex), &positions[ivertex * DIM]);
    }
    mat_print(polytoop_facet_getnumvertices(facet), DIM, positions);
    free(positions);
  }

  /* Delete polytoop object: */
  polytoop_delete(polytoop);

  return 0;
}
