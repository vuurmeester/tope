#include <stdlib.h>
#include <stdio.h>

#include "../src/util.h"
#include "../src/math.h"

#include <polytoop.h>

#define DIM 3
#define NPOINTS 20



int main()
{
  int m;
  double points[NPOINTS * DIM];
  Polytoop* polytoop;

  /* Define a cube. */
  m = 0;

  points[m * DIM + 0] = 0.0;
  points[m * DIM + 1] = 0.0;
  points[m * DIM + 2] = 0.0;
  ++m;

  points[m * DIM + 0] = 0.0;
  points[m * DIM + 1] = 0.0;
  points[m * DIM + 2] = 1.0;
  ++m;

  points[m * DIM + 0] = 0.0;
  points[m * DIM + 1] = 1.0;
  points[m * DIM + 2] = 0.0;
  ++m;

  points[m * DIM + 0] = 0.0;
  points[m * DIM + 1] = 1.0;
  points[m * DIM + 2] = 1.0;
  ++m;

  points[m * DIM + 0] = 1.0;
  points[m * DIM + 1] = 0.0;
  points[m * DIM + 2] = 0.0;
  ++m;

  points[m * DIM + 0] = 1.0;
  points[m * DIM + 1] = 0.0;
  points[m * DIM + 2] = 1.0;
  ++m;

  points[m * DIM + 0] = 1.0;
  points[m * DIM + 1] = 1.0;
  points[m * DIM + 2] = 0.0;
  ++m;

  points[m * DIM + 0] = 1.0;
  points[m * DIM + 1] = 1.0;
  points[m * DIM + 2] = 1.0;
  ++m;

  /* Internal point: */
  points[m * DIM + 0] = 0.5;
  points[m * DIM + 1] = 0.5;
  points[m * DIM + 2] = 0.5;
  ++m;

  /* Duplicate point: */
  points[m * DIM + 0] = 1.0;
  points[m * DIM + 1] = 0.0;
  points[m * DIM + 2] = 0.0;
  ++m;

  /* Point on facet: */
  points[m * DIM + 0] = 0.0;
  points[m * DIM + 1] = 0.5;
  points[m * DIM + 2] = 0.5;
  ++m;

  /* Point on edge: */
  points[m * DIM + 0] = 1.0;
  points[m * DIM + 1] = 0.0;
  points[m * DIM + 2] = 0.5;
  ++m;

  /* External point, added later: */
  points[m * DIM + 0] = 0.5;
  points[m * DIM + 1] = 0.5;
  points[m * DIM + 2] = 2.0;
  ++m;

  /* Create polytoop object: */
  polytoop = polytoop_new();
  polytoop_setmerge(polytoop, 1);

  /* Initialize polytoop with all but the last: */
  polytoop_frompoints(polytoop, m - 1, 3, points);
  if (polytoop_getnumfacets(polytoop) != 6) {
    return -1;
  }

  /* Add the last point: */
  polytoop_addvertex(polytoop, &points[(m - 1) * DIM]);
  if (polytoop_getnumfacets(polytoop) != 9) {
    return -1;
  }

  /* Print the polytoop: */
  polytoop_print(polytoop);

  /* Delete polytoop object: */
  polytoop_delete(polytoop);

  return 0;
}
