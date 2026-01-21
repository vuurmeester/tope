#include <stdlib.h>
#include <stdio.h>

#include "../src/util.h"
#include "../src/math.h"

#include <tope.h>

#define DIM 3
#define NPOINTS 20



int main()
{
  double points[NPOINTS * DIM];

  /* Define a cube. */
  int m = 0;

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

  /* Initialize tope with all but the last: */
  Tope* tope = tope_frompoints(m - 1, 3, points);
  tope_merge(tope);

  /* Print the tope: */
  tope_print(tope);

  if (tope_getnumfacets(tope) != 6) {
    return -1;
  }

  /* Add the last point: */
  tope_addvertex(tope, &points[(m - 1) * DIM]);

  /* Print the tope: */
  tope_print(tope);

  if (tope_getnumfacets(tope) != 9) {
    return -1;
  }

  /* Delete tope object: */
  tope_delete(tope);

  return 0;
}
