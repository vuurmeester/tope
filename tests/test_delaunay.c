#include <stdio.h>
#include <stdlib.h>

#include "../src/math.h"
#include "../src/util.h"

#include <tope.h>

#define DIM 2
#define NPOINTS 5



int main(void)
{
  int m;
  int indices[DIM + 1];
  double points[NPOINTS * DIM];
  double weights[DIM + 1];
  double xi[DIM];
  Tope* tope;

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

  /* Create tope object: */
  tope = tope_delaunay(m, 2, points);

  if (tope_getnumfacets(tope) != 4) {
    return -1;
  }

  /* Interpolation: */
  xi[0] = 1;
  xi[1] = 4;
  tope_interpolate(tope, xi, indices, weights);

  /* Sum of weights should be 1: */
  if (fabs(vec_sum(DIM + 1, weights) - 1.0) > 1.0e-6) {
    return -1;
  }

  /* Print the tope: */
  tope_print(tope);

  /* Delete tope object: */
  tope_delete(tope);

  return 0;
}
