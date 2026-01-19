#include <stdlib.h>
#include <stdio.h>

#include "../src/util.h"
#include "../src/math.h"

#include <tope.h>

#define DIM 8
#define NPOINTS (2 << DIM)



int main()
{
  double points[NPOINTS * DIM];
  Tope* tope;

  /* Define a cube. */
  for (int i = 0; i < NPOINTS; ++i) {
    for (int j = 0; j < DIM; ++j) {
      points[i * DIM + j] = ((i >> (DIM - 1 - j)) & 1) - 0.5;
    }
  }

  /* Initialize tope with all but the last: */
  tope = tope_frompoints(NPOINTS, DIM, points);
  tope_merge(tope);

  /* Print the tope: */
  tope_print(tope);

  if (tope_getnumfacets(tope) != 2 * DIM) {
    return -1;
  }

  /* Delete tope object: */
  tope_delete(tope);

  return 0;
}
