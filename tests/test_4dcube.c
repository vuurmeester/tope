#include <stdlib.h>
#include <stdio.h>

#include "../src/util.h"
#include "../src/math.h"

#include <tope.h>

#define DIM 4
#define NPOINTS (2 << DIM)



int main()
{
  int m;
  double points[NPOINTS * DIM];
  Tope* tope;

  /* Define a cube. */
  for (m = 0; m < NPOINTS; ++m) {
    for (int i = 0; i < DIM; ++i) {
      points[m * DIM + i] = ((m >> (DIM - 1 - i)) & 1) - 0.5;
    }
  }

  /* Initialize tope with all but the last: */
  tope = tope_frompoints(NPOINTS, DIM, points, true);

  /* Print the tope: */
  tope_print(tope);

  if (tope_getnumfacets(tope) != 2 * DIM) {
    return -1;
  }

  /* Delete tope object: */
  tope_delete(tope);

  return 0;
}
