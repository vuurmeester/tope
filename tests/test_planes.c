#include "../src/math.h"

#include <polytoop.h>

#define DIM 3
#define NPLANES 5



int main()
{
  double normals[(NPLANES + 1) * DIM];
  double distances[NPLANES + 1];

  /* Define a pyramid. */
  int m = 0;

  /* Bottom: */
  normals[m * DIM + 0] = 0.0;
  normals[m * DIM + 1] = 0.0;
  normals[m * DIM + 2] = -1.0;
  distances[m] = 0.0;
  ++m;

  for (; m < NPLANES; ++m) {
    double angle = (double)(m - 1) / (double)(NPLANES - 1) * M_2PI;
    normals[m * DIM + 0] = M_SQRT1_2 * cos(angle);
    normals[m * DIM + 1] = M_SQRT1_2 * sin(angle);
    normals[m * DIM + 2] = M_SQRT1_2;
    distances[m] = M_SQRT1_2;
  }

  /* Add a redundant plane: */
  normals[m * DIM + 0] = 0.0;
  normals[m * DIM + 1] = 0.0;
  normals[m * DIM + 2] = 1.0;
  distances[m] = 2.0;
  ++m;

  /* Compute polytoop from pyramid planes: */
  Polytoop* polytoop = polytoop_fromplanes(m, DIM, normals, distances);

  if (polytoop) {
    /* Print the polytoop: */
    polytoop_print(polytoop);

    /* Delete polytoop object: */
    polytoop_delete(polytoop);
  }

  return 0;
}
