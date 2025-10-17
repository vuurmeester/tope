#include "../src/math.h"

#include <polytoop.h>

#define DIM 3
#define NPLANES 5



int main(void)
{
  Polytoop* polytoop;
  double normals[(NPLANES + 1) * DIM];
  double distances[NPLANES + 1];

  /* Define a pyramid. */
  int m = 0;
  double xi[DIM];

  /* Bottom: */
  normals[m * DIM + 0] = 0.0;
  normals[m * DIM + 1] = 0.0;
  normals[m * DIM + 2] = -1.0;
  distances[m] = 0.0;
  ++m;

  double height = 4.0;  // nice 3 4 5 triangles
  double offset = 3.0;
  double steepness = atan2(height, offset);
  for (; m < NPLANES; ++m) {
    double angle = (double)(m - 1) / (double)(NPLANES - 1) * M_2PI;
    normals[m * DIM + 0] = sin(steepness) * cos(angle);
    normals[m * DIM + 1] = sin(steepness) * sin(angle);
    normals[m * DIM + 2] = cos(steepness);
    distances[m] = offset * sin(steepness);
  }

  /* Add a redundant plane: */
  normals[m * DIM + 0] = 0.0;
  normals[m * DIM + 1] = 0.0;
  normals[m * DIM + 2] = 1.0;
  distances[m] = 5.0;
  ++m;

  /* Compute polytoop from pyramid planes: */
  xi[0] = 0.0;
  xi[1] = 0.0;
  xi[2] = 0.25;
  polytoop = polytoop_fromplanes(m, DIM, normals, distances, xi);

  if (polytoop) {
    /* Print the polytoop: */
    polytoop_print(polytoop);

    /* Delete polytoop object: */
    polytoop_delete(polytoop);
  }

  return 0;
}
