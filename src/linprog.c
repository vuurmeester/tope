#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "math.h"
#include "util.h"

#define EPS 1.0e-9



int linprog(int m, int n, double const* A, double const* b, double const* c,
            double* x)
{
  int ret = 0;
  double* dirs = alloca(m * n * sizeof(double));
  double* d = alloca(m * sizeof(double));
  double normc = vec_norm(n, c);
  double* movedir = alloca(n * sizeof(double));

  while (true) {
    // d = A * x - b:
    mat_vecmul(m, n, A, x, d);
    vec_sub(m, d, b);
    double normd = vec_norm(m, d);
    if (d[vec_maxindex(m, d)] > EPS * normd) {
      // infeasible
      return -2;
    }

    // Determine active set:
    int nactive = 0;
    for (int i = 0; i < m; ++i) {
      if (d[i] > -EPS * normd) {
        memcpy(dirs + nactive * n, A + i * n, n * sizeof(double));
        ++nactive;
      }
    }

    // Orthogonalize move direction:
    memcpy(movedir, c, n * sizeof(double));
    vec_scale(n, movedir, 1.0 / normc);
    int i = 0;
    for (; i < nactive; ++i) {
      // Find pivot row:
      double maxip = -HUGE_VAL;
      int pivot = -1;
      for (int j = i; j < nactive; ++j) {
        double ip = vec_dot(n, dirs + j * n, movedir);
        if (ip > maxip) {
          maxip = ip;
          pivot = j;
        }
      }

      if (maxip < EPS) {
        break;
      }

      // Perform pivot:
      if (pivot > i) {
        memswp(dirs + i * n, dirs + pivot * n, n * sizeof(double));
      }

      double normsq = vec_nrmsq(n, dirs + i * n);
      // Orthogonalize:
      vec_adds(n, movedir, dirs + i * n,
               -vec_dot(n, movedir, dirs + i * n) / normsq);
      for (int j = i + 1; j < nactive; ++j) {
        vec_adds(n, dirs + j * n, dirs + i * n,
                 -vec_dot(n, dirs + j * n, dirs + i * n) / normsq);
        assert(fabs(vec_dot(n, dirs + i * n, dirs + j * n)) < EPS);
      }
    }

    double movedir_norm = vec_norm(n, movedir);
    if (movedir_norm < EPS) {
      // Stuck at optimum.
      break;
    }
    vec_scale(n, movedir, 1.0 / movedir_norm);

    // Find most limiting plane:
    int minindex = -1;
    double minmovedist = HUGE_VAL;
    for (int i = 0; i < m; ++i) {
      if (d[i] > -EPS * normd) {
        // active dir
        continue;
      }
      double c = vec_dot(n, movedir, A + i * n);
      if (c < EPS * vec_norm(n, A + i * n)) {
        // non interfering
        continue;
      }
      double movedist = -d[i] / c;
      if (movedist < minmovedist) {
        minmovedist = movedist;
        minindex = i;
      }
    }

    // Unbounded?
    if (minindex == -1) {
      return -1;
    }

    if (minmovedist <= EPS * normd) {
      break;
    }

    // Move to most limiting plane:
    vec_adds(n, x, movedir, minmovedist);
  }

  return ret;
}
