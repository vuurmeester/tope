#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "math.h"
#include "util.h"

#define EPS 1.0e-6



int linprog(int m, int n, double const* A, double const* b, double const* c,
            double* x)
{
  int stride = 2 * m + n;
  double* mat = alloca(stride * stride * sizeof(double));
  double* y = alloca(m * sizeof(double));       // dual
  double* s = alloca(m * sizeof(double));       // slack
  double* dz = alloca(stride * sizeof(double)); // step

  vec_set(m, y, 1.0);
  vec_set(m, s, 1.0);

  int niter = 0;
  while (true) {
    ++niter;

    // A x - b + s:
    mat_vecmul(m, n, A, x, dz);
    vec_sub(m, dz, b);
    vec_add(m, dz, s);

    // A^T y - c:
    mat_matmul(1, m, n, y, A, dz + m);
    vec_sub(n, dz + m, c);

    // Desired duality gap:
    double nu = 0.01 * vec_dot(m, s, y) / (double)m;
    if (nu < EPS) {
      break;
    }

    // s o y - nu:
    for (int i = 0; i < m; ++i) {
      dz[m + n + i] = s[i] * y[i] - nu;
    }

    // Derivative:
    memset(mat, 0, stride * stride * sizeof(double));
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < n; ++j) {
        mat[i * stride + j] = A[i * n + j];
        mat[(m + j) * stride + i + n] = A[i * n + j];
      }
      mat[i * stride + n + m + i] = 1.0;
      mat[(m + n + i) * stride + n + i] = s[i];
      mat[(m + n + i) * stride + n + m + i] = y[i];
    }

    // Solve mat dz = -err (Newton-Raphson):
    gauss(stride, mat, dz);
    vec_neg(stride, dz);

    // Adjust stepsize to remain feasible:
    double step = 1.0;
    for (int i = 0; i < m; ++i) {
      if (y[i] + step * dz[n + i] < 0.01 * y[i]) {
        step = -0.99 * y[i] / dz[n + i];
      }
      if (s[i] + step * dz[n + m + i] < 0.01 * s[i]) {
        step = -0.99 * s[i] / dz[n + m + i];
      }
    }

    // Perform step:
    vec_adds(n, x, dz, step);
    vec_adds(m, y, dz + n, step);
    vec_adds(m, s, dz + n + m, step);
  }

  printf("%d\n", niter);

  return 0;
}
