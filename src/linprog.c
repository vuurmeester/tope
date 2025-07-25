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

typedef struct {
  int m;
  int n;
  double const* A;
  double const* d;
} Data;



void apply_matrix(int stride, double const* z, double* w, void const* pdata)
{
  Data const* data = pdata;
  int m = data->m;
  int n = data->n;
  double const* A = data->A;
  double const* d = data->d;

  //   0    A    I
  //  A^T   0    0
  //   I    0  S^-1 Y
  //
  // z = [y; x; s]
  // w = [A x + s; A^T y; I y + S^-1 Y s]

  // A x + D s:
  mat_vecmul(m, n, A, z + m, w);
  for (int i = 0; i < m; ++i) {
    w[i] += z[m + n + i] * d[i];
  }

  // A^T y:
  mat_matmul(1, m, n, z, A, w + m);

  // D y + s:
  for (int i = 0; i < m; ++i) {
    w[m + n + i] = z[i] * d[i] + z[m + n + i];
  }
}



int linprog(int m, int n, double const* A, double const* b, double const* c,
            double* x)
{
  int stride = 2 * m + n;
  double* y = alloca(m * sizeof(double));        // dual
  double* s = alloca(m * sizeof(double));        // slack
  double* err = alloca(stride * sizeof(double)); // error
  double* dz = alloca(stride * sizeof(double));  // step
  double* d = alloca(m * sizeof(double));        // precon diagonal

  vec_set(m, y, 1.0);
  vec_set(m, s, 1.0);

  int niter = 0;
  while (true) {
    ++niter;

    // A x - b + s:
    mat_vecmul(m, n, A, x, err);
    vec_sub(m, err, b);
    vec_add(m, err, s);

    // A^T y - c:
    mat_matmul(1, m, n, y, A, err + m);
    vec_sub(n, err + m, c);

    // Desired duality gap (one-tenth of current average duality gap):
    double nu = 0.1 * vec_dot(m, s, y) / (double)m;
    if (nu < EPS) {
      break;
    }

    // S^-1 (s o y) - nu S^-1 e:
    for (int i = 0; i < m; ++i) {
      err[m + n + i] = y[i] - nu / s[i];
    }

    // Preconditioner (creates unit matrix on lower right m x m submatrix):
    for (int i = 0; i < m; ++i) {
      d[i] = sqrt(s[i] / y[i]);
    }

    // Apply preconditioner to RHS:
    for (int i = 0; i < m; ++i) {
      err[m + n + i] *= d[i];
    }

    // Solve mat dz = -err (Newton-Raphson):
    Data data = {.m = m, .n = n, .A = A, .d = d};
    memset(dz, 0, stride * sizeof(double));
    cr(stride, apply_matrix, err, dz, 1e-24, &data);

    // Apply preconditioner:
    for (int i = 0; i < m; ++i) {
      dz[m + n + i] *= d[i];
    }

    vec_neg(stride, dz);

    // Adjust stepsize to remain feasible:
    double step = 1.0;
    for (int i = 0; i < m; ++i) {
      if (y[i] + step * dz[i] < 0.1 * y[i]) {
        step = -0.9 * y[i] / dz[i];
      }
      if (s[i] + step * dz[m + n + i] < 0.1 * s[i]) {
        step = -0.9 * s[i] / dz[m + n + i];
      }
    }

    // Perform step:
    vec_adds(m, y, dz, step);
    vec_adds(n, x, dz + m, step);
    vec_adds(m, s, dz + m + n, step);
  }

#ifndef NDEBUG
  printf("%d\n", niter);
#endif

  return 0;
}
