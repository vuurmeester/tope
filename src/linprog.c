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



int linprog(
  int m,
  int n,
  double const* A,
  double const* b,
  double const* c,
  double* x
) {
  int stride = m + n + m;
  double* y = alloca(m * sizeof(double));        /* dual */
  double* s = alloca(m * sizeof(double));        /* slack */
  double* err = alloca(stride * sizeof(double)); /* error */
  double* mat = alloca(stride * stride * sizeof(double));

  vec_set(m, y, 1.0);
  vec_set(m, s, 1.0);

  int niter = 0;
  while (true) {
    ++niter;

    /* Desired duality gap (one-tenth of current average duality gap): */
    double nu = 0.1 * vec_dot(m, s, y) / (double)m;
    if (nu < EPS) {
      break;
    }

    /*    0   A     1    */
    /*   A^T  0     0    */
    /*    1   0   S^-1 Y */
    memset(mat, 0x00, stride * stride * sizeof(double));
    for (int i = 0; i < m; ++i) {
      mat[i * stride + m + n + i] = 1.0;
      mat[(m + n + i) * stride + i] = 1.0;
      mat[(m + n + i) * stride + m + n + i] = y[i] / s[i];
      for (int j = 0; j < n; ++j) {
        mat[i       * stride + m + j] = A[i * n + j];
        mat[(m + j) * stride + i    ] = A[i * n + j];
      }
    }

    /* A x - b + s: */
    mat_vecmul(m, n, A, x, err);
    vec_sub(m, err, b);
    vec_add(m, err, s);

    /* A^T y - c: */
    mat_matmul(1, m, n, y, A, err + m);
    vec_sub(n, err + m, c);

    /* S^-1 (s o y) - nu S^-1 e: */
    for (int i = 0; i < m; ++i) {
      err[m + n + i] = y[i] - nu / s[i];
    }

    /* Solve M dz = -err (Newton-Raphson): */
    gauss(stride, mat, err);
    vec_neg(stride, err);

    /* Adjust stepsize to remain feasible (y, s > 0): */
    double step = 1.0;
    for (int i = 0; i < m; ++i) {
      if (y[i] + step * err[i] < 0.0) {
        step = -y[i] / err[i];
      }
      if (s[i] + step * err[m + n + i] < 0.0) {
        step = -s[i] / err[m + n + i];
      }
    }
    step *= 0.98;

    /* Perform step: */
    vec_adds(m, y, err, step);
    vec_adds(n, x, err + m, step);
    vec_adds(m, s, err + m + n, step);

    if (step * vec_norm(stride, err) <= EPS) {
      break;
    }
  }

#ifndef NDEBUG
  printf("%d\n", niter);
#endif

  return 0;
}
