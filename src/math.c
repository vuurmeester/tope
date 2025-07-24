#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "math.h"
#include "util.h"



void vec_add(int n, double* x, double const* y)
{
  while (n--) {
    x[n] += y[n];
  }
}



void vec_sub(int n, double* x, double const* y)
{
  while (n--) {
    x[n] -= y[n];
  }
}



void vec_scale(int n, double* x, double scale)
{
  while (n--) {
    x[n] *= scale;
  }
}



void vec_set(int n, double* x, double scalar)
{
  while (n--) {
    x[n] = scalar;
  }
}



void vec_reset(int n, double* x)
{
  memset(x, 0, n * sizeof(double));
}



void vec_normalize(int n, double* x)
{
  double nrmsq = vec_nrmsq(n, x);
  if (nrmsq > 0.0) {
    vec_scale(n, x, 1.0 / sqrt(nrmsq));
  }
}



double vec_sum(int n, double const* x)
{
  double result = 0.0;
  while (n--) {
    result += x[n];
  }
  return result;
}



double vec_nrmsq(int n, double const* x)
{
  double result = 0.0;
  while (n--) {
    result += x[n] * x[n];
  }
  return result;
}



double vec_norm(int n, double const* x)
{
  return sqrt(vec_nrmsq(n, x));
}



void vec_adds(int n, double* x, double const* y, double scale)
{
  while (n--) {
    x[n] += scale * y[n];
  }
}



void vec_neg(int n, double* x)
{
  while (n--) {
    x[n] = -x[n];
  }
}



double vec_dot(int n, double const* x, double const* y)
{
  double result = 0.0;
  while (n--) {
    result += x[n] * y[n];
  }
  return result;
}



int vec_minindex(int n, double const* x)
{
  int result = 0;
  for (int i = 1; i < n; ++i) {
    if (x[i] < x[result]) {
      result = i;
    }
  }
  return result;
}



int vec_maxindex(int n, double const* x)
{
  int result = 0;
  for (int i = 1; i < n; ++i) {
    if (x[i] > x[result]) {
      result = i;
    }
  }
  return result;
}



void vector_sprint(int n, double const* x, char* str)
{
  mat_sprint(1, n, x, str);
}



void vector_fprint(int n, double const* x, FILE* outstream)
{
  char buffer[2048];
  vector_sprint(n, x, buffer);
  fprintf(outstream, "%s", buffer);
}



void vec_print(int n, double const* x)
{
  vector_fprint(n, x, stdout);
}



void mat_vecmul(int m, int n, double const* mat, double const* x, double* y)
{
  for (int i = 0; i < m; ++i) {
    y[i] = vec_dot(n, mat + i * n, x);
  }
}



void mat_matmul(int m, int n, int o, double const* mat1, double const* mat2, double* result)
{
  memset(result, 0, m * o * sizeof(double));
  for (int i = 0; i < m; ++i) {
    for (int k = 0; k < n; ++k) {
      for (int j = 0; j < o; ++j) {
        result[i * n + j] += mat1[i * n + k] * mat2[k * o + j];
      }
    }
  }
}



void mat_sprint(int m, int n, double const* mat, char* str)
{
  char format[256];
  char buffer[32];

  /* Determine the largest entry: */
  double maxel = 0.0;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      if (fabs(mat[i * n + j]) > maxel) {
        maxel = fabs(mat[i * n + j]);
      }
    }
  }

  /* Determine the width of the widest entry: */
  int width = 0;
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      double number = mat[i * n + j];
      if (fabs(number) == 0.0) {
        number = 0.0;
      }

      int candidate = sprintf(buffer, "%0.6g", number);
      if (candidate > width) {
        width = candidate;
      }
    }
  }

  width += 2; /* text separation */

  sprintf(format, "%%%d.6g", width);
  /* Now actually print the data: */
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < n; ++j) {
      double number = mat[i * n + j];
      if (fabs(number) == 0.0) {
        number = 0.0;
      }
      str += sprintf(str, format, number);
    }
    str += sprintf(str, "\n");
  }
}



void mat_print(int m, int n, double const* mat)
{
  mat_fprint(m, n, mat, stdout);
}



void mat_fprint(int m, int n, double const* mat, FILE* outstream)
{
  int size = (1 + m + n + m * n) * 32;
  char* buffer = malloc(size);
  memset(buffer, 0, size);
  mat_sprint(m, n, mat, buffer);
  fprintf(outstream, "%s", buffer);
  free(buffer);
}



void matrix_reset(int m, int n, double* mat)
{
  vec_reset(m * n, mat);
}



void mat_unit(int m, double* mat)
{
  /* Clear the matrix: */
  matrix_reset(m, m, mat);

  /* Fill diagonal with ones: */
  for (int i = 0; i < m; ++i) {
    mat[i * m + i] = 1.0;
  }
}



void boundingbox(int npoints, int ndims, double const* vertices,
                 int* minindices, int* maxindices, double* minima,
                 double* maxima)
{
  /* Initialize minima and maxima: */
  for (int j = 0; j < ndims; ++j) {
    minima[j] = HUGE_VAL;
    maxima[j] = -HUGE_VAL;
  }

  /* Loop over vertices: */
  for (int i = 0; i < npoints; ++i) {
    /* Loop over dimensions: */
    for (int j = 0; j < ndims; ++j) {
      if (vertices[i * ndims + j] < minima[j]) {
        /* Update minimum of dimension j: */
        minindices[j] = i;
        minima[j] = vertices[i * ndims + j];
      }

      if (vertices[i * ndims + j] > maxima[j]) {
        /* Update maximum of dimension j: */
        maxindices[j] = i;
        maxima[j] = vertices[i * ndims + j];
      }
    }
  }
}



void analysesimplex(int npoints, int ndims, double* points, double* volume,
                    double* centroid)
{
  assert(npoints <= ndims + 1);

  /* Accumulate centroid: */
  vec_reset(ndims, centroid);
  for (int i = 0; i < npoints; ++i) {
    vec_add(ndims, centroid, &points[i * ndims]);
  }
  vec_scale(ndims, centroid, 1.0 / (double)npoints);

  /* Subtract last point from all points: */
  for (int i = 0; i < npoints - 1; ++i) {
    vec_sub(ndims, &points[i * ndims], &points[(npoints - 1) * ndims]);
  }
  memset(&points[(npoints - 1) * ndims], 0, ndims * sizeof(double));

  *volume = 1.0;
  for (int i = 0; i < npoints - 1; ++i) {
    /* Get pivot: */
    double maxnormsq = 0.0;
    int pivot = -1;
    for (int j = i; j < npoints - 1; ++j) {
      double normsq = vec_nrmsq(ndims, &points[j * ndims]);
      if (normsq > maxnormsq) {
        maxnormsq = normsq;
        pivot = j;
      }
    }

    /* Update volume: */
    double nrm = sqrt(maxnormsq);
    *volume *= nrm / (double)(i + 1);
    assert(nrm > 0.0);

    /* Perform pivot: */
    if (pivot > i) {
      memswp(&points[i * ndims], &points[pivot * ndims],
             ndims * sizeof(double));
    }

    /* Normalize: */
    vec_scale(ndims, &points[i * ndims], 1.0 / nrm);

    /* Orthogonalize (subtract projection of row i from those below): */
    for (int j = i + 1; j < npoints - 1; ++j) {
      double fac = vec_dot(ndims, &points[i * ndims], &points[j * ndims]);
      vec_adds(ndims, &points[j * ndims], &points[i * ndims], -fac);
    }
  }
}



double gauss(int n, double* A, double* b)
{
  double det = 1.0;
  for (int i = 0; i < n; ++i) {
    // Find pivot element (largest magnitude in column i):
    int pivot = i;
    for (int j = i + 1; j < n; ++j) {
      if (fabs(A[j * n + i]) > fabs(A[pivot * n + i])) {
        pivot = j;
      }
    }

    // Perform pivot:
    if (pivot > i) {
      memswp(A + i * n, A + pivot * n, n * sizeof(double));
      memswp(b + i, b + pivot, sizeof(double));
      det *= -1.0;
    }

    // Guard singular:
    if (A[i * n + i] == 0.0) {
      return 0.0;
    }

    // Update determinant:
    det *= A[i * n + i];

    // Update U:
    for (int j = i + 1; j < n; ++j) {
      double mult = A[j * n + i] / A[i * n + i];
      vec_adds(n - i - 1, A + j * n + i + 1, A + i * n + i + 1, -mult);
      b[j] -= mult * b[i];
    }
  }

  // Backsubstitution:
  for (int i = n - 1; i >= 0; --i) {
    b[i] = (b[i] - vec_dot(n - i - 1, b + i + 1, A + i * n + i + 1)) / A[i * n + i];
  }

  return det;
}
