#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "math.h"
#include "util.h"



void vec_add(int n, double* x, double const* y)
{
  int i;
  for (i = 0; i < n; ++i) {
    x[i] += y[i];
  }
}



void vec_sub(int n, double* x, double const* y)
{
  int i;
  for (i = 0; i < n; ++i) {
    x[i] -= y[i];
  }
}



void vec_scale(int n, double* x, double scale)
{
  int i;
  for (i = 0; i < n; ++i) {
    x[i] *= scale;
  }
}



void vec_set(int n, double* x, double scalar)
{
  int i;
  for (i = 0; i < n; ++i) {
    x[i] = scalar;
  }
}



void vec_reset(int n, double* x)
{
  vec_set(n, x, 0.0);
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
  int i;
  double result;

  result = 0.0;
  for (i = 0; i < n; ++i) {
    result += x[i];
  }

  return result;
}



double vec_nrmsq(int n, double const* x)
{
  double result;
  int i;

  result = 0.0;
  for (i = 0; i < n; ++i) {
    result += x[i] * x[i];
  }

  return result;
}



double vec_norm(int n, double const* x)
{
  return sqrt(vec_nrmsq(n, x));
}



void vec_adds(int n, double* x, double const* y, double scale)
{
  int i;
  for (i = 0; i < n; ++i) {
    x[i] += scale * y[i];
  }
}



void vec_neg(int n, double* x)
{
  int i;
  for (i = 0; i < n; ++i) {
    x[i] = -x[i];
  }
}



double vec_dot(int n, double const* x, double const* y)
{
  int i;
  double result;

  result = 0.0;
  for (i = 0; i < n; ++i) {
    result += x[i] * y[i];
  }

  return result;
}



int vec_minindex(int n, double const* x)
{
  int result = 0;
  int i;

  for (i = 1; i < n; ++i) {
    if (x[i] < x[result]) {
      result = i;
    }
  }

  return result;
}



int vec_maxindex(int n, double const* x)
{
  int result = 0;
  int i;

  for (i = 1; i < n; ++i) {
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



void mat_sprint(int m, int n, double const* mat, char* str)
{
  int width;
  int i;
  int j;
  int candidate;
  char format[256];
  char buffer[32];
  double maxel;
  double number;

  /* Determine the largest entry: */
  maxel = 0.0;
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) {
      if (fabs(mat[i * n + j]) > maxel) {
        maxel = fabs(mat[i * n + j]);
      }
    }
  }

  /* Determine the width of the widest entry: */
  width = 0;
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) {
      number = mat[i * n + j];
      if (fabs(number) == 0.0) {
        number = 0.0;
      }

      candidate = sprintf(buffer, "%0.6g", number);
      if (candidate > width) {
        width = candidate;
      }
    }
  }

  width += 2; /* text separation */

  sprintf(format, "%%%d.6g", width);
  /* Now actually print the data: */
  for (i = 0; i < m; ++i) {
    for (j = 0; j < n; ++j) {
      number = mat[i * n + j];
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
  int size;
  char* buffer;

  size = (1 + m + n + m * n) * 32;
  buffer = malloc(size);
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
  int i;

  /* Clear the matrix: */
  matrix_reset(m, m, mat);

  /* Fill diagonal with ones: */
  for (i = 0; i < m; ++i) {
    mat[i * m + i] = 1.0;
  }
}



void boundingbox(int npoints, int ndims, double const* vertices,
                 int* minindices, int* maxindices, double* minima,
                 double* maxima)
{
  int i;
  int j;

  /* Initialize minima and maxima: */
  for (j = 0; j < ndims; ++j) {
    minima[j] = HUGE_VAL;
    maxima[j] = -HUGE_VAL;
  }

  /* Loop over vertices: */
  for (i = 0; i < npoints; ++i) {
    /* Loop over dimensions: */
    for (j = 0; j < ndims; ++j) {
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
  int i;
  int j;

  assert(npoints <= ndims + 1);

  /* Accumulate centroid: */
  vec_reset(ndims, centroid);
  for (i = 0; i < npoints; ++i) {
    vec_add(ndims, centroid, &points[i * ndims]);
  }
  vec_scale(ndims, centroid, 1.0 / (double)npoints);

  /* Subtract last point from all points: */
  for (i = 0; i < npoints - 1; ++i) {
    vec_sub(ndims, &points[i * ndims], &points[(npoints - 1) * ndims]);
  }
  memset(&points[(npoints - 1) * ndims], 0, ndims * sizeof(double));

  *volume = 1.0;
  for (i = 0; i < npoints - 1; ++i) {
    /* Get pivot: */
    double maxnormsq = 0.0;
    int pivot = -1;
    for (j = i; j < npoints - 1; ++j) {
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
    for (j = i + 1; j < npoints - 1; ++j) {
      double fac = vec_dot(ndims, &points[i * ndims], &points[j * ndims]);
      vec_adds(ndims, &points[j * ndims], &points[i * ndims], -fac);
    }
  }
}
