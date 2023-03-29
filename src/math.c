#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "math.h"
#include "util.h"



void vector_add(int n, double* x, double const* y)
{
  int i;
  for (i = 0; i < n; ++i) {
    x[i] += y[i];
  }
}



void vector_subtract(int n, double* x, double const* y)
{
  int i;
  for (i = 0; i < n; ++i) {
    x[i] -= y[i];
  }
}



void vector_scale(int n, double* x, double scale)
{
  int i;
  for (i = 0; i < n; ++i) {
    x[i] *= scale;
  }
}



void vector_set(int n, double* x, double scalar)
{
  int i;
  for (i = 0; i < n; ++i) {
    x[i] = scalar;
  }
}



void vector_reset(int n, double* x)
{
  vector_set(n, x, 0.0);
}



void vector_normalize(int n, double* x)
{
  double nrmsq = vector_normsq(n, x);
  if (nrmsq > 0.0) {
    vector_scale(n, x, 1.0 / sqrt(nrmsq));
  }
}



double vector_sum(int n, double const* x)
{
  int i;
  double result;

  result = 0.0;
  for (i = 0; i < n; ++i) {
    result += x[i];
  }

  return result;
}



double vector_normsq(int n, double const* x)
{
  double result;
  int i;

  result = 0.0;
  for (i = 0; i < n; ++i) {
    result += x[i] * x[i];
  }

  return result;
}



void vector_add_scaled(int n, double* x, double const* y, double scale)
{
  int i;
  for (i = 0; i < n; ++i) {
    x[i] += scale * y[i];
  }
}



void vector_negate(int n, double* x)
{
  int i;
  for (i = 0; i < n; ++i) {
    x[i] = -x[i];
  }
}



double vector_ip(int n, double const* x, double const* y)
{
  int i;
  double result;

  result = 0.0;
  for (i = 0; i < n; ++i) {
    result += x[i] * y[i];
  }

  return result;
}



int vector_minindex(int n, double const* x)
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



void vector_sprint(int n, double const* x, char* str)
{
  matrix_sprint(1, n, x, str);
}



void vector_fprint(int n, double const* x, FILE* outstream)
{
  char buffer[2048];
 vector_sprint(n, x, buffer);
  fprintf(outstream, "%s", buffer);
}



void vector_print(int n, double const* x)
{
  vector_fprint(n, x, stdout);
}



void matrix_sprint(int m, int n, double const* mat, char* str)
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

  width += 2;  /* text separation */

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



void matrix_print(int m, int n, double const* mat)
{
  matrix_fprint(m, n, mat, stdout);
}



void matrix_fprint(int m, int n, double const* mat, FILE* outstream)
{
  int size;
  char* buffer;

  size = (1 + m + n + m * n) * 32;
  buffer = malloc(size);
  memset(buffer, 0, size);
  matrix_sprint(m, n, mat, buffer);
  fprintf(outstream, "%s", buffer);
  free(buffer);
}



void matrix_reset(int m, int n, double* mat)
{
  vector_reset(m * n, mat);
}



void matrix_unit(int m, double* mat)
{
  int i;

  /* Clear the matrix: */
  matrix_reset(m, m, mat);

  /* Fill diagonal with ones: */
  for (i = 0; i < m; ++i) {
    mat[i * m + i] = 1.0;
  }
}



void boundingbox(int npoints,
  int ndims,
  double const* vertices,
  int* minindices,
  int* maxindices,
  double* minima,
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



void analysesimplex(int npoints,
  int ndims,
  double const* points,
  double* volume,
  double* centroid,
  double* span)
{
  int i;
  double* lqdcmp;
  int* p;

  /* Accumulate centroid: */
  vector_reset(ndims, centroid);
  for (i = 0; i < npoints; ++i) {
    vector_add(ndims, centroid, &points[i * ndims]);
  }
  vector_scale(ndims, centroid, 1.0 / (double)npoints);

  /* Determine directions: */
  for (i = 0; i < npoints - 1; ++i) {
    memcpy(&span[i * ndims], &points[(i + 1) * ndims], ndims * sizeof(double));
    vector_subtract(ndims, &span[i * ndims], &points[0 * ndims]);
  }

  /* LQ decomposition: */
  lqdcmp = alloca((npoints - 1) * ndims * sizeof(double));
  memcpy(lqdcmp, span, (npoints - 1) * ndims * sizeof(double));
  p = alloca((npoints - 1) * sizeof(int));
  lqdc(npoints - 1, ndims, lqdcmp, p);

  /* Volume: */
  *volume = fabs(lqdcmp[0 * ndims + 0]);
  for (i = 1; i < npoints - 1; ++i) {
    *volume *= fabs(lqdcmp[i * ndims + i]) / (double)(i + 1);
  }

  /* Directions: */
  lqformq(npoints - 1, ndims, lqdcmp, span);
}



double lqdc(int m, int n, double* matrix, int* p)
{
  int i;
  int j;
  int pivot;
  int mindim;
  double det;
  double alfa;
  double alfa1;
  double alfamax;
  double beta;
  double ip;

  /* Initialize determinant and permutation vector p: */
  det = 1.0;
  if (p) {
    for (i = 0; i < m; ++i) {
      p[i] = i;
    }
  }

  mindim = m < n ? m : n;
  alfamax = 0.0;
  for (i = 0; i < mindim; ++i) {
    /* Squared norm of row i: */
    alfa = vector_normsq(n - i, &matrix[i * n + i]);

    if (p) {
      /* Get largest row norm on top: */
      pivot = i;

      /* Find pivot row: */
      for (j = i + 1; j < m; ++j) {
        alfa1 = vector_normsq(n - i, &matrix[j * n + i]);
        if (alfa1 > alfa) {
          alfa = alfa1;
          pivot = j;
        }
      }

      /* Pivot: */
      if (pivot > i) {
        /* Flip parity for each pivot: */
        det *= -1.0;

        /* Swap rows i and pivot: */
        memswp(&matrix[pivot * n], &matrix[i * n], n * sizeof(double));

        /* Adjust permutation vector: */
        memswp(&p[i], &p[pivot], sizeof(int));
      }
    }

    if (alfa > alfamax) {
      alfamax = alfa;
    }
    else if (!(alfa > 1.0e-24 * alfamax)) {
      det = 0.0;
      break;
    }

    /* Norm of row i: */
    alfa = sqrt(alfa);

    /* Make sign opposite to x_1: */
    if (matrix[i * n + i] > 0.0) {
      alfa = -alfa;
    }

    beta = matrix[i * n + i] - alfa;

    /* Normalize v_2.. so that v1 = 1: */
    vector_scale(n - i - 1, &matrix[i * n + i + 1], 1.0 / beta);

    beta /= alfa;

    /* Transform lower triangle of submatrix (i:m,i:n): */
    matrix[i * n + i] = alfa;
    for (j = i + 1; j < m; ++j) {
      /* v^T x: */
      ip = matrix[j * n + i] + vector_ip(n - i - 1, &matrix[i * n + i + 1], &matrix[j * n + i + 1]);

      /* x - 2 v^T x / (v^T v) v: */
      ip *= beta;
      matrix[j * n + i] += ip;  /* v1 == 1 */
      vector_add_scaled(n - i - 1, &matrix[j * n + i + 1], &matrix[i * n + i + 1], ip);
    }

    /* Update determinant: */
    det *= -alfa;
  }

  return det;
}



void lqbs(int m,
               int n,
               double const* dcmp,
               int const* p,
               double const* b,
               double* x,
               double tol)
{
  int i;
  int mindim;
  double beta;
  double tolabsA00;
  double ip;

  /* Relative tolerance: */
  tolabsA00 = tol * fabs(dcmp[0]);

  /* Minimum dimension: */
  mindim = m < n ? m : n;

  /* Calculate x = L \ (P * b): */
  if (p) {
    for (i = 0; i < mindim; ++i) {
      x[i] = b[p[i]];
    }
  }
  else {
    memcpy(x, b, mindim * sizeof(double));
  }
  for (i = 0; i < mindim; ++i) {
    if (fabs(dcmp[i * n + i]) <= tolabsA00) {
      vector_reset(n - i, &x[i]);
      m = i;
      break;
    }
    x[i] -= vector_ip(i, &dcmp[i * n], x);
    x[i] /= dcmp[i * n + i];
  }

  /* Replace x <-- Q' * x: */
  for (i = mindim - 1; i >= 0; --i) {
    /* Determine beta: */
    beta = 1.0 + vector_normsq(n - i - 1, &dcmp[i * n + i + 1]);
    beta = 2.0 / beta;  /* 2 / (v^T v) */

    /* Determine inner product: */
    ip = x[i] + vector_ip(n - i - 1, &dcmp[i * n + i + 1], &x[i + 1]);

    /* Transform the vector: */
    ip *= beta;
    x[i] -= ip;
    vector_add_scaled(n - i - 1, &x[i + 1], &dcmp[i * n + i + 1], -ip);
  }
}



void lqformq(int m, int n, double const* dcmp, double* matq)
{
  int i;
  int j;
  double vtv;
  double beta;
  double ip;

  /* Unit matrix: */
  matrix_unit(n, matq);

  /* Apply orthogonal transformation Q_1 * ... * Q_{n - 1} (last to first): */
  for (i = m < n ? m - 1 : n - 1; i >= 0; --i) {
    /* Squared norm of Householder vector v: */
    vtv = 1.0 + vector_normsq(n - i - 1, &dcmp[i * n + i + 1]);

    /* Two times reciprocal of squared norm: */
    beta = 2.0 / vtv;

    /* Transform row i of Q which equals (1, 0, ..., 0): */
    matq[i * n + i] -= beta;  /* v_1 is implicitly 1 */
    vector_add_scaled(n - i - 1, &matq[i * n + i + 1], &dcmp[i * n + i + 1], -beta);

    /* Transform rows j = i + 1 .. n - 1: */
    for (j = i + 1; j < n; ++j) {
      /* Inner product of v with row j of Q: */
      ip = matq[j * n + i] + vector_ip(n - i - 1, &dcmp[i * n + i + 1], &matq[j * n + i + 1]);

      /* Transform row j of Q: */
      ip *= beta;
      matq[j * n + i] -= ip;  /* v_1 is implicitly 1 */
      vector_add_scaled(n - i - 1, &matq[j * n + i + 1], &dcmp[i * n + i + 1], -ip);
    }
  }
}
