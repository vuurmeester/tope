#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "math.h"
#include "util.h"

#define EPS 1.0e-12



int linprog_rn(int m, int n, double const* mata, double const* b,
               double const* c, double* x)
{
  int* basis = (int*)alloca(m * sizeof(int));
  for (int i = 0; i < m; ++i) {
    basis[i] = -1;
  }

  /* Everything non-basic: */
  int nart = m;

  /* Extract basis: */
  for (int j = 0; j < n; ++j) {
    /* Search for column with exactly one positive number and otherwise zero: */
    int jpos = -1;
    for (int i = 0; i < m; ++i) {
      if (mata[i * n + j] > 0.0) {
        if (jpos >= 0) {
          /* Second positive number, ergo not basic: */
          jpos = -1;
          break;
        }

        /* First positive number, could be basic. */
        jpos = i;
      } else if (mata[i * n + j] < 0.0) {
        /* Negative number, not basic: */
        jpos = -1;
        break;
      }
    }

    if (jpos >= 0) {
      /* Basic variable. */
      basis[jpos] = j;
      --nart; /* one less artificial variable required */
    }
  }

  /* nart could get negative when n > m */
  if (nart < 0) {
    nart = 0;
  }

  /* Tableau size: */
  int nrows = m + 2;
  int ncols = n + nart + 1;
  double* tableau = malloc(nrows * ncols * sizeof(double));

  /* Initial tableau:  n  nart  1
   *             m     A    I | b
   *             1    -c^T  0 | 0  minus phase 2 objective
   *             1     0    1 | 0  minus phase 1 objective
   */
  for (int i = 0; i < m; ++i) { /* A (m x n)*/
    for (int j = 0; j < n; ++j) {
      tableau[i * ncols + j] = mata[i * n + j];
    }
  }
  int iart = 0;
  for (int i = 0; i < m; ++i) { /* I (m x nart) */
    for (int j = 0; j < nart; ++j) {
      tableau[i * ncols + n + j] = 0.0;
    }
    if (basis[i] < 0) {
      tableau[i * ncols + n + iart] = 1.0;
      basis[i] = n + iart;
      ++iart;
    }
  }
  for (int i = 0; i < m; ++i) { /* b (m x 1)*/
    tableau[i * ncols + n + nart + 0] = b[i];
  }

  for (int j = 0; j < n; ++j) { /* -c^T (1 x n)*/
    tableau[m * ncols + j] = -c[j];
  }
  for (int j = 0; j < nart; ++j) { /* 0 (1 x nart) */
    tableau[m * ncols + n + j] = 0.0;
  }
  tableau[m * ncols + n + nart + 0] = 0.0; /* 0 (1 x 1) */

  for (int j = 0; j < n; ++j) { /* 0 (1 x n)*/
    tableau[(m + 1) * ncols + j] = 0.0;
  }
  for (int j = 0; j < nart; ++j) { /* 1 (1 x nart) */
    tableau[(m + 1) * ncols + n + j] = 1.0;
  }
  tableau[(m + 1) * ncols + n + nart + 0] = 0.0; /* 0 (1 x 1) */

  /* Make sure the basis coefficients are 1: */
  for (int i = 0; i < m; ++i) {
    vec_scale(ncols, &tableau[i * ncols], 1.0 / tableau[i * ncols + basis[i]]);
  }

  /* Make the tableau proper by subtracting each of the artificial rows from the
   * phase 1 objective: */
  for (int i = 0; i < m; ++i) {
    if (basis[i] >= n) {
      vec_sub(ncols, &tableau[(m + 1) * ncols], &tableau[i * ncols]);
    }
  }

  int ret = 0;
  int niter = 0;
  int isphase1 = nart > 0 ? 1 : 0; /* go directly to phase II if nart == 0 */
  while (1) {
    ++niter;

#ifdef PRINT
    printf("phase %s:\n", isphase1 ? "I" : "II");
    mat_print(nrows, ncols, tableau);
    printf("basis:");
    for (i = 0; i < m; ++i) {
      printf(" %d", basis[i]);
    }
    printf("\n");
#endif

    /* Entering variable (pivot column) is most negative entry in objective
     * function: */
    int cpiv = vec_minindex(n, &tableau[(m + isphase1) * ncols]);
    if (tableau[(m + isphase1) * ncols + cpiv] > -EPS) {
      if (isphase1) {
        /* No negative entry found in phase I objective. */
        if (tableau[(m + isphase1) * ncols + n + nart] < -EPS) {
          /* Could not get the artificial variables to zero: infeasible. */
          ret = -2;
          break;
        }

        /* Feasible point found, switch to phase II: */
        isphase1 = 0;
        continue;
      } else {
        /* No negative entry found in phase II objective, done. */
        break;
      }
    }

    /* Leaving / blocking variable (pivot row) is given by the minimum ratio of
     * rhs/pivotcolumn: */
    double minrat = HUGE_VAL;
    int rpiv = -1;
    for (int i = 0; i < m; ++i) {
      /* Denominator should be > 0: */
      if (tableau[i * ncols + cpiv] < EPS) {
        continue;
      }

      /* Ratio: */
      double rat = tableau[i * ncols + n + nart] / tableau[i * ncols + cpiv];
      if (rat < minrat) {
        minrat = rat;
        rpiv = i;
      }
    }

#ifdef PRINT
    printf("rpiv = %d, cpiv = %d\n\n", rpiv, cpiv);
#endif

    if (rpiv == -1) {
      /* No pivot row found (unbounded), exit. */
      ret = -1;
      break;
    }

    /* Perform pivot: */
    /* Scale row rpiv: */
    vec_scale(ncols, &tableau[rpiv * ncols],
              1.0 / tableau[rpiv * ncols + cpiv]);

    /* Subtract row rpiv from others: */
    for (int i = 0; i < nrows; ++i) {
      if (i == rpiv) {
        continue;
      }
      vec_adds(ncols, &tableau[i * ncols], &tableau[rpiv * ncols],
               -tableau[i * ncols + cpiv]);
    }

    /* cpiv enters, rpiv leaves basis: */
    basis[rpiv] = cpiv;
  }

  /* Extract solution from tableau: */
  memset(x, 0, n * sizeof(double));
  for (int i = 0; i < m; ++i) {
    x[basis[i]] = tableau[i * ncols + n + nart];
  }

  /* Clean up: */
  free(tableau);

  return ret;
}



int linprog_cn(int m, int n, double const* mata, double const* b,
               double const* c, double* x)
{
  /* Total size: */
  int ncols = n + n + m;

  /* Allocate augmented system (add mineq slack variables): */
  double* aprime = alloca(m * ncols * sizeof(double));
  double* bprime = alloca(m * sizeof(double));
  double* cprime = alloca(ncols * sizeof(double));
  double* xprime = alloca(ncols * sizeof(double));

  /*         n    n  m
   *  A' =   A   -A  I  m
   */
  for (int i = 0; i < m; ++i) { /* [A (m x n), -A (m x n)]*/
    for (int j = 0; j < n; ++j) {
      aprime[i * ncols + j] = mata[i * n + j];
      aprime[i * ncols + n + j] = -mata[i * n + j];
    }
  }
  for (int i = 0; i < m; ++i) { /* I (m x m)*/
    for (int j = 0; j < m; ++j) {
      aprime[i * ncols + n + n + j] = 0.0;
    }
    aprime[i * ncols + n + n + i] = 1.0;
  }

  /* Make sure all b_i >= 0: */
  for (int i = 0; i < m; ++i) {
    if (b[i] < 0.0) {
      bprime[i] = -b[i];
      vec_neg(ncols, &aprime[i * ncols]);
    } else {
      bprime[i] = b[i];
    }
  }

  /* c' = [c; -c; 0]: */
  memcpy(cprime, c, n * sizeof(double));
  memcpy(cprime + n, c, n * sizeof(double));
  vec_neg(n, cprime + n);
  memset(cprime + n + n, 0, m * sizeof(double));

  /* x' = 0_{ncols x 1}: */
  memset(xprime, 0, ncols * sizeof(double));

  /* Solve: */
  int ret = linprog_rn(m, ncols, aprime, bprime, cprime, xprime);

  /* Copy answer: */
  memcpy(x, xprime, n * sizeof(double)); /* x+ */
  vec_sub(n, x, &xprime[n]);             /* x+ - x- */

  return ret;
}
