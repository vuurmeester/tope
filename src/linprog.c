#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "util.h"
#include "math.h"

#define EPS 1.0e-12



int linprog_rn(int m, int n, double const* mata, double const* b, double const* c, double* x)
{
  int i;
  int j;
  int rpiv;
  int cpiv;
  int nrows;
  int ncols;
  int jpos;
  int ret;
  int niter;
  int iart;
  int nart;
  int isphase1;
  int* basis;
  double rat;
  double minrat;
  double* tableau;

  basis = (int*)alloca(m * sizeof(int));
  for (i = 0; i < m; ++i) {
    basis[i] = -1;
  }

  /* Everything non-basic: */
  nart = m;

  /* Extract basis: */
  for (j = 0; j < n; ++j) {
    /* Search for column with exactly one positive number and otherwise zero: */
    jpos = -1;
    for (i = 0; i < m; ++i) {
      if (mata[i * n + j] > 0.0) {
        if (jpos >= 0) {
          /* Second positive number, ergo not basic: */
          jpos = -1;
          break;
        }

        /* First positive number, could be basic. */
        jpos = i;
      }
      else if (mata[i * n + j] < 0.0) {
        /* Negative number, not basic: */
        jpos = -1;
        break;
      }
    }

    if (jpos >= 0) {
      /* Basic variable. */
      basis[jpos] = j;
      --nart;  /* one less artificial variable required */
    }
  }

  /* nart could get negative when n > m */
  if (nart < 0) {
    nart = 0;
  }

  /* Tableau size: */
  nrows = m + 2;
  ncols = n + nart + 1;
  tableau = malloc(nrows * ncols * sizeof(double));

  /* Initial tableau:  n  nart  1
   *             m     A    I | b
   *             1    -c^T  0 | 0  minus phase 2 objective
   *             1     0    1 | 0  minus phase 1 objective
   */
  for (i = 0; i < m; ++i) {  /* A (m x n)*/
    for (j = 0; j < n; ++j) {
      tableau[i * ncols + j] = mata[i * n + j];
    }
  }
  iart = 0;
  for (i = 0; i < m; ++i) {  /* I (m x nart) */
    for (j = 0; j < nart; ++j) {
      tableau[i * ncols + n + j] = 0.0;
    }
    if (basis[i] < 0) {
      tableau[i * ncols + n + iart] = 1.0;
      basis[i] = n + iart;
      ++iart;
    }
  }
  for (i = 0; i < m; ++i) {  /* b (m x 1)*/
    tableau[i * ncols + n + nart + 0] = b[i];
  }

  for (j = 0; j < n; ++j) {  /* -c^T (1 x n)*/
    tableau[m * ncols + j] = -c[j];
  }
  for (j = 0; j < nart; ++j) {  /* 0 (1 x nart) */
    tableau[m * ncols + n + j] = 0.0;
  }
  tableau[m * ncols + n + nart + 0] = 0.0;  /* 0 (1 x 1) */

  for (j = 0; j < n; ++j) {  /* 0 (1 x n)*/
    tableau[(m + 1) * ncols + j] = 0.0;
  }
  for (j = 0; j < nart; ++j) {  /* 1 (1 x nart) */
    tableau[(m + 1) * ncols + n + j] = 1.0;
  }
  tableau[(m + 1) * ncols + n + nart + 0] = 0.0;  /* 0 (1 x 1) */

  /* Make sure the basis coefficients are 1: */
  for (i = 0; i < m; ++i) {
    vector_scale(ncols, &tableau[i * ncols], 1.0 / tableau[i * ncols + basis[i]]);
  }

  /* Make the tableau proper by subtracting each of the artificial rows from the
   * phase 1 objective: */
  for (i = 0; i < m; ++i) {
    if (basis[i] >= n) {
      vector_subtract(ncols, &tableau[(m + 1) * ncols], &tableau[i * ncols]);
    }
  }

  ret = 0;
  niter = 0;
  isphase1 = nart > 0 ? 1 : 0;  /* go directly to phase II if nart == 0 */
  while (1) {
    ++niter;

#ifdef PRINT
    printf("phase %s:\n", isphase1 ? "I" : "II");
    matrix_print(nrows, ncols, tableau);
    printf("basis:");
    for (i = 0; i < m; ++i) {
      printf(" %d", basis[i]);
    }
    printf("\n");
#endif

    /* Entering variable (pivot column) is most negative entry in objective function: */
    cpiv = vector_minindex(n, &tableau[(m + isphase1) * ncols]);
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
      }
      else {
        /* No negative entry found in phase II objective, done. */
        break;
      }
    }

    /* Leaving / blocking variable (pivot row) is given by the minimum ratio of rhs/pivotcolumn: */
    minrat = DBL_MAX;
    rpiv = -1;
    for (i = 0; i < m; ++i) {
      /* Denominator should be > 0: */
      if (tableau[i * ncols + cpiv] < EPS) {
        continue;
      }

      /* Ratio: */
      rat = tableau[i * ncols + n + nart] / tableau[i * ncols + cpiv];
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
    vector_scale(ncols, &tableau[rpiv * ncols], 1.0 / tableau[rpiv * ncols + cpiv]);

    /* Subtract row rpiv from others: */
    for (i = 0; i < nrows; ++i) {
      if (i == rpiv) {
        continue;
      }
      vector_add_scaled(ncols, &tableau[i * ncols], &tableau[rpiv * ncols], -tableau[i * ncols + cpiv]);
    }

    /* cpiv enters, rpiv leaves basis: */
    basis[rpiv] = cpiv;
  }

  /* Extract solution from tableau: */
  vector_reset(n, x);
  for (i = 0; i < m; ++i) {
    x[basis[i]] = tableau[i * ncols + n + nart];
  }

  /* Clean up: */
  free(tableau);

  return ret;
}



int linprog_cn(int meq, int mineq, int n, double const* mata, double const* b, double const* c, double* x)
{
  int i;  /* generic counter */
  int j;  /* generic counter */
  int ret;  /* return value */
  int nrows;  /* total number of equations */
  int ncols;  /* total number of variables */
  double* aprime;  /* augmented matrix with slack variables */
  double* bprime;  /* augmented rhs (negative for equality constraints) */
  double* cprime;  /* augmented objective function (zero for slack variables) */
  double* xprime;  /* augmented solution (with slack variables) */

  /* Total size: */
  nrows = meq + mineq;
  ncols = n + n + mineq;

  /* Allocate augmented system (add mineq slack variables): */
  aprime = alloca(nrows * ncols * sizeof(double));
  bprime = alloca(nrows * sizeof(double));
  cprime = alloca(ncols * sizeof(double));
  xprime = alloca(ncols * sizeof(double));

  /*         n    n  mineq
   *  A' =  A1  -A1      0  meq
   *        A2  -A2      I  mineq
   */
  for (i = 0; i < nrows; ++i) {  /* [A (m x n), -A (m x n)]*/
    for (j = 0; j < n; ++j) {
      aprime[i * ncols + j] = mata[i * n + j];
      aprime[i * ncols + n + j] = -mata[i * n + j];
    }
  }
  for (i = 0; i < meq; ++i) {  /* 0 (meq x mineq) */
    for (j = 0; j < mineq; ++j) {
      aprime[i * ncols + n + n + j] = 0.0;
    }
  }
  for (i = 0; i < mineq; ++i) {  /* I (mineq x mineq)*/
    for (j = 0; j < mineq; ++j) {
      aprime[(meq + i) * ncols + n + n + j] = 0.0;
    }
    aprime[(meq + i) * ncols + n + n + i] = 1.0;
  }

  /* Make sure all b_i >= 0: */
  for (i = 0; i < nrows; ++i) {
    if (b[i] < 0.0) {
      bprime[i] = -b[i];
      vector_negate(ncols, &aprime[i * ncols]);
    }
    else {
      bprime[i] = b[i];
    }
  }

  /* c' = [c; -c; 0]: */
  memcpy(cprime, c, n * sizeof(double));
  memcpy(&cprime[n], c, n * sizeof(double));
  vector_negate(n, &cprime[n]);
  vector_reset(mineq, &cprime[n + n]);

  /* x' = 0_{ncols x 1}: */
  vector_reset(ncols, xprime);

  /* Solve: */
  ret = linprog_rn(nrows, ncols, aprime, bprime, cprime, xprime);

  /* Copy answer: */
  memcpy(x, xprime, n * sizeof(double));  /* x+ */
  vector_subtract(n, x, &xprime[n]);  /* x+ - x- */

  return ret;
}
