#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../src/math.h"
#include "../src/util.h"



int main()
{
  int itest;
  int i;
  int j;
  double d;

  double start = polytoop_gettime();
  int nunb = 0;
  int nopt = 0;
  int ntests = 1000;
  for (itest = 0; itest < ntests; ++itest) {
    int n = random_getint(2, 5);     /* dimension */
    int m = random_getint(n, 8 * n); /* number of planes */

    /* Normals plus extra plane, plus extra dimension: */
    double* A = malloc((m + 1) * (n + 1) * sizeof(double));
    double* b = malloc((m + 1) * sizeof(double));
    double* c = malloc((n + 1) * sizeof(double));
    double* x = malloc((n + 1) * sizeof(double));
    double* center = malloc(n * sizeof(double));

    /* Random center: */
    for (j = 0; j < n; ++j) {
      center[j] = 5.0 * (random_getdouble() - 0.5);
    }

    /* Distance from center: */
    d = 0.1 + random_getdouble();

    /* Planes: */
    for (i = 0; i < m; ++i) {
      /* Normal: */
      for (j = 0; j < n; ++j) {
        A[i * (n + 1) + j] = random_getdouble() - 0.5;
      }
      vec_normalize(n, A + i * (n + 1));

      /* Distance from origin: */
      b[i] = vec_dot(n, A + i * (n + 1), center) + d;
    }

    /* Extra dimension: */
    for (i = 0; i < m; ++i) {
      A[i * (n + 1) + n] = 1.0;
    }

    /* Extra plane to ensure boundedness: */
    memset(A + m * (n + 1), 0, n * sizeof(double));
    A[m * (n + 1) + n] = 1.0;
    b[m] = 2.0; /* ceiling at 'height' 2 */

    /* Optimization direction: (0, ..., 0, 1)^T: */
    memset(c, 0, n * sizeof(double));
    c[n] = 1;

    /* Do the optimization: */
    memset(x, 0, (n + 1) * sizeof(double));
    linprog(m + 1, n + 1, A, b, c, x);

    /* Check result: */
    if (x[n] > b[m]) {
      /* Should be below ceiling. */
      return -1;
    }
    if (fabs(x[n] - b[m]) < 1.0e-3) {
      /* Fine, unbounded. */
      ++nunb;
    } else {
      ++nopt;
      /* x should equal center: */
      vec_sub(n, x, center);
      if (vec_norm(n, x) > 1.0e-3) {
        return -1;
      }
    }

    /* Clean up iteration: */
    free(center);
    free(x);
    free(c);
    free(b);
    free(A);
  }

  printf("t = %g\n", polytoop_gettime() - start);
  printf("%d unbounded\n", nunb);
  printf("%d optimized\n", nopt);

  return 0;
}
