#pragma once

#define _USE_MATH_DEFINES
#include <math.h>
#undef _USE_MATH_DEFINES

#include <stdio.h>

#ifndef M_PI
#define M_PI 3.1415926535898
#endif

#ifndef M_2PI
#define M_2PI (2.0 * M_PI)
#endif

#ifndef M_SQRT2
#define M_SQRT2 1.4142135623731
#endif

#ifndef M_SQRT1_2
#define M_SQRT1_2 (0.5 * M_SQRT2)
#endif

/** Add a vector to another. */
void vec_add(int n, double* x, double const* y);

/** Subtract a vector from another. */
void vec_sub(int n, double* x, double const* y);

/** Add a scaled vector to another. */
void vec_adds(int n, double* x, double const* y, double scale);

/** Scale a vector. */
void vec_scale(int n, double* x, double scale);

/** Negate a vector. */
void vec_neg(int n, double* x);

/** Set all components of a vector. */
void vec_set(int n, double* x, double scalar);

/** Reset vector. */
void vec_reset(int n, double* x);

/** Inner product of two vectors. */
double vec_dot(int n, double const* x, double const* y);

/** Squared norm of a vector. */
double vec_nrmsq(int n, double const* x);

/** Norm of a vector. */
double vec_norm(int n, double const* x);

/** Normalize the vector. */
void vec_normalize(int n, double* x);

/** Sum of elements of a vector. */
double vec_sum(int n, double const* x);

/** Maximum element of a vector. */
int vec_maxindex(int n, double const* x);

/** Minimum element of a vector. */
int vec_minindex(int n, double const* x);

/** Print a vector to stdout. */
void vec_print(int n, double const* x);

/** Matrix-vector multiplication. */
void mat_vecmul(int m, int n, double const* mat, double const* x, double* y);

/** (m x n) times (n x o) matrix-matrix multiplication. */
void mat_matmul(
    int m,
    int n,
    int o,
    double const* mat1,
    double const* mat2,
    double* result
);

/** Print a matrix to a buffer. */
void mat_sprint(int m, int n, double const* mat, char* str);

/** Print a matrix to stdout. */
void mat_print(int m, int n, double const* mat);

/** Print a matrix. */
void mat_fprint(int m, int n, double const* mat, FILE* outstream);

/** Axis-aligned bounding box of a point set. */
void boundingbox(
    int npoints,
    int ndims,
    double const* points,
    int* minindices,
    int* maxindices,
    double* minima,
    double* maxima
);

/** Compute volume, centroid and span of multidimensional simplex. */
void analysesimplex(
    int npoints,
    int ndims,
    double* points,
    double* volume,
    double* centroid
);

/** Solve linear program.
    Maximimize c^T x
      subject to
    A x <= b (m x n).

    return value ==  0: optimum found
    return value == -1: unbounded
    return value == -2: infeasible
  */
int linprog(
    int m,
    int n,
    double const* A,
    double const* b,
    double const* c,
    double* x
);



/** Solve A x = b  nonsingular square system (A clobbered, b <- x).
    Return determinant. Gaussian LU eliminiation, partial pivoting. */
double gauss(int n, double* A, double* b);

/** Conjugate residual method for symmetric matrices. */
void cr(
    int n,
    void (*applymatrix)(int n, double const* x, double* y, void const* data),
    double const* b,
    double* x,
    double tol,
    void const* data
);