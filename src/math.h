#pragma once

#include <stdio.h>

#define _USE_MATH_DEFINES
#include <math.h>
#undef _USE_MATH_DEFINES

#ifndef M_2PI
#define M_2PI (2.0 * M_PI)
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

/** Reset all components of a vector. */
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

/** Print a matrix to a buffer. */
void mat_sprint(int m, int n, double const* mat, char* str);

/** Print a matrix to stdout. */
void mat_print(int m, int n, double const* mat);

/** Print a matrix. */
void mat_fprint(int m, int n, double const* mat, FILE* outstream);

/** Axis-aligned bounding box of a point set. */
void boundingbox(int npoints, int ndims, double const* points, int* minindices,
                 int* maxindices, double* minima, double* maxima);

/** Compute volume, centroid and span of multidimensional simplex. */
void analysesimplex(int npoints, int ndims, double* points,
                    double* volume, double* centroid);

/** LQ = PA decomposition. Returns determinant. */
double lqdc(int m, int n, double* mat, int* p);

/** Solve linear program with Dantzig's simplex method.
    Restricted normal form.
    Maximimize c^T x
    subject to
    A x == b (m x n), b >= 0, x >= 0.

    return value ==  0: optimum found
    return value == -1: unbounded
    return value == -2: infeasible
*/
int linprog_rn(int m, int n, double const* mata, double const* b,
               double const* c, double* x);

/** Solve linear program with Dantzig's simplex method.
    Canonical form.
    Maximimize c^T x
      subject to
    A1 x == b1 (meq x n), A2 x <= b2 (mineq x n).
    With A = [A1; A2], b = [b1; b2].

    return value ==  0: optimum found
    return value == -1: unbounded
    return value == -2: infeasible
  */
int linprog_cn(int meq, int mineq, int n, double const* mata, double const* b,
               double const* c, double* x);
