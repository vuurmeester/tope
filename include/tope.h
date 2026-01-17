#pragma once

#include <stdint.h>
#include <stdbool.h>

#ifdef _MSC_VER
  #ifdef BUILDING_TOPE
    #define TOPE_API __declspec(dllexport)
  #else
    #define TOPE_API __declspec(dllimport)
  #endif
#else
  #ifdef BUILDING_TOPE
    #define TOPE_API __attribute__ ((visibility("default")))
  #else
    #define TOPE_API
  #endif
#endif

typedef struct _Tope Tope;
typedef struct _tope_Facet tope_Facet;
typedef struct _tope_Vertex tope_Vertex;

/** Delete tope object. */
TOPE_API void tope_delete(Tope* tope);

/** Create tope from planes and interior point. */
TOPE_API Tope* tope_fromplanes(
  int nplanes,
  int dim,
  double const* normals,
  double const* distances,
  double const* xc
);

/** Create tope from point cloud. */
TOPE_API Tope* tope_frompoints(int npoints, int dim, double const* points, bool merge);

/** Delaunay triangulation. */
TOPE_API Tope* tope_delaunay(int npoints, int dim, double const* points);

/** Add vertex to existing tope. */
TOPE_API void tope_addvertex(Tope* tope, double const* point);

/** Print tope. */
TOPE_API void tope_print(Tope* tope);

/** Number of facets. */
TOPE_API int tope_getnumfacets(Tope* tope);

/** Number of ridges. */
TOPE_API int tope_getnumridges(Tope* tope);

/** Number of vertices. */
TOPE_API int tope_getnumvertices(Tope* tope);

/** Interpolate (when tope is Delaunay). */
TOPE_API void tope_interpolate(
  Tope* tope,
  double const* xi,
  int* indices,
  double* weights
);

/** First facet. */
TOPE_API tope_Facet* tope_firstfacet(Tope* tope);

/** Next facet. */
TOPE_API tope_Facet* tope_facet_nextfacet(tope_Facet* facet);

/** Facet normal. */
TOPE_API void tope_facet_getnormal(Tope* tope, tope_Facet* facet, double* normal);

/** Facet centroid. */
TOPE_API void tope_facet_getcentroid(Tope* tope, tope_Facet* facet, double* centroid);

/** Facet plane offset. */
TOPE_API double tope_facet_getoffset(Tope* tope, tope_Facet* facet);

/** Facet volume. */
TOPE_API double tope_facet_getvolume(Tope* tope, tope_Facet* facet);

/** Number of facet vertices. */
TOPE_API int tope_facet_getnumvertices(Tope* tope, tope_Facet* facet);

/** Vertex position. */
TOPE_API void tope_vertex_getposition(Tope* tope, tope_Vertex* vertex, double* position);

/** Memory in use. */
TOPE_API uint64_t tope_bytes_used(Tope* tope);
