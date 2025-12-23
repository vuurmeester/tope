#pragma once

#include <stdint.h>

typedef struct _Tope Tope;
typedef struct _tope_Facet tope_Facet;
typedef struct _tope_Vertex tope_Vertex;

/** Delete tope object. */
void tope_delete(Tope* tope);

/** Create tope from planes and interior point. */
Tope* tope_fromplanes(
  int nplanes,
  int dim,
  double const* normals,
  double const* distances,
  double const* xc
);

/** Create tope from point cloud. */
Tope* tope_frompoints(int npoints, int dim, double const* points);

/** Delaunay triangulation. */
Tope* tope_delaunay(int npoints, int dim, double const* points);

/** Add vertex to existing tope. */
void tope_addvertex(Tope* tope, double const* point);

/** Print tope. */
void tope_print(Tope* tope);

/** Number of facets. */
int tope_getnumfacets(Tope* tope);

/** Number of ridges. */
int tope_getnumridges(Tope* tope);

/** Number of vertices. */
int tope_getnumvertices(Tope* tope);

/** Interpolate (when tope is Delaunay). */
void tope_interpolate(
  Tope* tope,
  double const* xi,
  int* indices,
  double* weights
);

/** First facet. */
tope_Facet* tope_firstfacet(Tope* tope);

/** First vertex. */
tope_Vertex* tope_firstvertex(Tope* tope);

/** Next facet. */
tope_Facet* tope_facet_nextfacet(tope_Facet* facet);

/** Facet normal. */
void tope_facet_getnormal(Tope* tope, tope_Facet* facet, double* normal);

/** Facet centroid. */
void tope_facet_getcentroid(Tope* tope, tope_Facet* facet, double* centroid);

/** Facet plane offset. */
double tope_facet_getoffset(Tope* tope, tope_Facet* facet);

/** Facet volume. */
double tope_facet_getvolume(Tope* tope, tope_Facet* facet);

/** Number of facet vertices. */
int tope_facet_getnumvertices(Tope* tope, tope_Facet* facet);

/** Next vertex. */
tope_Vertex* tope_vertex_nextvertex(tope_Vertex* vertex);

/** Vertex position. */
void tope_vertex_getposition(Tope* tope, tope_Vertex* vertex, double* position);

/** Memory used. */
uint64_t tope_bytes_used(Tope* tope);
