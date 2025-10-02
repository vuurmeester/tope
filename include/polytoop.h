#pragma once

typedef struct _Polytoop Polytoop;
typedef struct _polytoop_Facet polytoop_Facet;
typedef struct _polytoop_Vertex polytoop_Vertex;

/** Delete polytoop object. */
void polytoop_delete(Polytoop* polytoop);

/** Create polytoop from planes and interior point. */
Polytoop* polytoop_fromplanes(
    int nplanes,
    int dim,
    double const* normals,
    double const* distances,
    double const* xc
);

/** Create polytoop from point cloud. */
Polytoop* polytoop_frompoints(int npoints, int dim, double const* points);

/** Delaunay triangulation. */
Polytoop* polytoop_delaunay(int npoints, int dim, double const* points);

/** Add vertex to existing polytoop (increment). */
void polytoop_addvertex(Polytoop* polytoop, double const* point);

/** Print polytoop. */
void polytoop_print(Polytoop* polytoop);

/** Number of facets. */
int polytoop_getnumfacets(Polytoop* polytoop);

/** Number of ridges. */
int polytoop_getnumridges(Polytoop* polytoop);

/** Total number of vertices. */
int polytoop_getnumvertices(Polytoop* polytoop);

/** Interpolate if polytoop is Delaunay grid. */
void polytoop_interpolate(
    Polytoop* polytoop,
    double const* xi,
    int* indices,
    double* weights
);

/** First facet. */
polytoop_Facet* polytoop_firstfacet(Polytoop* polytoop);

/** First vertex. */
polytoop_Vertex* polytoop_firstvertex(Polytoop* polytoop);

/** Next facet. */
polytoop_Facet* polytoop_facet_nextfacet(polytoop_Facet* facet);

/** Facet normal. */
void polytoop_facet_getnormal(polytoop_Facet* facet, double* normal);

/** Facet centroid. */
void polytoop_facet_getcentroid(polytoop_Facet* facet, double* centroid);

/** Facet plane offset. */
double polytoop_facet_getoffset(polytoop_Facet* facet);

/** Facet volume. */
double polytoop_facet_getvolume(polytoop_Facet* facet);

/** Number of facet vertices. */
int polytoop_facet_getnumvertices(polytoop_Facet* facet);

/** Next vertex. */
polytoop_Vertex* polytoop_vertex_nextvertex(polytoop_Vertex* vertex);

/** Vertex position. */
void polytoop_vertex_getposition(polytoop_Vertex* vertex, double* position);
