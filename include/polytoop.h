#ifndef __polytoop_h__
#define __polytoop_h__

#ifdef __cplusplus
extern "C" {
#endif

typedef struct _Polytoop Polytoop;
typedef struct _polytoop_Facet polytoop_Facet;
typedef struct _polytoop_Vertex polytoop_Vertex;

/** Create polytoop object. */
Polytoop* polytoop_new(void);

/** Delete polytoop object. */
void polytoop_delete(Polytoop* polytoop);

/** Set facet merging. */
void polytoop_setmerge(Polytoop* polytoop, int merge);

/** Clear polytoop: */
void polytoop_clear(Polytoop* polytoop);

/** Create polytoop from planes. */
void polytoop_fromplanes(Polytoop* polytoop, int nplanes, int dim, double* normals, double* distances);

/** Create polytoop from point cloud. */
void polytoop_frompoints(Polytoop* polytoop, int npoints, int dim, double* points);

/** Delaunay triangulation. */
void polytoop_delaunay(Polytoop* polytoop, int npoints, int dim, double* points);

/** Add vertex to existing polytoop (increment). */
void polytoop_addvertex(Polytoop* polytoop, double* point);

/** Print polytoop. */
void polytoop_print(Polytoop* polytoop);

/** Number of facets. */
int polytoop_getnumfacets(Polytoop* polytoop);

/** Number of ridges. */
int polytoop_getnumridges(Polytoop* polytoop);

/** Total number of vertices. */
int polytoop_getnumvertices(Polytoop* polytoop);

/** Get vertex index. */
int polytoop_getvertexindex(Polytoop* polytoop, polytoop_Facet* facet, int ivertex);

/** Interpolate if polytoop is Delaunay grid. */
void polytoop_interpolate(Polytoop* polytoop, double const* xi, int* indices, double* weights);

/** First facet. */
polytoop_Facet* polytoop_firstfacet(Polytoop* polytoop);

/** First vertex. */
polytoop_Vertex* polytoop_firstvertex(Polytoop* polytoop);



/** Next facet. */
polytoop_Facet* polytoop_facet_nextfacet(polytoop_Facet* facet);

/** Number of neighbouring facets. */
int polytoop_facet_getnumneighbours(polytoop_Facet* facet);

/** Get facet neighbour. */
polytoop_Facet* polytoop_facet_getneighbour(polytoop_Facet* facet, int i);

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

/** Get facet vertex. */
polytoop_Vertex* polytoop_facet_getvertex(polytoop_Facet* facet, int ivertex);



/** Next vertex. */
polytoop_Vertex* polytoop_vertex_nextvertex(polytoop_Vertex* vertex);

/** Vertex position. */
void polytoop_vertex_getposition(polytoop_Vertex* vertex, double* position);

#ifdef __cplusplus
}
#endif

#endif  /* __polytoop_h__ */
