#pragma once

#include <tope.h>

#include <stdbool.h>

#include "allocator.h"
#include "list.h"

typedef struct _Ridge Ridge;
typedef struct _Point Point;
typedef tope_Vertex Vertex;
typedef tope_Facet Facet;

typedef struct _HashMap {
  u32 cap;  /* capacity (2^n) */
  u32 len;  /* number of elements */
  Ridge** ridges;  /* entries */
  u32* hashes;  /* hash cache */
} HashMap;

struct _Tope {
  Allocator alc;  /* allocator */
  int dim;  /* tope dimension */
  bool isdelaunay;  /* Delaunay triangulation? */
  double* shift;  /* offset applied to input points */
  double* scales;  /* scale applied to input points */
  double* center;  /* interior point */
  int nfacets;  /* facet count */
  Facet* firstfacet;  /* first facet */
  Facet* lastfacet;  /* last facet */
  int nridges;  /* ridge count */
  int nverts;  /* vertex count */
  int maxindex;  /* max vertex index */
  HashMap newridges;  /* used to keep track of new ridges */
};

struct _tope_Facet {
  Facet* next;  /* intrusive doubly linked list */
  Facet* prev;

  double size;
  Point* outsidehead; /* first entry in visible points list */
  Point* outsidetail; /* last entry in visible points list */
  double* centroid;
  double* normal;
  List* ridges;  /* ridge list */
};

struct _Ridge {
  Facet* facets[2];  /* 2 adjacent facets */
  List* verts;  /* vertex list */
  double size;
  double* centroid;  /* d-dimensional vector */
  double* normal1;  /* d-dimensional normal vector */
  double* normal2;  /* d-dimensional normal vector */
};

struct _tope_Vertex {
  int index;
  int nridges;        /* the number of ridges attached to this vertex */
  double position[1]; /* d-vector */
};

struct _Point {
  Point* next;

  int index;
  double height;  /* height above facet which contains this point in 'outsidepoints' */
  double* pos;  /* d-vector */
};
