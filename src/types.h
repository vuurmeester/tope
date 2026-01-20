#pragma once

#include <tope.h>

#include <stdbool.h>

#include "allocator.h"

typedef struct _Ridge Ridge;
typedef struct _Point Point;
typedef tope_Vertex Vertex;
typedef tope_Facet Facet;

/* Non-intrusive singly linked list. */
typedef struct _List List;
struct _List {
  List* next;  /* next in list */
  void* val;  /* generic pointer */
};

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
  HashMap newridges;  /* used to keep track of new ridges */
  Ridge** horizonridges;  /* used to keep track of horizon ridges */
  int horizonridges_len;
  int horizonridges_cap;
  Facet** newfacets;  /* used to keep track of new facets */
  int newfacets_len;
  int newfacets_cap;
};

struct _tope_Facet {
  Facet* next;  /* intrusive doubly linked list */
  Facet* prev;

  double volume;
  Point* outsidehead; /* first entry in visible points list */
  Point* outsidetail; /* last entry in visible points list */
  double* centroid;
  double* normal;
  List* verts;  /* vertex list */
  List* ridges;  /* ridge list */
};

struct _Ridge {
  Facet* facets[2];    /* 2 adjacent facets */
  Vertex* vertices[1]; /* d - 1 adjacent vertices */
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
