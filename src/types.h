#pragma once

#include <tope.h>
#include <stdbool.h>

#include "allocator.h"

typedef struct _Ridge Ridge;
typedef struct _Point Point;
typedef tope_Vertex Vertex;
typedef tope_Facet Facet;

typedef struct _List List;
struct _List {
  List* tail;
  Vertex* vert;
  Ridge* ridge;
};

typedef struct _HashMap {
  u32 cap;     /* current capacity */
  u32 len;     /* number of elements */
  Ridge** ridges; /* entries */
  u32* hashes;
} HashMap;

struct _Tope {
  Allocator alc;
  int dim;
  bool isdelaunay;
  double* shift;
  double* scales;
  double* center;
  int nfacets;
  Facet* firstfacet;
  Facet* lastfacet;
  int nridges;
  int nverts;
  HashMap newridges;
  Ridge** horizonridges;
  int horizonridges_len;
  int horizonridges_cap;
  Facet** newfacets;
  int newfacets_len;
  int newfacets_cap;
};

struct _tope_Facet {
  Facet* next;
  Facet* prev;

  double volume;
  double dist;        /* distance from origin */
  Point* outsidehead; /* visible points list */
  Point* outsidetail; /* last entry in visible points list */
  bool visible;
  double* centroid;
  double* normal;
  List* rvs;  /* ridge/vertex list */
};

struct _Ridge {
  double* vdn;        /* volume, distance, normal */
  Facet* facets[2];   /* 2 adjacent facets */
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
  double height;
  double* pos;
};
