#pragma once

#include <tope.h>
#include <stdbool.h>

#include "allocator.h"

typedef struct _Ridge Ridge;
typedef struct _Point Point;

typedef struct _HashMap {
  u32 cap;     /* current capacity */
  u32 len;     /* number of elements */
  u32* ridges; /* entries */
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
  u32 firstfacet;
  u32 lastfacet;
  int nridges;
  u32 firstridge;
  u32 lastridge;
  int nverts;
  u32 firstvertex;
  HashMap newridges;
  u32* horizonridges;
  int horizonridges_len;
  int horizonridges_cap;
  u32* newfacets;
  int newfacets_len;
  int newfacets_cap;
};

struct _tope_Facet {
  u32 next;
  u32 prev;

  Tope* tope;
  double volume;
  double dist;        /* distance from origin */
  Point* outsidehead; /* visible points list */
  Point* outsidetail; /* last entry in visible points list */
  bool visible;
  u32 hcentroid;
  u32 hnormal;
  u32 hridges;
  u32 hvertices;
};

struct _Ridge {
  u32 next;
  u32 prev;

  u32 hvdn;        /* handle to volume, distance and normal */
  u32 facets[2];   /* 2 adjacent facets */
  u32 vertices[1]; /* d - 1 adjacent vertices */
};

struct _tope_Vertex {
  u32 next;
  u32 prev;

  Tope* tope;
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
