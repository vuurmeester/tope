#pragma once

#include <polytoop.h>

#include "allocator.h"

typedef struct _Ridge Ridge;
typedef struct _Point Point;

typedef struct _HashMap {
  uint32_t cap; /* current capacity */
  uint32_t len; /* number of elements */
  u32* ridges;  /* entries */
  uint32_t* hashes;
} HashMap;

struct _Polytoop {
  Allocator alc;
  int dim;
  int isdelaunay;
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

struct _polytoop_Facet {
  u32 next;
  u32 prev;

  Polytoop* polytoop;
  double volume;
  double dist;        /* distance from origin */
  Point* outsidehead; /* visible points list */
  Point* outsidetail; /* last entry in visible points list */
  int visible;
  double centroid[1];   /* d-vector
  /*double* normal;*/   /* outward pointing plane normal */
  /*u32* ridges;   */   /* d adjacent ridges */
  /*u32* vertices; */   /* d adjacent vertices */
};

struct _Ridge {
  u32 next;
  u32 prev;

  double volume;
  u32 facets[2];   /* 2 adjacent facets */
  u32 vertices[1]; /* d - 1 adjacent vertices */
};

struct _polytoop_Vertex {
  u32 next;
  u32 prev;

  Polytoop* polytoop;
  int index;
  int nridges; /* the number of ridges attached to this vertex */
  double position[1];
};

struct _Point {
  Point* next;

  int index;
  double height;
  double* pos;
};
