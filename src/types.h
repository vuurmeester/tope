#pragma once

#include <tope.h>
#include <stdbool.h>

#include "allocator.h"

typedef struct _Ridge Ridge;
typedef struct _Point Point;
typedef tope_Vertex Vertex;
typedef tope_Facet Facet;

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
  Ridge* firstridge;
  Ridge* lastridge;
  int nverts;
  Vertex* firstvertex;
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

  Tope* tope;
  double volume;
  double dist;        /* distance from origin */
  Point* outsidehead; /* visible points list */
  Point* outsidetail; /* last entry in visible points list */
  bool visible;
  double* centroid;
  double* normal;
  Ridge** ridges;
  Vertex** vertices;
};

struct _Ridge {
  Ridge* next;
  Ridge* prev;

  double* vdn;        /* volume, distance, normal */
  Facet* facets[2];   /* 2 adjacent facets */
  Vertex* vertices[1]; /* d - 1 adjacent vertices */
};

struct _tope_Vertex {
  Vertex* next;
  Vertex* prev;

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
