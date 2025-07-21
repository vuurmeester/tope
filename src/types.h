#pragma once

#include <polytoop.h>

#include "allocator.h"

typedef struct _Ridge Ridge;
typedef struct _Point Point;

typedef struct _HashMap {
  int cap;        /* current capacity */
  int len;        /* number of elements */
  Ridge** ridges; /* entries */
  unsigned* hashes;
} HashMap;

struct _Polytoop {
  Allocator alc;
  int dim;
  int isdelaunay;
  double* shift;
  double* scales;
  double* center;
  int nfacets;
  polytoop_Facet* firstfacet;
  polytoop_Facet* lastfacet;
  int nridges;
  Ridge* firstridge;
  Ridge* lastridge;
  int nverts;
  polytoop_Vertex* firstvertex;
  HashMap newridges;
  Ridge** horizonridges;
  int horizonridges_len;
  int horizonridges_cap;
  polytoop_Facet** newfacets;
  int newfacets_len;
  int newfacets_cap;
};

struct _polytoop_Facet {
  polytoop_Facet* next;
  polytoop_Facet* prev;

  Polytoop* polytoop;
  double volume;
  double* centroid;
  double* normal;             /* outward pointing plane normal */
  double dist;                /* distance from origin */
  Ridge** ridges;             /* d adjacent ridges */
  polytoop_Vertex** vertices; /* d adjacent vertices */
  Point* outsidehead;         /* visible points list */
  Point* outsidetail;         /* last entry in visible points list */
  int visible;
};

struct _Ridge {
  Ridge* next;
  Ridge* prev;

  double volume;
  double* centroid;
  double* normal;
  polytoop_Facet* facets[2];  /* 2 adjacent facets */
  polytoop_Vertex** vertices; /* d - 1 adjacent vertices */
};

struct _polytoop_Vertex {
  polytoop_Vertex* next;
  polytoop_Vertex* prev;

  Polytoop* polytoop;
  int index;
  double* position;
  int nridges; /* the number of ridges attached to this vertex */
};

struct _Point {
  Point* next;

  int index;
  double height;
  double* pos;
};
