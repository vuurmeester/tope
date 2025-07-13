#pragma once

#include <polytoop.h>

#include "allocator.h"

typedef struct _Ridge Ridge;
typedef struct _Point Point;
typedef struct _HashEntry HashEntry;

struct _HashEntry {
  unsigned hash;
  void* key;
  void* value;
  HashEntry* nextinchain;
  HashEntry* nextinlist;
  HashEntry* previnlist;
};

typedef struct _HashMap {
  int cap;                            /* current capacity */
  unsigned mask;                      /* for efficient modulo */
  int len;                            /* number of elements */
  void* data;                         /* userdata */
  unsigned (*hashfunc)(void*, void*); /* hash function */
  int (*compar)(void*, void*, void*); /* key comparison function */
  HashEntry** entries;                /* chain pointers */
  Allocator* allocator;
  HashEntry* first; /* first entry in linked list */
  HashEntry* last;  /* last entry in linked list */
} HashMap;

typedef struct _Array {
  int len;
  int cap;
  void** values;
} Array;

struct _Polytoop {
  Allocator* allocator;
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
  int merge;
  HashMap* newridges;
};

struct _polytoop_Facet {
  polytoop_Facet* next;
  polytoop_Facet* prev;

  Polytoop* polytoop;
  double volume;
  double* centroid;
  double* normal;     /* outward pointing plane normal */
  Array ridges;       /* d+ adjacent ridges */
  Array vertices;     /* d+ adjacent vertices */
  Point* outsidehead; /* visible points list */
  Point* outsidetail; /* last entry in visible points list */
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

  int d;
  int index;
  double height;
  double* pos;
};
