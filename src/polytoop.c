#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "hashmap.h"
#include "math.h"
#include "util.h"

#include <polytoop.h>

#define EPS 1.0e-9

typedef struct _Ridge Ridge;
typedef struct _Point Point;

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



static Array array_new(void)
{
  Array arr;
  arr.len = 0;
  arr.cap = 0;
  arr.values = NULL;
  return arr;
}



static void array_delete(Array* arr, Allocator* alc)
{
  allocator_free(alc, arr->values, arr->cap * sizeof(void*));
  memset(arr, 0, sizeof(Array));
}



static void array_append(Array* arr, void* value, Allocator* alc)
{
  if (arr->len == arr->cap) {
    int newcap = arr->cap > 0 ? 2 * arr->cap : 4;
    arr->values =
        allocator_realloc(alc, arr->values, arr->cap * sizeof(void*),
                          newcap * sizeof(void*));
    arr->cap = newcap;
  }
  arr->values[arr->len] = value;
  ++arr->len;
}



static void array_remove(Array* arr, int i)
{
  --arr->len;
  if (i < arr->len) {
    arr->values[i] = arr->values[arr->len];
  }
}



static int array_find(Array arr, void* value)
{
  int i;
  for (i = 0; i < arr.len; ++i) {
    if (arr.values[i] == value) {
      return i;
    }
  }
  return -1;
}



/* Create and initialize a vertex struct with a position and id: */
static polytoop_Vertex* vertex_new(Polytoop* polytoop, double* pos, int index)
{
  polytoop_Vertex* vertex;

  /* Allocate vertex: */
  vertex = allocator_alloc(polytoop->allocator, sizeof(polytoop_Vertex));

  /* Prepend to list: */
  vertex->next = polytoop->firstvertex;
  vertex->prev = NULL;
  if (vertex->next) {
    vertex->next->prev = vertex;
  }
  polytoop->firstvertex = vertex;
  ++polytoop->nverts;

  /* Fill in other variables: */
  vertex->polytoop = polytoop;
  vertex->index = index;
  vertex->position =
      allocator_alloc(polytoop->allocator, polytoop->dim * sizeof(double));
  memcpy(vertex->position, pos, polytoop->dim * sizeof(double));
  vertex->nridges = 0;

  return vertex;
}



static void vertex_free(Polytoop* polytoop, polytoop_Vertex* vertex)
{
  allocator_free(polytoop->allocator, vertex->position,
                 polytoop->dim * sizeof(double));
  allocator_free(polytoop->allocator, vertex, sizeof(polytoop_Vertex));
}



static void vertex_remove(Polytoop* polytoop, polytoop_Vertex* vertex)
{
  polytoop_Vertex* next;
  polytoop_Vertex* prev;

  next = vertex->next;
  prev = vertex->prev;
  if (next) {
    next->prev = prev;
  }
  if (prev) {
    prev->next = next;
  } else {
    polytoop->firstvertex = next;
  }
  --polytoop->nverts;

  vertex_free(polytoop, vertex);
}



/* Create and initialize a ridge struct with vertices: */
static Ridge* create_ridge(Polytoop* polytoop, polytoop_Vertex** vertices)
{
  int i;
  Ridge* ridge;

  /* Create ridge: */
  ridge = allocator_alloc(polytoop->allocator, sizeof(Ridge));
  ridge->next = NULL;
  ridge->prev = NULL;
  ridge->volume = 0.0;
  ridge->centroid = NULL;
  ridge->normal = NULL;
  ridge->facets[0] = NULL;
  ridge->facets[1] = NULL;
  ridge->vertices = allocator_alloc(
      polytoop->allocator, (polytoop->dim - 1) * sizeof(polytoop_Vertex*));
  memcpy(ridge->vertices, vertices,
         (polytoop->dim - 1) * sizeof(polytoop_Vertex*));

  /* Associate vertices with ridge: */
  for (i = 0; i < polytoop->dim - 1; ++i) {
    ++vertices[i]->nridges;
  }

  /* Append ridge to polytoop list: */
  ridge->prev = polytoop->lastridge;
  ridge->next = NULL;
  polytoop->lastridge = ridge;
  if (ridge->prev) {
    ridge->prev->next = ridge;
  } else {
    polytoop->firstridge = ridge;
  }
  ++polytoop->nridges;

  return ridge;
}



/* Create and initialize a facet struct: */
static polytoop_Facet* facet_new(Polytoop* polytoop)
{
  polytoop_Facet* facet;

  facet = allocator_alloc(polytoop->allocator, sizeof(polytoop_Facet));
  ++polytoop->nfacets;

  /* Append to list: */
  facet->prev = polytoop->lastfacet;
  if (facet->prev) {
    facet->prev->next = facet;
  } else {
    polytoop->firstfacet = facet;
  }
  polytoop->lastfacet = facet;
  facet->next = NULL;

  facet->polytoop = polytoop;
  facet->volume = 0.0;
  facet->centroid =
      allocator_alloc(polytoop->allocator, polytoop->dim * sizeof(double));
  vec_reset(polytoop->dim, facet->centroid);
  facet->normal =
      allocator_alloc(polytoop->allocator, polytoop->dim * sizeof(double));
  vec_reset(polytoop->dim, facet->normal);
  facet->ridges = array_new();
  facet->vertices = array_new();
  facet->outsidehead = NULL;
  facet->outsidetail = NULL;
  facet->visible = 0;

  return facet;
}



static void facet_free(Polytoop* polytoop, polytoop_Facet* facet)
{
  array_delete(&facet->vertices, polytoop->allocator);
  array_delete(&facet->ridges, polytoop->allocator);
  allocator_free(polytoop->allocator, facet->normal,
                 polytoop->dim * sizeof(double));
  allocator_free(polytoop->allocator, facet->centroid,
                 polytoop->dim * sizeof(double));
  allocator_free(polytoop->allocator, facet, sizeof(polytoop_Facet));
}



static void facet_remove(Polytoop* polytoop, polytoop_Facet* facet)
{
  polytoop_Facet* next;
  polytoop_Facet* prev;
  next = facet->next;
  prev = facet->prev;

  if (next) {
    next->prev = prev;
  } else {
    assert(polytoop->lastfacet == facet);
    polytoop->lastfacet = prev;
  }

  if (prev) {
    prev->next = next;
  } else {
    assert(polytoop->firstfacet == facet);
    polytoop->firstfacet = next;
  }
  --polytoop->nfacets;
  facet_free(polytoop, facet);
}



static void facet_addoutside(Polytoop* polytoop, polytoop_Facet* facet,
                             Point* point)
{
  if (!facet->outsidehead) {
    /* Initialize list: */
    facet->outsidehead = point;
    facet->outsidetail = point;
    point->next = NULL;
  } else if (point->height > facet->outsidehead->height) {
    /* Prepend: */
    point->next = facet->outsidehead;
    facet->outsidehead = point;
  } else {
    /* Append: */
    facet->outsidetail->next = point;
    facet->outsidetail = point;
    point->next = NULL;
  }

  /* Move facet to front: */
  if (facet->prev) {
    /* Unlink from facet list: */
    if (facet->next) {
      facet->next->prev = facet->prev;
    } else {
      assert(polytoop->lastfacet == facet);
      polytoop->lastfacet = facet->prev;
    }
    facet->prev->next = facet->next;

    if (polytoop->firstfacet->outsidehead &&
        polytoop->firstfacet->outsidehead->height >
            facet->outsidehead->height) {
      /* Move to second: */
      facet->next = polytoop->firstfacet->next;
      if (facet->next) {
        facet->next->prev = facet;
      }
      facet->prev = polytoop->firstfacet;
      facet->prev->next = facet;
    } else {
      /* Move to front: */
      facet->next = polytoop->firstfacet;
      facet->next->prev = facet;
      polytoop->firstfacet = facet;
      facet->prev = NULL;
    }
  }
}



static void ridge_free(Polytoop* polytoop, Ridge* ridge)
{
  if (polytoop->isdelaunay) {
    allocator_free(polytoop->allocator, ridge->normal,
                   (polytoop->dim - 1) * sizeof(double));
    allocator_free(polytoop->allocator, ridge->centroid,
                   (polytoop->dim - 1) * sizeof(double));
  }
  allocator_free(polytoop->allocator, ridge->vertices,
                 (polytoop->dim - 1) * sizeof(polytoop_Vertex*));
  allocator_free(polytoop->allocator, ridge, sizeof(Ridge));
}



static void ridge_remove(Polytoop* polytoop, Ridge* ridge)
{
  if (ridge->prev) {
    ridge->prev->next = ridge->next;
  } else {
    assert(polytoop->firstridge == ridge);
    polytoop->firstridge = ridge->next;
  }
  if (ridge->next) {
    ridge->next->prev = ridge->prev;
  } else {
    assert(polytoop->lastridge == ridge);
    polytoop->lastridge = ridge->prev;
  }
  --polytoop->nridges;

  ridge_free(polytoop, ridge);
}



/* Integrate facet into polytoop (compute volume, create ridges, etc): */
static void addfacet(Polytoop* polytoop, polytoop_Facet* facet, Ridge** ridges)
{
  int i;
  int j;
  double ip;
  double norm;
  double* dirs;
  Ridge* ridge;

  /* Facet centroid: */
  vec_reset(polytoop->dim, facet->centroid);
  for (i = 0; i < polytoop->dim; ++i) {
    vec_add(polytoop->dim, facet->centroid,
            ((polytoop_Vertex*)facet->vertices.values[i])->position);
  }
  vec_scale(polytoop->dim, facet->centroid, 1.0 / (double)polytoop->dim);

  /* Vertex span: */
  dirs = alloca((polytoop->dim - 1) * polytoop->dim * sizeof(double));
  for (i = 0; i < polytoop->dim - 1; ++i) {
    /* Vertex position difference: */
    memcpy(&dirs[i * polytoop->dim],
           ((polytoop_Vertex*)facet->vertices.values[i + 1])->position,
           polytoop->dim * sizeof(double));
    vec_sub(polytoop->dim, &dirs[i * polytoop->dim],
            ((polytoop_Vertex*)facet->vertices.values[0])->position);
  }

  /* Basis and volume: */
  facet->volume = 1.0;
  for (i = 0; i < polytoop->dim - 1; ++i) {
    double maxnrmsq;
    int pivot;

    /* Pivot largest row on top: */
    pivot = i;
    maxnrmsq = vec_nrmsq(polytoop->dim, &dirs[i * polytoop->dim]);
    for (j = i + 1; j < polytoop->dim - 1; ++j) {
      double nrmsq;
      nrmsq = vec_nrmsq(polytoop->dim, &dirs[j * polytoop->dim]);
      if (nrmsq > maxnrmsq) {
        maxnrmsq = nrmsq;
        pivot = j;
      }
    }
    if (pivot > i) {
      memswp(&dirs[i * polytoop->dim], &dirs[pivot * polytoop->dim],
             polytoop->dim * sizeof(double));
    }

    /* Normalize: */
    norm = sqrt(maxnrmsq);
    vec_scale(polytoop->dim, &dirs[i * polytoop->dim], 1.0 / norm);

    /* Orthogonalize: */
    for (j = i + 1; j < polytoop->dim - 1; ++j) {
      ip = vec_dot(polytoop->dim, &dirs[j * polytoop->dim],
                   &dirs[i * polytoop->dim]);
      vec_adds(polytoop->dim, &dirs[j * polytoop->dim],
               &dirs[i * polytoop->dim], -ip);
    }

    /* Accumulate area: */
    facet->volume *= norm / (double)(i + 1);
  }

  /* Compute normal: */
  memcpy(facet->normal, facet->centroid, polytoop->dim * sizeof(double));
  vec_sub(polytoop->dim, facet->normal, polytoop->center);
  for (i = 0; i < polytoop->dim - 1; ++i) {
    ip = vec_dot(polytoop->dim, facet->normal, &dirs[i * polytoop->dim]);
    vec_adds(polytoop->dim, facet->normal, &dirs[i * polytoop->dim], -ip);
  }
  vec_normalize(polytoop->dim, facet->normal);

  /* Assign ridges: */
  for (i = 0; i < polytoop->dim; ++i) {
    /* Get ridge: */
    ridge = ridges[i];

    /* Remember facet. */
    if (ridge->facets[0] == NULL) {
      ridge->facets[0] = facet;
    } else {
      assert(ridge->facets[1] == NULL);
      ridge->facets[1] = facet;
    }

    array_append(&facet->ridges, ridge, polytoop->allocator);
  }
}

static unsigned hashvertset(void* ptr, void* data)
{
  int i;
  int d;
  int nverts;
  unsigned ret;
  polytoop_Vertex** verts;

  /* Dimension: */
  d = *(int*)data;

  /* Vertices: */
  nverts = d - 1;
  verts = ptr;

  /* Sum vertex id's: */
  ret = 0;
  for (i = 0; i < nverts; ++i) {
    ret += verts[i]->index;
  }

  return ret;
}



/* Compare sets of d - 1 sorted vertices (for ordering ridges): */
static int compvertsets(void* ptr1, void* ptr2, void* data)
{
  int i;
  int d;
  polytoop_Vertex** vertset1;
  polytoop_Vertex** vertset2;

  d = *(int*)data;
  vertset1 = ptr1;
  vertset2 = ptr2;

  for (i = 0; i < d - 1; ++i) {
    if (vertset1[i] < vertset2[i]) {
      return -1;
    }
    if (vertset1[i] > vertset2[i]) {
      return 1;
    }
  }

  return 0;
}



/* Create initial simplex, add outside vertices, etc: */
static void initialsimplex(Polytoop* polytoop, int npoints, Point** points)
{
  int i;
  int j;
  int k;
  int minindex;
  int maxindex;
  int* p;
  double ip;
  double nrmsq;
  double maxnrmsq;
  double prevmaxnormsq;
  double* prevpoint;
  double* span;
  double* vec;
  polytoop_Facet* facet;
  polytoop_Vertex** vertices;
  polytoop_Vertex** ridgeverts;
  Ridge** facetridges;

  /* Allocations: */
  span = malloc((npoints - 1) * polytoop->dim * sizeof(double));
  p = malloc(npoints * sizeof(int));
  vec = alloca(polytoop->dim * sizeof(double));
  prevpoint = alloca(polytoop->dim * sizeof(double));

  /* Unit box center: */
  vec_set(polytoop->dim, prevpoint, 0.5);

  /* Point cloud extreme points (estimate): */
  minindex = -1;
  prevmaxnormsq = -HUGE_VAL;
  while (1) {
    maxindex = -1;
    maxnrmsq = -HUGE_VAL;
    for (i = 0; i < npoints; ++i) {
      nrmsq = 0.0;
      for (j = 0; j < polytoop->dim; ++j) {
        double d;
        d = points[i]->pos[j] - prevpoint[j];
        nrmsq += d * d;
      }
      if (nrmsq > maxnrmsq) {
        maxnrmsq = nrmsq;
        maxindex = i;
      }
    }
    if (maxnrmsq > prevmaxnormsq) {
      minindex = maxindex;
      memcpy(prevpoint, points[maxindex]->pos, polytoop->dim * sizeof(double));
      prevmaxnormsq = maxnrmsq;
    } else {
      if (minindex > maxindex) {
        i = maxindex;
        maxindex = minindex;
        minindex = i;
      }
      break;
    }
  }

  /* Place extreme points up front: */
  memswp(&points[0], &points[minindex], sizeof(Point*));
  memswp(&points[1], &points[maxindex], sizeof(Point*));

  /* Initialize span (positions w.r.t. first point): */
  for (i = 1; i < npoints; ++i) {
    memcpy(&span[(i - 1) * polytoop->dim], points[i]->pos,
           polytoop->dim * sizeof(double));
    vec_sub(polytoop->dim, &span[(i - 1) * polytoop->dim], points[0]->pos);
  }

  /* Perform LQ decomposition on span: */
  lqdc(npoints - 1, polytoop->dim, span, &p[1]);

  /* Due to pivoting, the permutation vector will represent the point order
     which leads to the greatest volume. */

  /* Adjust permutation vector to reflect point order (as opposed to 'span'
   * order): */
  p[0] = 0;
  for (i = 1; i < npoints; ++i) {
    ++p[i];
  }

  /* Initialize polytoop center to initial simplex centroid: */
  vec_reset(polytoop->dim, polytoop->center);
  for (i = 0; i < polytoop->dim + 1; ++i) {
    vec_add(polytoop->dim, polytoop->center, points[p[i]]->pos);
  }
  vec_scale(polytoop->dim, polytoop->center, 1.0 / (double)(polytoop->dim + 1));

  /* Create initial simplex vertices: */
  vertices = alloca((polytoop->dim + 1) * sizeof(polytoop_Vertex*));
  for (i = 0; i < polytoop->dim + 1; ++i) {
    vertices[i] = vertex_new(polytoop, points[p[i]]->pos, points[p[i]]->index);
  }

  facetridges = alloca(polytoop->dim * sizeof(Ridge*));
  ridgeverts = alloca((polytoop->dim - 1) * sizeof(polytoop_Vertex*));

  hashmap_clear(polytoop->newridges);

  /* Create facets: */
  for (i = 0; i < polytoop->dim + 1; ++i) {
    /* Create facet: */
    facet = facet_new(polytoop);

    /* Add all vertices except for one: */
    for (j = 0; j < polytoop->dim + 1; ++j) {
      if (j != i) {
        /* Add vertex position: */
        array_append(&facet->vertices, vertices[j], polytoop->allocator);
      }
    }

    /* Add ridges: */
    for (j = 0; j < polytoop->dim; ++j) {
      int size;

      /* Ridge verts: */
      size = 0;
      for (k = 0; k < polytoop->dim; ++k) {
        if (k != j) {
          ridgeverts[size++] = facet->vertices.values[k];
        }
      }

      facetridges[j] = hashmap_retrieve(polytoop->newridges, ridgeverts);
      if (!facetridges[j]) {
        facetridges[j] = create_ridge(polytoop, ridgeverts);
        hashmap_insert(polytoop->newridges, facetridges[j]->vertices,
                       facetridges[j]);
      }
    }

    /* Vertex array filled, compute and integrate the facet: */
    addfacet(polytoop, facet, facetridges);
  }

  /* Create outside sets (assign each remaining vertex to a facet which it
   * 'sees'): */
  for (i = polytoop->dim + 1; i < npoints; ++i) {
    /* Determine maximum distance: */
    for (facet = polytoop->firstfacet; facet != NULL; facet = facet->next) {
      /* Vector from facet centroid to vertex: */
      memcpy(vec, points[p[i]]->pos, polytoop->dim * sizeof(double));
      vec_sub(polytoop->dim, vec, facet->centroid);

      /* Orthogonal distance to facet: */
      ip = vec_dot(polytoop->dim, facet->normal, vec);

      /* If above facet, add it to outside set and move on. */
      if (ip > EPS) {
        points[p[i]]->height = ip;
        facet_addoutside(polytoop, facet, points[p[i]]);
        break;
      }
    }
  }

  /* Clean up: */
  free(p);
  free(span);
}



/* Add point to polytoop: */
static void addpoint(Polytoop* polytoop, polytoop_Facet* facet, Point* apex)
{
  int ifacet;
  int iridge;
  int iridge2;
  int ivertex;
  double ip;
  double* vec;
  polytoop_Vertex* vertex;
  polytoop_Facet* neighbour;
  Ridge* ridge;
  Ridge* horizonridge;
  Point* outsidepoints;
  Point* point;
  polytoop_Facet* visiblelist;
  Array newfacets;
  Array horizonridges;
  polytoop_Vertex** ridgeverts;
  Ridge** facetridges;
  Ridge* newridge;

  /* Allocations: */
  newfacets = array_new();
  horizonridges = array_new();
  vec = alloca(polytoop->dim * sizeof(double));

  /* Remove facet from polytoop list: */
  if (facet->prev) {
    facet->prev->next = facet->next;
  } else {
    polytoop->firstfacet = facet->next;
  }
  if (facet->next) {
    facet->next->prev = facet->prev;
  } else {
    polytoop->lastfacet = facet->prev;
  }
  --polytoop->nfacets;

  /* Add facet to visible list: */
  visiblelist = facet;
  facet->prev = NULL;
  facet->next = NULL;
  facet->visible = 1;

  /* Initialize outside points list: */
  outsidepoints = NULL;

  while (visiblelist) {
    /* Pop a value: */
    facet = visiblelist;
    visiblelist = facet->next;

    /* Loop over ridges to visit neighbours: */
    for (iridge = 0; iridge < facet->ridges.len; ++iridge) {
      /* Ridge i: */
      ridge = facet->ridges.values[iridge];

      /* Retrieve neighbour: */
      neighbour = NULL;
      if (facet == ridge->facets[0]) {
        /* Current facet is 0, neighbour is 1: */
        ridge->facets[0] = NULL; /* forget reference to this facet */
        neighbour = ridge->facets[1];
      } else {
        /* Current facet is 1, neighbour is 0: */
        ridge->facets[1] = NULL; /* forget reference to this facet */
        neighbour = ridge->facets[0];
      }

      if (!neighbour) {
        /* This neighbour was already removed, remove ridge. */

        /* Remove reference to this ridge from adjacent vertices: */
        for (ivertex = 0; ivertex < polytoop->dim - 1; ++ivertex) {
          vertex = ridge->vertices[ivertex];
          --vertex->nridges;
          if (vertex->nridges == 0) {
            /* Remove the vertex as well: */
            vertex_remove(polytoop, vertex);
          }
        }

        /* Remove the ridge from the polytoop container: */
        ridge_remove(polytoop, ridge);
      } else if (!neighbour->visible) {
        /* Neighbour is untested. */

        /* Vector from neighbour centroid to apex: */
        memcpy(vec, apex->pos, polytoop->dim * sizeof(double));
        vec_sub(polytoop->dim, vec, neighbour->centroid);

        /* Height of apex above neighbour: */
        ip = vec_dot(polytoop->dim, vec, neighbour->normal);

        /* Decide if neighbour is visible: */
        if (ip > EPS) {
          /* Remove neighbour from polytoop list: */
          if (neighbour->prev) {
            neighbour->prev->next = neighbour->next;
          } else {
            polytoop->firstfacet = neighbour->next;
          }
          if (neighbour->next) {
            neighbour->next->prev = neighbour->prev;
          } else {
            polytoop->lastfacet = neighbour->prev;
          }
          --polytoop->nfacets;

          /* Prepend neighbour to visible list: */
          neighbour->next = visiblelist;
          neighbour->prev = NULL; /* singly linked list */
          visiblelist = neighbour;
          neighbour->visible = 1;
        } else {
          /* Neighbour not visible. Ridge belongs to horizon: */
          array_append(&horizonridges, ridge, polytoop->allocator);
        }
      }
    }

    /* Remember outside vertices of investigated facet: */
    if (facet->outsidehead) {
      facet->outsidetail->next = outsidepoints;
      outsidepoints = facet->outsidehead;
    }

    /* Deallocate investigated facet: */
    facet_free(polytoop, facet);
  }

  /* Add apex vertex to polytoop: */
  vertex = vertex_new(polytoop, apex->pos, apex->index);

  ridgeverts = alloca((polytoop->dim - 1) * sizeof(polytoop_Vertex*));
  facetridges = alloca(polytoop->dim * sizeof(Ridge*));

  hashmap_clear(polytoop->newridges);

  /* Form new facets: */
  for (iridge = 0; iridge < horizonridges.len; ++iridge) {
    /* The horizon ridge: */
    horizonridge = horizonridges.values[iridge];

    /* New facet: */
    facet = facet_new(polytoop);

    /* Facet vertices are horizon vertices + apex: */
    for (ivertex = 0; ivertex < polytoop->dim - 1; ++ivertex) {
      array_append(&facet->vertices, horizonridge->vertices[ivertex],
                   polytoop->allocator);
    }
    array_append(&facet->vertices, vertex, polytoop->allocator);

    /* 1 existing (horizon) ridge and dim - 1 new ridges: */
    facetridges[0] = horizonridge;
    for (ivertex = 0; ivertex < polytoop->dim - 1; ++ivertex) {
      int size;
      int jvertex;

      /* Add all but one of the horizon ridge vertices: */
      size = 0;
      for (jvertex = 0; jvertex < polytoop->dim - 1; ++jvertex) {
        if (jvertex != ivertex) {
          ridgeverts[size++] = horizonridge->vertices[jvertex];
        }
      }

      /* Always add apex: */
      ridgeverts[size++] = vertex;

      /* Get or create new ridge: */
      newridge = hashmap_retrieve(polytoop->newridges, ridgeverts);
      if (!newridge) {
        /* Create new ridge: */
        newridge = create_ridge(polytoop, ridgeverts);
        hashmap_insert(polytoop->newridges, newridge->vertices, newridge);
      }
      facetridges[ivertex + 1] = newridge;
    }

    addfacet(polytoop, facet, facetridges);

    if (polytoop->merge) {
      /* Retrieve neighbour: */
      if (horizonridge->facets[0] == facet) {
        neighbour = horizonridge->facets[1];
      } else {
        neighbour = horizonridge->facets[0];
      }

      /* Normal vector difference: */
      memcpy(vec, facet->normal, polytoop->dim * sizeof(double));
      vec_sub(polytoop->dim, vec, neighbour->normal);

      if (vec_nrmsq(polytoop->dim, vec) < EPS * EPS) {
        double newvol;
        Ridge* ridge2;

        /* Parallel, facet is to be removed. */

        /* Add facet ridges to neighbour: */
        for (iridge2 = 0; iridge2 < facet->ridges.len; ++iridge2) {
          ridge2 = facet->ridges.values[iridge2];
          if (ridge2 != horizonridge) {
            array_append(&neighbour->ridges, ridge2, polytoop->allocator);
            if (ridge2->facets[0] == facet) {
              ridge2->facets[0] = neighbour;
            } else {
              ridge2->facets[1] = neighbour;
            }
          }
        }

        /* Add apex vertex to neighbour: */
        array_append(&neighbour->vertices, vertex, polytoop->allocator);

        /* Add neighbour centroid and volume: */
        newvol = facet->volume + neighbour->volume;
        vec_scale(polytoop->dim, neighbour->centroid,
                  neighbour->volume / newvol);
        vec_adds(polytoop->dim, neighbour->centroid, facet->centroid,
                 facet->volume / newvol);
        neighbour->volume = newvol;

        /* Disassociate the dividing ridge from the neighbour: */
        array_remove(&neighbour->ridges,
                     array_find(neighbour->ridges, horizonridge));

        /* Disassociate vertices from the dividing ridge: */
        for (ivertex = 0; ivertex < polytoop->dim - 1; ++ivertex) {
          --horizonridge->vertices[ivertex]->nridges;
          assert(horizonridge->vertices[ivertex]->nridges > 0);
        }

        /* Remove the dividing ridge: */
        ridge_remove(polytoop, horizonridge);
        array_remove(&horizonridges, array_find(horizonridges, horizonridge));
        --iridge;

        /* Remove the facet and replace by neighbour: */
        facet_remove(polytoop, facet);
        facet = neighbour; /* so it is added to newfacets */
      }
    }

    /* Add new facet to newfacets array: */
    array_append(&newfacets, facet, polytoop->allocator);
  }

  /* Assign outside verts: */
  point = outsidepoints;
  while (point) {
    Point* next;

    /* Remember next in list: */
    next = point->next;

    /* Find visible facet: */
    for (ifacet = 0; ifacet < newfacets.len; ++ifacet) {
      facet = newfacets.values[ifacet];
      memcpy(vec, point->pos, polytoop->dim * sizeof(double));
      vec_sub(polytoop->dim, vec, facet->centroid);
      ip = vec_dot(polytoop->dim, vec, facet->normal);
      if (ip > EPS) {
        /* Outside. */
        point->height = ip;
        facet_addoutside(polytoop, facet, point);
        break;
      }
    }

    /* Next in list: */
    point = next;
  }

  array_delete(&horizonridges, polytoop->allocator);
  array_delete(&newfacets, polytoop->allocator);
}



static void build(Polytoop* polytoop, int npoints, Point** points)
{
  assert(npoints >= polytoop->dim + 1);

  initialsimplex(polytoop, npoints, points);

  /* Now, all vertices are either on the initial simplex OR
     added to the outside set of one of its facets OR
     inside the initial simplex. */
  while (polytoop->firstfacet->outsidehead) {
    Point* point;

    /* Pop head: */
    point = polytoop->firstfacet->outsidehead;
    polytoop->firstfacet->outsidehead = point->next;

    /* Add vertex to polytoop: */
    addpoint(polytoop, polytoop->firstfacet, point);
  }
}



Polytoop* polytoop_new()
{
  Polytoop* polytoop;

  polytoop = malloc(sizeof(Polytoop));

  polytoop->allocator = allocator_new();

  /* Member variables: */
  polytoop->dim = 0;
  polytoop->isdelaunay = 0;
  polytoop->shift = NULL;
  polytoop->scales = NULL;
  polytoop->center = NULL;

  /* Facet linked list: */
  polytoop->nfacets = 0;
  polytoop->firstfacet = NULL;
  polytoop->lastfacet = NULL;

  /* Ridges linked list: */
  polytoop->nridges = 0;
  polytoop->firstridge = NULL;
  polytoop->lastridge = NULL;

  /* Vertex linked list: */
  polytoop->nverts = 0;
  polytoop->firstvertex = NULL;

  polytoop->merge = 0;

  polytoop->newridges = hashmap_new(hashvertset, compvertsets);
  hashmap_setdata(polytoop->newridges, &polytoop->dim);

  return polytoop;
}



void polytoop_delete(Polytoop* polytoop)
{
  polytoop_clear(polytoop);
  hashmap_delete(polytoop->newridges);
  allocator_delete(polytoop->allocator);
  free(polytoop);
}



void polytoop_setmerge(Polytoop* polytoop, int merge)
{
  if (merge) {
    polytoop->merge = 1;
  } else {
    polytoop->merge = 0;
  }
}



void polytoop_clear(Polytoop* polytoop)
{
#ifndef NDEBUG
  polytoop_Facet* facet;
  Ridge* ridge;
  polytoop_Vertex* vertex;

  /* Deallocate facets: */
  facet = polytoop->firstfacet;
  while (facet) {
    polytoop_Facet* next = facet->next;
    facet_free(polytoop, facet);
    facet = next;
  }

  /* Deallocate ridges: */
  ridge = polytoop->firstridge;
  while (ridge) {
    Ridge* next = ridge->next;
    ridge_free(polytoop, ridge);
    ridge = next;
  }

  /* Deallocate vertices: */
  vertex = polytoop->firstvertex;
  while (vertex) {
    polytoop_Vertex* next = vertex->next;
    vertex_free(polytoop, vertex);
    vertex = next;
  }
  allocator_free(polytoop->allocator, polytoop->center,
                 polytoop->dim * sizeof(double));
  allocator_free(polytoop->allocator, polytoop->scales,
                 polytoop->dim * sizeof(double));
  allocator_free(polytoop->allocator, polytoop->shift,
                 polytoop->dim * sizeof(double));
#endif

  hashmap_clear(polytoop->newridges);
  polytoop->nfacets = 0;
  polytoop->firstfacet = NULL;
  polytoop->lastfacet = NULL;
  polytoop->nridges = 0;
  polytoop->firstridge = NULL;
  polytoop->lastridge = NULL;
  polytoop->nverts = 0;
  polytoop->firstvertex = NULL;
  polytoop->center = NULL;
  polytoop->scales = NULL;
  polytoop->shift = NULL;
  polytoop->isdelaunay = 0;
  polytoop->dim = 0;
  allocator_clear(polytoop->allocator);
}



void polytoop_fromplanes(Polytoop* polytoop, int n, int d, double* normals,
                         double* dists)
{
  int i;
  int ifacet;
  int ivert;
  int nfacets;
  int* p;
  double dist;
  double* mata;
  double* b;
  double* c;
  double* x;
  double* points;
  double* verts;
  double* dcmp;
  double* matq;
  double* centroid;
  double* normal;
  Polytoop* recpolytoop;
  polytoop_Facet* facet;

  /* Default result: */
  polytoop_clear(polytoop);

  /* Allocations: */
  mata = malloc(n * (d + 1) * sizeof(double));
  b = malloc(n * sizeof(double));
  c = malloc((d + 1) * sizeof(double));
  x = malloc((d + 1) * sizeof(double));

  /* Copy data to augmented matrix: */
  for (i = 0; i < n; ++i) {
    memcpy(&mata[i * (d + 1)], &normals[i * d], d * sizeof(double));
    mata[i * (d + 1) + d] = 1.0;
    b[i] = dists[i];
  }

  /* Gradient: */
  vec_reset(d, c);
  c[d] = 1.0;

  /* Feasible start point: */
  vec_reset(d, x);
  x[d] = b[vec_minindex(n, b)] - 1.0;

  /* Find interior point. */
  linprog_cn(0, n, d + 1, mata, b, c, x);
  if (x[d] <= 0.0) {
    /* Infeasible. */
    goto earlycleanup;
  }

  /* Construct reciprocal points: */
  points = malloc(n * d * sizeof(double));
  for (i = 0; i < n; ++i) {
    memcpy(&points[i * d], &normals[i * d], d * sizeof(double));
    dist = dists[i] - vec_dot(d, x, &normals[i * d]);
    assert(dist > 0.0);
    vec_scale(d, &points[i * d], 1.0 / dist);
  }

  /* Construct reciprocal polytoop: */
  recpolytoop = polytoop_new();
  polytoop_frompoints(recpolytoop, n, d, points);
  free(points);

  /* Some allocations for backtransforming from reciprocal space: */
  verts = malloc(d * d * sizeof(double));
  dcmp = malloc((d - 1) * d * sizeof(double));
  matq = malloc(d * d * sizeof(double));
  centroid = malloc(d * sizeof(double));
  normal = malloc(d * sizeof(double));
  p = malloc((d - 1) * sizeof(int));

  /* Convert reciprocal polytoop to polytoop: */
  nfacets = polytoop_getnumfacets(recpolytoop);
  points = malloc(nfacets * d * sizeof(double));
  for (facet = recpolytoop->firstfacet, ifacet = 0; facet != NULL;
       facet = facet->next, ++ifacet) {
    /* Get vertices and centroid: */
    vec_reset(d, centroid);
    for (ivert = 0; ivert < d; ++ivert) {
      polytoop_vertex_getposition(polytoop_facet_getvertex(facet, ivert),
                                  &verts[ivert * d]);
      vec_add(d, centroid, &verts[ivert * d]);
    }
    vec_scale(d, centroid, 1.0 / (double)d);

    /* Vertices relative to first vertex: */
    for (ivert = 0; ivert < d - 1; ++ivert) {
      memcpy(&dcmp[ivert * d], &verts[(ivert + 1) * d], d * sizeof(double));
      vec_sub(d, &dcmp[ivert * d], &verts[0 * d]);
    }

    /* LQ decomposition of vertices: */
    lqdc(d - 1, d, dcmp, p);

    /* Q-matrix:*/
    lqformq(d - 1, d, dcmp, matq);

    /* Last vector of Q is normal: */
    memcpy(normal, &matq[(d - 1) * d], d * sizeof(double));

    /* Plane distance from origin: */
    dist = vec_dot(d, normal, centroid);
    if (dist < 0.0) {
      vec_neg(d, normal);
      dist *= -1.0;
    }

    /* Construct vertex position: */
    memcpy(&points[ifacet * d], normal, d * sizeof(double));
    vec_scale(d, &points[ifacet * d], 1.0 / dist);
    vec_add(d, &points[ifacet * d], x);
  }

  polytoop_frompoints(polytoop, nfacets, d, points);

  /* Clean up: */
  free(p);
  free(normal);
  free(centroid);
  free(matq);
  free(dcmp);
  free(verts);
  polytoop_delete(recpolytoop);
  free(points);
earlycleanup:
  free(x);
  free(c);
  free(b);
  free(mata);
}



void polytoop_frompoints(Polytoop* polytoop, int npoints, int dim,
                         double* orgpoints)
{
  int idim;
  int ipoint;
  int* minindices;
  int* maxindices;
  double* minima;
  double* maxima;
  Point** points;

  /* Clear polytoop: */
  polytoop_clear(polytoop);

  /* Initialize member variables: */
  polytoop->dim = dim;
  polytoop->isdelaunay = 0;
  polytoop->shift =
      allocator_alloc(polytoop->allocator, polytoop->dim * sizeof(double));
  polytoop->scales =
      allocator_alloc(polytoop->allocator, polytoop->dim * sizeof(double));
  polytoop->center =
      allocator_alloc(polytoop->allocator, polytoop->dim * sizeof(double));

  /* Bounding box: */
  minindices = alloca(polytoop->dim * sizeof(int));
  maxindices = alloca(polytoop->dim * sizeof(int));
  minima = alloca(polytoop->dim * sizeof(double));
  maxima = alloca(polytoop->dim * sizeof(double));
  boundingbox(npoints, polytoop->dim, orgpoints, minindices, maxindices, minima,
              maxima);

  /* Scales: */
  memcpy(polytoop->scales, maxima, polytoop->dim * sizeof(double));
  vec_sub(polytoop->dim, polytoop->scales, minima);

  /* Translation: */
  memcpy(polytoop->shift, minima, polytoop->dim * sizeof(double));

  /* Points array: */
  points = malloc(npoints * (sizeof(Point*)));
  for (ipoint = 0; ipoint < npoints; ++ipoint) {
    points[ipoint] = malloc(sizeof(Point));
    points[ipoint]->next = NULL;
    points[ipoint]->d = polytoop->dim;
    points[ipoint]->index = ipoint;
    points[ipoint]->height = 0.0;
    points[ipoint]->pos = malloc(polytoop->dim * sizeof(double));
    for (idim = 0; idim < polytoop->dim; ++idim) {
      points[ipoint]->pos[idim] =
          (orgpoints[ipoint * polytoop->dim + idim] - polytoop->shift[idim]) /
          polytoop->scales[idim];
    }
  }

  /* Build the polytoop: */
  build(polytoop, npoints, points);

  /* Clean up: */
  for (ipoint = npoints - 1; ipoint >= 0; --ipoint) {
    free(points[ipoint]->pos);
    free(points[ipoint]);
  }
  free(points);
}



void polytoop_delaunay(Polytoop* polytoop, int npoints, int dim,
                       double* orgpoints)
{
  int i;
  int j;
  int k;
  int ipoint;
  int* minindices;
  int* maxindices;
  double recd;
  Ridge* ridge;
  Ridge** ridges;
  polytoop_Facet* facet;
  Point** points;

  /* Reciprocal of dimension: */
  recd = 1.0 / (double)dim;

  /* Clear polytoop: */
  polytoop_clear(polytoop);

  /* Initialize member variables: */
  polytoop->dim = dim + 1;
  polytoop->isdelaunay = 1;
  polytoop->shift =
      allocator_alloc(polytoop->allocator, polytoop->dim * sizeof(double));
  polytoop->scales =
      allocator_alloc(polytoop->allocator, polytoop->dim * sizeof(double));
  polytoop->center =
      allocator_alloc(polytoop->allocator, polytoop->dim * sizeof(double));

  /* Bounding box: */
  minindices = alloca(dim * sizeof(int));
  maxindices = alloca(dim * sizeof(int));
  boundingbox(npoints, dim, orgpoints, minindices, maxindices, polytoop->shift,
              polytoop->scales);
  polytoop->shift[dim] = 0.0;

  /* Scales: */
  vec_sub(dim, polytoop->scales, polytoop->shift);
  polytoop->scales[dim] = 1.0;

  /* Points array: */
  points = malloc((npoints + 1) * sizeof(Point*));
  for (ipoint = 0; ipoint < npoints + 1; ++ipoint) {
    points[ipoint] = allocator_alloc(polytoop->allocator, sizeof(Point));
    points[ipoint]->next = NULL;
    points[ipoint]->d = polytoop->dim;
    points[ipoint]->index = ipoint;
    points[ipoint]->height = 0.0;
    points[ipoint]->pos =
        allocator_alloc(polytoop->allocator, polytoop->dim * sizeof(double));

    if (ipoint < npoints) {
      /* Transformed point: */
      for (i = 0; i < dim; ++i) {
        points[ipoint]->pos[i] =
            (orgpoints[ipoint * dim + i] - polytoop->shift[i]) /
            polytoop->scales[i];
      }

      /* Create paraboloid in extra dimension: */
      points[ipoint]->pos[dim] = recd * vec_nrmsq(dim, points[ipoint]->pos);
    }
  }

  /* Add point above paraboloid to guarantee full dimensionality: */
  vec_reset(dim, points[npoints]->pos);
  for (ipoint = 0; ipoint < npoints; ++ipoint) {
    vec_add(dim, points[npoints]->pos, points[ipoint]->pos);
  }
  vec_scale(dim, points[npoints]->pos, 1.0 / (double)npoints);
  points[npoints]->pos[dim] = 2.0;

  /* Construct convex polytoop of paraboloid: */
  build(polytoop, npoints + 1, points);

  /* Clean up: */
  for (ipoint = npoints; ipoint >= 0; --ipoint) {
    allocator_free(polytoop->allocator, points[ipoint]->pos,
                   polytoop->dim * sizeof(double));
    allocator_free(polytoop->allocator, points[ipoint], sizeof(Point));
  }
  free(points);

  /* Remove upper delaunay surfaces: */
  ridges = malloc(polytoop->nridges * sizeof(Ridge*));
  ridge = polytoop->firstridge;
  for (i = 0, ridge = polytoop->firstridge; i < polytoop->nridges;
       ++i, ridge = ridge->next) {
    ridges[i] = ridge;
  }

  for (i = 0; i < polytoop->nridges; ++i) {
    ridge = ridges[i];
    for (j = 0; j < 2; ++j) {
      facet = ridge->facets[j];

      if (facet != NULL && facet->normal[polytoop->dim - 1] > 0.0) {
        /* Upper delaunay. */

        /* Disassociate facet from ridges: */
        ridge->facets[j] = NULL;

        for (k = 0; k < polytoop->dim; ++k) {
          Ridge* ridgek;

          ridgek = facet->ridges.values[k];
          if (ridgek->facets[0] == facet) {
            ridgek->facets[0] = NULL;
          } else if (ridgek->facets[1] == facet) {
            ridgek->facets[1] = NULL;
          }
        }

        /* Remove facet: */
        facet_remove(polytoop, facet);
      }
    }

    if (ridge->facets[0] == NULL && ridge->facets[1] == NULL) {
      /* Both adjacent facets are gone. */

      /* Disassociate from adjacent vertices: */
      for (j = 0; j < polytoop->dim - 1; ++j) {
        --ridge->vertices[j]->nridges;
        if (ridge->vertices[j]->nridges == 0) {
          /* Remove vertex as well: */
          vertex_remove(polytoop, ridge->vertices[j]);
        }
      }

      /* Remove ridge: */
      ridge_remove(polytoop, ridge);
      ridges[i] = NULL; /* forget reference */
    }
  }
  free(ridges);
}



void polytoop_addvertex(Polytoop* polytoop, double* point)
{
  int idim;
  double dist;
  double* pos;
  double* vec;
  Point apex;
  polytoop_Facet* maxfacet;
  polytoop_Facet* facet;

  /* Allocations: */
  pos = alloca(polytoop->dim * sizeof(double));
  vec = alloca(polytoop->dim * sizeof(double));

  /* Initialize apex: */
  apex.next = NULL;
  apex.d = polytoop->dim;
  apex.index = polytoop->nverts;
  apex.height = 0.0;
  apex.pos = pos;

  /* Transformed position: */
  memcpy(pos, point, polytoop->dim * sizeof(double));
  vec_sub(polytoop->dim, pos, polytoop->shift);
  for (idim = 0; idim < polytoop->dim; ++idim) {
    if (polytoop->scales[idim] > 0.0) {
      pos[idim] /= polytoop->scales[idim];
    }
  }

  /* Find facet that can see the point: */
  maxfacet = NULL;
  for (facet = polytoop->firstfacet; facet != NULL; facet = facet->next) {
    memcpy(vec, pos, polytoop->dim * sizeof(double));
    vec_sub(polytoop->dim, vec, facet->centroid);
    dist = vec_dot(polytoop->dim, vec, facet->normal);
    if (dist > EPS) {
      maxfacet = facet;
      break;
    }
  }

  if (maxfacet != NULL) {
    /* Add vertex to polytoop: */
    addpoint(polytoop, maxfacet, &apex);
  }
}



void polytoop_print(Polytoop* polytoop)
{
  int i;
  int j;
  int d;
  double* sv;
  polytoop_Vertex* vertex;
  polytoop_Facet* facet;

  printf("%d facets\n", polytoop->nfacets);
  printf("%d ridges\n", polytoop->nridges);
  printf("%d vertices\n", polytoop->nverts);

  d = polytoop->dim;
  sv = alloca(d * sizeof(double));
  vec_reset(d, sv);
  for (facet = polytoop->firstfacet, i = 0; facet != NULL;
       facet = facet->next, ++i) {

    vec_adds(d, sv, facet->normal, facet->volume);

    printf("facet %d\n", i + 1);
    printf("  volume: %g\n", facet->volume);
    printf("  centroid: ");
    vec_print(polytoop->dim, facet->centroid);
    printf("\n");
    printf("  normal: ");
    vec_print(polytoop->dim, facet->normal);
    printf("\n");

    printf("  vertices:\n");
    for (j = 0; j < facet->vertices.len; ++j) {
      vertex = facet->vertices.values[j];
      vec_print(polytoop->dim, vertex->position);
    }
    printf("\n");
  }
  assert(i == polytoop->nfacets);
  printf("  sv:\n");
  vec_print(d, sv);
  printf("\n");
}



int polytoop_getnumfacets(Polytoop* polytoop) { return polytoop->nfacets; }



int polytoop_getnumridges(Polytoop* polytoop) { return polytoop->nridges; }



int polytoop_getnumvertices(Polytoop* polytoop) { return polytoop->nverts; }



int polytoop_getvertexindex(Polytoop* polytoop, polytoop_Facet* facet,
                            int ivertex)
{
  polytoop_Vertex* vertex;

  (void)polytoop;
  vertex = facet->vertices.values[ivertex];
  return vertex->index;
}



void polytoop_interpolate(Polytoop* polytoop, double const* xi, int* indices,
                          double* weights)
{
  int i;
  int maxindex;
  int iridge;
  int ivertex;
  double h;
  double volume;
  double hmax;
  double totalweight;
  double* verts;
  double* matq;
  double* normal;
  double* vec;
  double* centroid;
  double* xiprime;
  polytoop_Facet* currentfacet;
  Ridge* ridge;

  assert(polytoop->isdelaunay);

  /* Some space for various tasks: */
  int dim = polytoop->dim - 1;
  xiprime = alloca(dim * sizeof(double));
  verts = alloca(dim * dim * sizeof(double));
  matq = alloca(dim * dim * sizeof(double));
  normal = alloca(dim * sizeof(double));
  vec = alloca(dim * sizeof(double));
  centroid = alloca(dim * sizeof(double));

  /* Transform xi to local coordinates: */
  for (i = 0; i < dim; ++i) {
    xiprime[i] = (xi[i] - polytoop->shift[i]) / polytoop->scales[i];
  }

  /* Start with first facet: */
  assert(polytoop->firstfacet != NULL);
  currentfacet = polytoop->firstfacet;

  while (1) {
    /* In the current facet, find the ridge where xi is highest above: */
    hmax = -HUGE_VAL;
    maxindex = -1;
    totalweight = 0.0;
    for (iridge = 0; iridge < polytoop->dim; ++iridge) {
      /* Retrieve ridge: */
      ridge = currentfacet->ridges.values[iridge];

      /* On demand construction of ridge normal and centroid (d - 1): */
      if (ridge->normal == NULL) {
        /* Matrix of vertex coordinates: */
        for (ivertex = 0; ivertex < dim; ++ivertex) {
          memcpy(&verts[ivertex * dim],
                 ridge->vertices[ivertex]->position,
                 dim * sizeof(double));
        }

        /* Analyse ridge simplex: */
        ridge->centroid =
            allocator_alloc(polytoop->allocator, dim * sizeof(double));
        analysesimplex(dim, dim, verts,
                       &ridge->volume, ridge->centroid);

        /* Construct normal: */
        ridge->normal =
            allocator_alloc(polytoop->allocator, dim * sizeof(double));
        memcpy(ridge->normal, ridge->centroid, dim * sizeof(double));
        vec_sub(dim, ridge->normal, currentfacet->centroid);
        for (i = 0; i < dim - 1; ++i) {
          double fac = vec_dot(dim, ridge->normal, &verts[i * dim]);
          vec_adds(dim, ridge->normal, &verts[i * dim], -fac);
        }
        vec_normalize(dim, ridge->normal);
      }

      /* Retrieve normal, centroid, volume: */
      memcpy(normal, ridge->normal, dim * sizeof(double));
      memcpy(centroid, ridge->centroid, dim * sizeof(double));
      volume = ridge->volume;

      /* Facet centroid to ridge centroid: */
      memcpy(vec, centroid, dim * sizeof(double));
      vec_sub(dim, vec, currentfacet->centroid);

      /* Normal wants to be outward pointing: */
      if (vec_dot(dim, vec, normal) < 0.0) {
        vec_neg(dim, normal);
      }

      /* Ridge centroid to interpolation point: */
      memcpy(vec, xiprime, dim * sizeof(double));
      vec_sub(dim, vec, centroid);

      /* Height of interpolation point above ridge: */
      h = vec_dot(dim, vec, normal);

      /* Weight for this ridge: */
      weights[iridge] = -h * volume;
      totalweight += weights[iridge];

      /* Index for this ridge: */
      indices[iridge] =
          ((polytoop_Vertex*)
               currentfacet->vertices.values[dim - iridge])
              ->index;

      /* Keep track of maximum height (if not boundary ridge): */
      if (h > hmax && ridge->facets[0] != NULL && ridge->facets[1] != NULL) {
        hmax = h;
        maxindex = iridge;
      }
    }

    /* If hmax insignificant, stop: */
    if (hmax < EPS) {
      break;
    }

    /* Else, goto next facet (cross ridge): */
    ridge = currentfacet->ridges.values[maxindex];
    if (ridge->facets[0] == currentfacet) {
      currentfacet = ridge->facets[1];
    } else {
      assert(ridge->facets[1] == currentfacet);
      currentfacet = ridge->facets[0];
    }
  }

  if (totalweight > 0.0) {
    vec_scale(polytoop->dim, weights, 1.0 / totalweight);
  }
}



polytoop_Facet* polytoop_firstfacet(Polytoop* polytoop)
{
  return polytoop->firstfacet;
}



polytoop_Vertex* polytoop_firstvertex(Polytoop* polytoop)
{
  return polytoop->firstvertex;
}



polytoop_Facet* polytoop_facet_nextfacet(polytoop_Facet* facet)
{
  return facet->next;
}



int polytoop_facet_getnumneighbours(polytoop_Facet* facet)
{
  return facet->ridges.len;
}



polytoop_Facet* polytoop_facet_getneighbour(polytoop_Facet* facet, int i)
{
  Ridge* ridge;
  assert(i < facet->ridges.len);
  ridge = facet->ridges.values[i];
  if (ridge->facets[0] == facet) {
    return ridge->facets[1];
  }
  assert(ridge->facets[1] == facet);
  return ridge->facets[0];
}



void polytoop_facet_getnormal(polytoop_Facet* facet, double* normal)
{
  int i;

  /* Scale normal: */
  memcpy(normal, facet->normal, facet->polytoop->dim * sizeof(double));
  for (i = 0; i < facet->polytoop->dim; ++i) {
    normal[i] /= facet->polytoop->scales[i];
  }

  /* Normalize: */
  vec_normalize(facet->polytoop->dim, normal);
}



void polytoop_facet_getcentroid(polytoop_Facet* facet, double* centroid)
{
  int i;

  memcpy(centroid, facet->centroid, facet->polytoop->dim * sizeof(double));
  for (i = 0; i < facet->polytoop->dim; ++i) {
    centroid[i] *= facet->polytoop->scales[i];
  }
  vec_add(facet->polytoop->dim, centroid, facet->polytoop->shift);
}



double polytoop_facet_getoffset(polytoop_Facet* facet)
{
  double* centroid;
  double* normal;

  centroid = alloca(facet->polytoop->dim * sizeof(double));
  normal = alloca(facet->polytoop->dim * sizeof(double));
  polytoop_facet_getcentroid(facet, centroid);
  polytoop_facet_getnormal(facet, normal);
  return vec_dot(facet->polytoop->dim, centroid, normal);
}



double polytoop_facet_getvolume(polytoop_Facet* facet) { return facet->volume; }



int polytoop_facet_getnumvertices(polytoop_Facet* facet)
{
  return facet->polytoop->dim;
}



polytoop_Vertex* polytoop_facet_getvertex(polytoop_Facet* facet, int ivertex)
{
  return facet->vertices.values[ivertex];
}



polytoop_Vertex* polytoop_vertex_nextvertex(polytoop_Vertex* vertex)
{
  return vertex->next;
}



void polytoop_vertex_getposition(polytoop_Vertex* vertex, double* position)
{
  int i;
  int dim;

  dim = vertex->polytoop->dim;
  if (vertex->polytoop->isdelaunay) {
    --dim;
  }

  memcpy(position, vertex->position, dim * sizeof(double));
  for (i = 0; i < dim; ++i) {
    position[i] *= vertex->polytoop->scales[i];
  }
  vec_add(dim, position, vertex->polytoop->shift);
}
