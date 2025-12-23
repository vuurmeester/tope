#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "hashmap.h"
#include "math.h"
#include "util.h"

#include <tope.h>

#define EPS 1.0e-9



/* Create and initialize a vertex struct with a position and id: */
static Vertex* vertex_create(Tope* tope, double const* pos, int index)
{
  ++tope->nverts;
  Allocator* alc = &tope->alc;

  /* Allocate vertex: */
  Vertex* vertex = allocator_alloc(
    alc,
    sizeof(Vertex) + (tope->dim - 1) * sizeof(double)
  );

  /* Fill in other variables: */
  vertex->index = index;
  vertex->nridges = 0;
  memcpy(vertex->position, pos, tope->dim * sizeof(double));

  return vertex;
}



static void vertex_remove(Tope* tope, Vertex* vertex)
{
  --tope->nverts;
  allocator_free(
    &tope->alc,
    vertex,
    sizeof(Vertex) + (tope->dim - 1) * sizeof(double)
  );
}



/* Create and initialize a ridge struct with vertices: */
static Ridge* ridge_create(Tope* tope, Vertex** vertices)
{
  Allocator* alc = &tope->alc;
  int d = tope->dim;

  /* Create ridge: */
  Ridge* ridge = allocator_alloc(alc, sizeof(Ridge) + (d - 2) * sizeof(Vertex*));
  ridge->vdn = NULL;
  ridge->facets[0] = NULL;
  ridge->facets[1] = NULL;
  memcpy(ridge->vertices, vertices, (d - 1) * sizeof(Vertex*));

  /* Associate vertices with ridge: */
  for (int i = 0; i < d - 1; ++i) {
    ++vertices[i]->nridges;
  }

  /* Append ridge to tope list: */
  ++tope->nridges;

  return ridge;
}



static void ridge_remove(Tope* tope, Ridge* ridge)
{
  Allocator* alc = &tope->alc;
  int d = tope->dim;
  --tope->nridges;

  /* Remove reference to this ridge from adjacent vertices: */
  for (int ivertex = 0; ivertex < tope->dim - 1; ++ivertex) {
    Vertex* vertex = ridge->vertices[ivertex];
    --vertex->nridges;
    if (vertex->nridges == 0) {
      /* Remove the vertex as well: */
      vertex_remove(tope, vertex);
    }
  }

  if (ridge->vdn != NULL) {
    assert(tope->isdelaunay);
    allocator_free(alc, ridge->vdn, (d + 1) * sizeof(double));
  }
  allocator_free(
    alc, 
    ridge,
    sizeof(Ridge) + (d - 2) * sizeof(Vertex*)
  );
}



/* Create and initialize a facet struct: */
static Facet* facet_create(Tope* tope)
{
  ++tope->nfacets;
 
  Allocator* alc = &tope->alc;
  int d = tope->dim;

  Facet* facet    = allocator_alloc(alc,     sizeof(Facet  ));
  facet->centroid = allocator_alloc(alc, d * sizeof(double ));
  facet->normal   = allocator_alloc(alc, d * sizeof(double ));
  facet->ridges   = allocator_alloc(alc, d * sizeof(Ridge* ));
  facet->verts = NULL;

  /* Append to list: */
  facet->prev = tope->lastfacet;
  if (facet->prev != NULL) {
    facet->prev->next = facet;
  }
  else {
    tope->firstfacet = facet;
  }
  tope->lastfacet = facet;
  facet->next = NULL;

  facet->volume = 0.0;
  facet->dist = 0.0;
  facet->outsidehead = NULL;
  facet->outsidetail = NULL;
  facet->visible = false;

  vec_reset(d, facet->centroid);
  vec_reset(d, facet->normal);
  memset(facet->ridges, 0x00, d * sizeof(Ridge* ));

  return facet;
}



static void facet_free(Tope* tope, Facet* facet)
{
  int d = tope->dim;
  Allocator* alc = &tope->alc;
  while (facet->verts) {
    List* next = facet->verts->tail;
    allocator_free(alc, facet->verts, sizeof(List));
    facet->verts = next;
  }
  allocator_free(alc, facet->ridges  , d * sizeof(Ridge*));
  allocator_free(alc, facet->normal  , d * sizeof(double));
  allocator_free(alc, facet->centroid, d * sizeof(double));
  allocator_free(alc, facet, sizeof(Facet));
}



static void facet_remove(Tope* tope, Facet* facet)
{
  if (facet->next != NULL) {
    facet->next->prev = facet->prev;
  }
  else {
    assert(tope->lastfacet == facet);
    tope->lastfacet = facet->prev;
  }

  if (facet->prev != NULL) {
    facet->prev->next = facet->next;
  }
  else {
    assert(tope->firstfacet == facet);
    tope->firstfacet = facet->next;
  }

  facet->next = NULL;
  facet->prev = NULL;
  --tope->nfacets;
}



static void facet_addoutside(Tope* tope, Facet* facet, Point* point)
{
  if (!facet->outsidehead) {
    /* Initialize list: */
    facet->outsidehead = point;
    facet->outsidetail = point;
    point->next = NULL;
  }
  else if (point->height > facet->outsidehead->height) {
    /* Prepend: */
    point->next = facet->outsidehead;
    facet->outsidehead = point;
  }
  else {
    /* Append: */
    facet->outsidetail->next = point;
    facet->outsidetail = point;
    point->next = NULL;
  }

  /* Move facet to front: */
  if (facet->prev != NULL) {
    /* Unlink from facet list: */
    if (facet->next != NULL) {
      facet->next->prev = facet->prev;
    }
    else {
      assert(tope->lastfacet == facet);
      tope->lastfacet = facet->prev;
    }
    facet->prev->next = facet->next;

    Facet* firstfacet = tope->firstfacet;
    if (firstfacet->outsidehead != NULL &&
        firstfacet->outsidehead->height > facet->outsidehead->height) {
      /* Move to second: */
      facet->next = firstfacet->next;
      if (facet->next != NULL) {
        facet->next->prev = facet;
      }
      facet->prev = tope->firstfacet;
      firstfacet->next = facet;
    }
    else {
      /* Move to front: */
      facet->next = tope->firstfacet;
      firstfacet->prev = facet;
      tope->firstfacet = facet;
      facet->prev = NULL;
    }
  }
}



/* Integrate facet into tope (compute volume, create ridges, etc): */
static void addfacet(Tope* tope, Facet* facet, Ridge** ridges)
{
  int d = tope->dim;

  /* Facet centroid: */
  vec_reset(d, facet->centroid);
  for (List* lst = facet->verts; lst != NULL; lst = lst->tail) {
    Vertex* vertex = lst->head;
    vec_add(d, facet->centroid, vertex->position);
  }
  vec_scale(d, facet->centroid, 1.0 / (double)d);

  /* Vertex span: */
  double* dirs = alloca((d - 1) * d * sizeof(double));
  List* lst = facet->verts;
  Vertex* vertex = lst->head;
  double* x0 = vertex->position;
  int i = 0;
  lst = lst->tail;
  for (int i = 0; lst != NULL; ++i) {
    /* Vertex position difference: */
    vertex = lst->head;
    memcpy(dirs + i * d, vertex->position, d * sizeof(double));
    vec_sub(d, dirs + i * d, x0);
    lst = lst->tail;
  }

  /* Basis and volume: */
  facet->volume = 1.0;
  for (int i = 0; i < d - 1; ++i) {
    /* Pivot largest row on top: */
    int pivot = i;
    double maxnrmsq = vec_nrmsq(d, dirs + i * d);
    for (int j = i + 1; j < d - 1; ++j) {
      double nrmsq = vec_nrmsq(d, dirs + j * d);
      if (nrmsq > maxnrmsq) {
        maxnrmsq = nrmsq;
        pivot = j;
      }
    }
    if (pivot > i) {
      memswp(dirs + i * d, dirs + pivot * d, d * sizeof(double));
    }

    /* Normalize: */
    double nrm = sqrt(maxnrmsq);
    vec_scale(d, dirs + i * d, 1.0 / nrm);

    /* Orthogonalize: */
    for (int j = i + 1; j < d - 1; ++j) {
      double ip = vec_dot(d, dirs + j * d, dirs + i * d);
      vec_adds(d, dirs + j * d, dirs + i * d, -ip);
    }

    /* Accumulate area: */
    facet->volume *= nrm / (double)(i + 1);
  }

  /* Compute normal: */
  memcpy(facet->normal, facet->centroid, d * sizeof(double));
  vec_sub(d, facet->normal, tope->center);
  for (int i = 0; i < d - 1; ++i) {
    double ip = vec_dot(d, facet->normal, dirs + i * d);
    vec_adds(d, facet->normal, dirs + i * d, -ip);
  }
  vec_normalize(d, facet->normal);

  /* Plane distance: */
  facet->dist = vec_dot(d, facet->normal, facet->centroid);

  /* Assign ridges: */
  memcpy(facet->ridges, ridges, d * sizeof(Ridge*));
  for (int i = 0; i < d; ++i) {
    /* Get ridge: */
    Ridge* ridge = facet->ridges[i];

    /* Remember facet. */
    if (ridge->facets[0] == NULL) {
      ridge->facets[0] = facet;
    }
    else {
      assert(ridge->facets[1] == NULL);
      ridge->facets[1] = facet;
    }
  }
}



/* Create initial simplex, add outside vertices, etc: */
static void initialsimplex(Tope* tope, int npoints, Point* points)
{
  int d = tope->dim;
  Allocator* alc = &tope->alc;

  /* Allocation: */
  double* span = malloc((npoints - 1) * d * sizeof(double));

  /* Point furthest from unit box center: */
  int maxindex = -1;
  double maxdist = 0.0;
  for (int i = 0; i < npoints; ++i) {
    double dist = 0.0;
    for (int j = 0; j < d; ++j) {
      double coord = points[i].pos[j] - 0.5;
      dist += coord * coord;
    }
    if (dist > maxdist) {
      maxindex = i;
      maxdist = dist;
    }
  }

  /* Place extreme point up front: */
  if (maxindex > 0) {
    memswp(&points[0], &points[maxindex], sizeof(Point*));
  }

  /* Initialize span (positions w.r.t. first point): */
  for (int i = 1; i < npoints; ++i) {
    memcpy(&span[(i - 1) * d], points[i].pos, d * sizeof(double));
    vec_sub(d, &span[(i - 1) * d], points[0].pos);
  }

  /* Permutation vector: */
  int* p = malloc(npoints * sizeof(int));
  for (int i = 0; i < npoints; ++i) {
    p[i] = i;
  }

  /* Due to pivoting, the permutation vector will represent the point order
     which leads to the greatest volume. */
  for (int i = 0; i < d; ++i) {
    /* Max norm: */
    double maxnormsq = 0.0;
    int pivot = -1;
    for (int j = i; j < npoints - 1; ++j) {
      double normsq = vec_nrmsq(d, &span[j * d]);
      if (normsq > maxnormsq) {
        pivot = j;
        maxnormsq = normsq;
      }
    }

    /* Perform pivot: */
    if (pivot > i) {
      memswp(&p[i + 1], &p[pivot + 1], sizeof(int));
      memswp(&span[i * d], &span[pivot * d], d * sizeof(double));
    }

    /* Orthogonalize: */
    for (int j = i + 1; j < npoints - 1; ++j) {
      double fac = vec_dot(d, span + i * d, span + j * d) / maxnormsq;
      vec_adds(d, span + j * d, span + i * d, -fac);
    }
  }

  /* Initialize tope center to initial simplex centroid: */
  vec_reset(d, tope->center);
  for (int i = 0; i < d + 1; ++i) {
    vec_add(d, tope->center, points[p[i]].pos);
  }
  vec_scale(d, tope->center, 1.0 / (double)(tope->dim + 1));

  /* Create initial simplex vertices: */
  Vertex** vertices = alloca((d + 1) * sizeof(Vertex*));
  for (int i = 0; i < d + 1; ++i) {
    vertices[i] = vertex_create(tope, points[p[i]].pos, points[p[i]].index);
  }

  Ridge** facetridges = alloca(d * sizeof(Ridge*));
  Vertex** ridgeverts = alloca((d - 1) * sizeof(Vertex*));

  hashmap_clear(&tope->newridges);

  /* Create facets: */
  for (int i = 0; i < d + 1; ++i) {
    /* Create facet: */
    Facet* facet = facet_create(tope);

    /* Add all vertices except for one: */
    List** pli = &facet->verts;
    for (int j = 0; j < d + 1; ++j) {
      if (i == j) {
        continue;
      }
      *pli = allocator_alloc(alc, sizeof(List));
      (*pli)->head = vertices[j];
      (*pli)->tail = NULL;
      pli = &(*pli)->tail;
    }

    /* Add ridges: */
    for (int j = 0; j < d; ++j) {
      /* Ridge verts: */
      List* lst = facet->verts;
      for (int k = 0; k < d - 1; ++k) {
        if (k == j) {
          lst = lst->tail;
        }
        ridgeverts[k] = lst->head;
        lst = lst->tail;
      }

      Ridge** pnewridge = hashmap_get(&tope->newridges, d, ridgeverts);
      if (*pnewridge == NULL) {
        /* new ridge apparently */
        *pnewridge = ridge_create(tope, ridgeverts);
      }
      facetridges[j] = *pnewridge;
    }

    /* Vertex array filled, compute and integrate the facet: */
    addfacet(tope, facet, facetridges);
  }

  /* Create outside sets (assign each remaining vertex to a facet which it
   * 'sees'): */
  for (int i = tope->dim + 1; i < npoints; ++i) {
    /* Determine maximum distance: */
    Facet* facet = tope->firstfacet;
    while (facet != NULL) {
      /* Distance to facet: */
      double h = vec_dot(d, facet->normal, points[p[i]].pos) - facet->dist;

      /* If above facet, add it to outside set and move to next point. */
      if (h > EPS) {
        points[p[i]].height = h;
        facet_addoutside(tope, facet, points + p[i]);
        break;
      }

      facet = facet->next;
    }
  }

  /* Clean up: */
  free(p);
  free(span);
}



/* Add point to tope: */
static void addpoint(Tope* tope, Facet* facet, Point* apex)
{
  int d = tope->dim;
  Allocator* alc = &tope->alc;

  facet_remove(tope, facet);
  facet->visible = true;

  /* Add facet to visible list: */
  Facet* visiblelist = facet;

  /* Initialize outside points list: */
  Point* outsidepoints = NULL;

  while (visiblelist != NULL) {
    /* Pop a value: */
    facet = visiblelist;
    visiblelist = facet->next;

    /* Loop over ridges to visit neighbours: */
    for (int i = 0; i < d; ++i) {
      /* Ridge i: */
      Ridge* ridge = facet->ridges[i];

      /* Retrieve neighbour: */
      Facet* neighbour = NULL;
      if (facet == ridge->facets[0]) {
        /* Current facet is 0, neighbour is 1: */
        ridge->facets[0] = NULL; /* forget reference to this facet */
        neighbour = ridge->facets[1];
      }
      else {
        assert(facet == ridge->facets[1]);
        /* Current facet is 1, neighbour is 0: */
        ridge->facets[1] = NULL; /* forget reference to this facet */
        neighbour = ridge->facets[0];
      }

      if (neighbour == NULL) {
        /* This neighbour was already removed, remove ridge. */

        /* Remove the ridge from the tope: */
        ridge_remove(tope, ridge);
      }
      else {
        if (!neighbour->visible) {
          /* Neighbour is untested. */

          /* Height of apex above neighbour: */
          double* normal = neighbour->normal;
          double h = vec_dot(d, apex->pos, normal) - neighbour->dist;

          /* Decide if neighbour is visible: */
          if (h > EPS) {
            /* Remove neighbour from tope list: */
            facet_remove(tope, neighbour);

            /* Prepend neighbour to visible list: */
            neighbour->next = visiblelist;
            visiblelist = neighbour;
            neighbour->visible = true;
          }
          else {
            /* Neighbour not visible. Ridge belongs to horizon: */
            if (tope->horizonridges_len == tope->horizonridges_cap) {
              tope->horizonridges_cap = 3 * tope->horizonridges_cap / 2 + 1;
              tope->horizonridges = realloc(
                tope->horizonridges,
                tope->horizonridges_cap * sizeof(Ridge*)
              );
            }
            tope->horizonridges[tope->horizonridges_len++] = ridge;
          }
        }
      }
    }

    /* Remember outside vertices of investigated facet: */
    if (facet->outsidehead) {
      facet->outsidetail->next = outsidepoints;
      outsidepoints = facet->outsidehead;
    }

    /* Deallocate investigated facet: */
    facet_free(tope, facet);
  }

  /* Add apex vertex to tope: */
  Vertex* vertex = vertex_create(tope, apex->pos, apex->index);
  Vertex** ridgeverts = alloca((d - 1) * sizeof(Vertex*));
  Ridge** facetridges = alloca(d * sizeof(Ridge*));

  hashmap_clear(&tope->newridges);

  /* Form new facets: */
  while (tope->horizonridges_len > 0) {
    /* New facet: */
    facet = facet_create(tope);

    /* Pop a horizon ridge: */
    Ridge* horizonridge = tope->horizonridges[--tope->horizonridges_len];

    /* Facet vertices are apex + horizon vertices: */
    List** pli = &facet->verts;
    *pli = allocator_alloc(alc, sizeof(List));
    (*pli)->head = vertex;
    (*pli)->tail = NULL;
    pli = &(*pli)->tail;
    for (int ivertex = 1; ivertex < d; ++ivertex) {
      *pli = allocator_alloc(alc, sizeof(List));
      (*pli)->head = horizonridge->vertices[ivertex - 1];
      (*pli)->tail = NULL;
      pli = &(*pli)->tail;
    }

    facetridges[0] = horizonridge;
    for (int iridge = 1; iridge < d; ++iridge) {
      /* Add all facet vertices except the one at index 'iridge': */
      List* li = facet->verts;
      for (int ivertex = 0; ivertex < d - 1; ++ivertex) {
        if (ivertex == iridge) {
          li = li->tail;
        }
        ridgeverts[ivertex] = li->head;
        li = li->tail;
      }

      /* Get or create new ridge: */
      Ridge** pnewridge = hashmap_get(&tope->newridges, d, ridgeverts);
      if (*pnewridge == NULL) {
        /* Create new ridge: */
        *pnewridge = ridge_create(tope, ridgeverts);
      }
      facetridges[iridge] = *pnewridge;
    }

    addfacet(tope, facet, facetridges);

    /* Add new facet to newfacets array: */
    if (tope->newfacets_len == tope->newfacets_cap) {
      tope->newfacets_cap = 3 * tope->newfacets_cap / 2 + 1;
      tope->newfacets = realloc(tope->newfacets, tope->newfacets_cap * sizeof(Facet*));
    }
    tope->newfacets[tope->newfacets_len++] = facet;
  }

  /* Assign outside verts: */
  facet = tope->newfacets[0];
  while (outsidepoints) {
    /* Remember next in list: */
    Point* next = outsidepoints->next;

    /* Height of outside point above facet: */
    double h = vec_dot(tope->dim, outsidepoints->pos, facet->normal) - facet->dist;

    /* Visit neighbours, skipping horizon ridge: */
    Ridge** ridges = facet->ridges;
    Facet* prev = NULL;
    for (int iridge = 1; iridge < d; ++iridge) {
      Ridge* ridge = ridges[iridge];
      Facet* neighbour = NULL;
      if (ridge->facets[0] == facet) {
        neighbour = ridge->facets[1];
      }
      else {
        assert(ridge->facets[1] == facet);
        neighbour = ridge->facets[0];
      }
      if (neighbour == prev) {
        continue;
      }

      /* Neighbour height: */
      double nh = vec_dot(d, outsidepoints->pos, neighbour->normal) - neighbour->dist;
      if (nh > h) {
        prev = facet;
        facet = neighbour;
        ridges = facet->ridges;
        h = nh;
        iridge = 0; /* reset iteration (new base facet) */
      }
    }

    if (h > EPS) {
      /* Outside. */
      outsidepoints->height = h;
      facet_addoutside(tope, facet, outsidepoints);
    }

    /* Next in list: */
    outsidepoints = next;
  }

  assert(tope->horizonridges_len == 0);
  tope->newfacets_len = 0;
}



static void build(Tope* tope, int npoints, Point* points)
{
  assert(npoints >= tope->dim + 1);

  initialsimplex(tope, npoints, points);

  /* Now, all vertices are either on the initial simplex OR
     added to the outside set of one of its facets OR
     inside the initial simplex. */
  Facet* firstfacet = tope->firstfacet;
  while (firstfacet->outsidehead) {
    /* Pop head: */
    Point* point = firstfacet->outsidehead;
    firstfacet->outsidehead = point->next;

    /* Add vertex to tope: */
    addpoint(tope, tope->firstfacet, point);
    firstfacet = tope->firstfacet;
  }
}



static Tope* tope_new()
{
  Tope* tope = calloc(1, sizeof(Tope));

  allocator_init(&tope->alc);
  hashmap_init(&tope->newridges);

  return tope;
}



void tope_delete(Tope* tope)
{
  free(tope->newfacets);
  free(tope->horizonridges);
  hashmap_destroy(&tope->newridges);
  allocator_destroy(&tope->alc);
  free(tope->center);
  free(tope->scales);
  free(tope->shift);
  free(tope);
}



Tope* tope_fromplanes(
  int n,
  int d,
  double const* normals,
  double const* dists,
  double const* xc
)
{
  /* Construct reciprocal points: */
  double* points = malloc(n * d * sizeof(double));
  for (int i = 0; i < n; ++i) {
    memcpy(&points[i * d], &normals[i * d], d * sizeof(double));
    double dist = dists[i] - vec_dot(d, xc, &normals[i * d]);
    assert(dist > 0.0);
    vec_scale(d, &points[i * d], 1.0 / dist);
  }

  /* Construct reciprocal tope: */
  Tope* rectope = tope_frompoints(n, d, points);
  free(points);

  /* Convert reciprocal tope to tope: */
  int nfacets = tope_getnumfacets(rectope);
  points = malloc(nfacets * d * sizeof(double));

  u32 ifacet = 0;
  for (Facet* facet = rectope->firstfacet; facet != NULL; facet = facet->next, ++ifacet) {
    /* Unscaled normal and distance: */
    double dist = facet->dist;
    for (int i = 0; i < d; ++i) {
      points[ifacet * d + i] = facet->normal[i] / rectope->scales[i];
      dist += points[ifacet * d + i] * rectope->shift[i];
    }

    /* Straight space point: */
    vec_scale(d, points + ifacet * d, 1.0 / dist);
    vec_add(d, points + ifacet * d, xc);
  }

  Tope* tope = tope_frompoints(nfacets, d, points);

  /* Clean up: */
  tope_delete(rectope);
  free(points);

  return tope;
}



Tope* tope_frompoints(int npoints, int dim, double const* orgpoints)
{
  int ipoint;
  int idim;
  int* minindices;
  int* maxindices;
  double* minima;
  double* maxima;
  double* positions;
  Point* points;

  Tope* tope = tope_new();

  /* Initialize member variables: */
  tope->dim = dim;
  tope->isdelaunay = false;
  tope->shift  = malloc(tope->dim * sizeof(double));
  tope->scales = malloc(tope->dim * sizeof(double));
  tope->center = malloc(tope->dim * sizeof(double));

  /* Bounding box: */
  minindices = alloca(tope->dim * sizeof(int));
  maxindices = alloca(tope->dim * sizeof(int));
  minima     = alloca(tope->dim * sizeof(double));
  maxima     = alloca(tope->dim * sizeof(double));
  boundingbox(
    npoints,
    tope->dim,
    orgpoints,
    minindices,
    maxindices,
    minima,
    maxima
  );

  /* Scales: */
  memcpy(tope->scales, maxima, tope->dim * sizeof(double));
  vec_sub(tope->dim, tope->scales, minima);

  /* Translation: */
  memcpy(tope->shift, minima, tope->dim * sizeof(double));

  /* Points array: */
  points = malloc(npoints * sizeof(Point));
  positions = malloc(npoints * dim * sizeof(double));
  for (ipoint = 0; ipoint < npoints; ++ipoint) {
    points[ipoint].next = NULL;
    points[ipoint].index = ipoint;
    points[ipoint].height = 0.0;
    points[ipoint].pos = positions + ipoint * dim;
    for (idim = 0; idim < tope->dim; ++idim) {
      points[ipoint].pos[idim] =
          (orgpoints[ipoint * tope->dim + idim] - tope->shift[idim]) /
          tope->scales[idim];
    }
  }

  /* Build the tope: */
  build(tope, npoints, points);

  /* Clean up: */
  free(positions);
  free(points);

  return tope;
}



Tope* tope_delaunay(int npoints, int dim, double const* orgpoints)
{
  Tope* tope = tope_new();

  /* Reciprocal of dimension: */
  double recd = 1.0 / (double)dim;

  /* Initialize member variables: */
  tope->dim = dim + 1;
  tope->isdelaunay = true;
  tope->shift  = malloc(tope->dim * sizeof(double));
  tope->scales = malloc(tope->dim * sizeof(double));
  tope->center = malloc(tope->dim * sizeof(double));

  /* Bounding box: */
  int* minindices = alloca(dim * sizeof(int));
  int* maxindices = alloca(dim * sizeof(int));
  boundingbox(
    npoints,
    dim,
    orgpoints,
    minindices,
    maxindices,
    tope->shift,
    tope->scales
  );
  tope->shift[dim] = 0.0;

  /* Scales: */
  vec_sub(dim, tope->scales, tope->shift);
  tope->scales[dim] = 1.0;

  /* Points array: */
  Point* points = malloc((npoints + 1) * sizeof(Point));
  double* positions = malloc((npoints + 1) * (dim + 1) * sizeof(double));
  for (int ipoint = 0; ipoint < npoints + 1; ++ipoint) {
    points[ipoint].next = NULL;
    points[ipoint].index = ipoint;
    points[ipoint].height = 0.0;
    points[ipoint].pos = positions + ipoint * (dim + 1);

    if (ipoint < npoints) {
      /* Transformed point: */
      for (int i = 0; i < dim; ++i) {
        points[ipoint].pos[i] =
          (orgpoints[ipoint * dim + i] - tope->shift[i]) /
          tope->scales[i];
      }

      /* Create paraboloid in extra dimension: */
      points[ipoint].pos[dim] = recd * vec_nrmsq(dim, points[ipoint].pos);
    }
  }

  /* Add point above paraboloid to guarantee full dimensionality: */
  vec_reset(dim, points[npoints].pos);
  for (int ipoint = 0; ipoint < npoints; ++ipoint) {
    vec_add(dim, points[npoints].pos, points[ipoint].pos);
  }
  vec_scale(dim, points[npoints].pos, 1.0 / (double)npoints);
  points[npoints].pos[dim] = 2.0;

  /* Construct convex tope of paraboloid: */
  build(tope, npoints + 1, points);

  /* Clean up: */
  free(positions);
  free(points);

  /* Find upper Delaunay facet: */
  Facet* facet = tope->firstfacet;
  while (facet != NULL) {
    Facet* next = facet->next;
    if (facet->normal[dim] > -EPS) {
      // Upper delaunay
      for (int i = 0; i <= dim; ++i) {
        Ridge* ridge = facet->ridges[i];
        if (ridge->facets[0] == facet) {
          ridge->facets[0] = NULL;
          if (ridge->facets[1] == NULL) {
            ridge_remove(tope, ridge);
          }
        }
        else {
          assert(ridge->facets[1] == facet);
          ridge->facets[1] = NULL;
          if (ridge->facets[0] == NULL) {
            ridge_remove(tope, ridge);
          }
        }
      }
      facet_remove(tope, facet);
      facet_free(tope, facet);
    }
    facet = next;
  }

  return tope;
}



void tope_addvertex(Tope* tope, double const* point)
{
  /* Allocations: */
  double* pos = alloca(tope->dim * sizeof(double));

  /* Initialize apex: */
  Point apex = {
    .next = NULL,
    .index = tope->nverts,
    .height = 0.0,
    .pos = pos,
  };

  /* Transformed position: */
  memcpy(pos, point, tope->dim * sizeof(double));
  vec_sub(tope->dim, pos, tope->shift);
  for (int idim = 0; idim < tope->dim; ++idim) {
    if (tope->scales[idim] > 0.0) {
      pos[idim] /= tope->scales[idim];
    }
  }

  /* Find facet that can see the point: */
  Facet* maxfacet = NULL;
  Facet* facet = tope->firstfacet;
  while (facet != NULL) {
    double h = vec_dot(tope->dim, pos, facet->normal) - facet->dist;
    if (h > EPS) {
      maxfacet = facet;
      break;
    }
    facet = facet->next;
  }

  if (maxfacet != NULL) {
    /* Add vertex to tope: */
    addpoint(tope, maxfacet, &apex);
  }
}



void tope_print(Tope* tope)
{
  int d = tope->dim;

  double* sv       = alloca(d * sizeof(double));
  double* normal   = alloca(d * sizeof(double));
  double* centroid = alloca(d * sizeof(double));
  double* position = alloca(d * sizeof(double));

  printf("%d facets\n", tope->nfacets);
  printf("%d ridges\n", tope->nridges);
  printf("%d vertices\n", tope->nverts);

  vec_reset(d, sv);

  int i = 0;
  Facet* facet = tope->firstfacet;
  while (facet != NULL) {
    tope_facet_getnormal(tope, facet, normal);
    double volume = tope_facet_getvolume(tope, facet);
    vec_adds(d, sv, normal, volume);
    tope_facet_getcentroid(tope, facet, centroid);

    printf("facet %d\n", i + 1);
    printf("  volume: %g\n", volume);
    printf("  centroid: ");
    vec_print(d, centroid);
    printf("\n");
    printf("  normal: ");
    vec_print(d, normal);
    printf("\n");

    printf("  vertices:\n");
    List* li = facet->verts;
    while (li) {
      Vertex* vertex = li->head;
      tope_vertex_getposition(tope, vertex, position);
      vec_print(d, position);
      li = li->tail;
    }
    printf("\n");

    ++i;
    facet = facet->next;
  }
  assert(i == tope->nfacets);
  printf("  sv:\n");
  vec_print(d, sv);
  printf("\n");
}



int tope_getnumfacets(Tope* tope)
{
  return tope->nfacets;
}



int tope_getnumridges(Tope* tope)
{
  return tope->nridges;
}



int tope_getnumvertices(Tope* tope)
{
  return tope->nverts;
}



void tope_interpolate(
  Tope* tope,
  double const* xi,
  int* indices,
  double* weights
)
{
  Allocator* alc = &tope->alc;

  assert(tope->isdelaunay);

  /* Some space for various tasks: */
  int d = tope->dim - 1;
  double* xiprime = alloca(d * sizeof(double));
  double* centroid = alloca(d * sizeof(double));
  double* verts = alloca(d * d * sizeof(double));

  /* Transform xi to local coordinates: */
  for (int i = 0; i < d; ++i) {
    xiprime[i] = (xi[i] - tope->shift[i]) / tope->scales[i];
  }

  /* Start with first facet: */
  assert(tope->firstfacet != NULL);
  Facet* currentfacet = tope->firstfacet;

  double totalweight = 0.0;
  while (1) {
    Ridge** ridges = currentfacet->ridges;

    /* In the current facet, find the ridge where xi is highest above: */
    totalweight = 0.0;
    double hmin = HUGE_VAL;
    Ridge* minridge = NULL;
    List* li = currentfacet->verts;
    for (int iridge = 0; iridge < tope->dim; ++iridge) {
      /* Retrieve ridge: */
      Ridge* ridge = ridges[iridge];
      Vertex* vertex = li->head;
      li = li->tail;

      /* On demand construction of ridge volume, distance and normal
       * (tope->dim - 1): */
      if (ridge->vdn == NULL) {
        ridge->vdn = allocator_alloc(alc, (2 + d) * sizeof(double));

        /* Matrix of vertex coordinates: */
        for (int ivertex = 0; ivertex < d; ++ivertex) {
          Vertex* vertex = ridge->vertices[ivertex];
          memcpy(&verts[ivertex * d], vertex->position, d * sizeof(double));
        }

        /* Analyse ridge simplex: */
        analysesimplex(d, d, verts, ridge->vdn + 0, centroid);

        /* Construct normal: */
        memcpy(ridge->vdn + 2, centroid, d * sizeof(double));
        vec_sub(d, ridge->vdn + 2, currentfacet->centroid);
        for (int i = 0; i < d - 1; ++i) {
          double fac = vec_dot(d, ridge->vdn + 2, &verts[i * d]);
          vec_adds(d, ridge->vdn + 2, &verts[i * d], -fac);
        }
        vec_normalize(d, ridge->vdn + 2);

        /* Ridge plane distance: */
        ridge->vdn[1] = vec_dot(d, ridge->vdn + 2, centroid);
      }

      /* Sign of normal: */
      double sign = 1.0;
      if (vec_dot(d, currentfacet->centroid, ridge->vdn + 2) < ridge->vdn[1]) {
        /* outward pointing normal */
        sign = -1.0;
      }

      /* Height of interpolation point above ridge: */
      double h = sign * (vec_dot(d, xiprime, ridge->vdn + 2) - ridge->vdn[1]);

      /* Weight for this ridge: */
      weights[iridge] = h * ridge->vdn[0];
      totalweight += weights[iridge];

      /* Index for this ridge: */
      indices[iridge] = vertex->index;

      /* Keep track of minimum height (if not boundary ridge): */
      if (h < hmin && ridge->facets[0] != NULL && ridge->facets[1] != NULL) {
        hmin = h;
        minridge = ridge;
      }
    }

    /* If hmin insignificant, stop: */
    if (hmin > -EPS) {
      break;
    }

    /* Else, goto next facet (cross ridge): */
    if (minridge->facets[0] == currentfacet) {
      currentfacet = minridge->facets[1];
    }
    else {
      assert(minridge->facets[1] == currentfacet);
      currentfacet = minridge->facets[0];
    }
  }

  if (totalweight > 0.0) {
    vec_scale(tope->dim, weights, 1.0 / totalweight);
  }
}



Facet* tope_firstfacet(Tope* tope)
{
  return tope->firstfacet;
}



Facet* tope_facet_nextfacet(Facet* facet)
{
  return facet->next;
}



void tope_facet_getnormal(Tope* tope, Facet* facet, double* normal)
{
  int d = tope->dim;
  memcpy(normal, facet->normal, d * sizeof(double));
  for (int i = 0; i < d; ++i) {
    normal[i] /= tope->scales[i];
  }
  vec_normalize(d, normal);
}



void tope_facet_getcentroid(Tope* tope, Facet* facet, double* centroid)
{
  int d = tope->dim;
  double* fcentroid = facet->centroid;
  memcpy(centroid, fcentroid, d * sizeof(double));
  for (int i = 0; i < d; ++i) {
    centroid[i] *= tope->scales[i];
  }
  vec_add(tope->dim, centroid, tope->shift);
}



double tope_facet_getoffset(Tope* tope, Facet* facet)
{
  double* centroid = alloca(tope->dim * sizeof(double));
  double* normal = alloca(tope->dim * sizeof(double));
  tope_facet_getcentroid(tope, facet, centroid);
  tope_facet_getnormal(tope, facet, normal);
  return vec_dot(tope->dim, centroid, normal);
}



double tope_facet_getvolume(Tope* tope, Facet* facet)
{
  int d = tope->dim;
  double* s = tope->scales;
  double* n = facet->normal;
  double volume = facet->volume;

  /* Transformed normal: */
  double scale = 0.0;
  for (int i = 0; i < d; ++i) {
    scale += n[i] * n[i] / (s[i] * s[i]);
  }

  /* Determinant of S: */
  double prod = 1.0;
  for (int i = 0; i < d; ++i) {
    prod *= s[i];
  }

  volume *= prod * sqrt(scale);

  return volume;
}



int tope_facet_getnumvertices(Tope* tope, Facet* facet)
{
  (void)facet;
  return tope->dim;
}



void tope_vertex_getposition(Tope* tope, Vertex* vertex, double* position)
{
  double* scales = tope->scales;
  double* shift = tope->shift;
  int i = tope->dim;
  i -= tope->isdelaunay ? 1 : 0;
  while (i--) {
    position[i] = vertex->position[i] * scales[i] + shift[i];
  }
}



u64 tope_bytes_used(Tope* tope)
{
  return tope->alc.used;
}