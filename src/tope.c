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

typedef tope_Vertex Vertex;
typedef tope_Facet Facet;



/* Create and initialize a vertex struct with a position and id: */
static u32 vertex_new(Tope* tope, double const* pos, int index)
{
  Allocator* alc = &tope->alc;

  /* Allocate vertex: */
  u32 hvertex = allocator_alloc(
      alc, sizeof(Vertex) - sizeof(double) + tope->dim * sizeof(double)
  );
  Vertex* vertex = allocator_mem(alc, hvertex);

  /* Prepend to list: */
  vertex->next = tope->firstvertex;
  vertex->prev = UINT32_MAX;
  if (vertex->next != UINT32_MAX) {
    Vertex* next = allocator_mem(alc, vertex->next);
    next->prev = hvertex;
  }
  tope->firstvertex = hvertex;
  ++tope->nverts;

  /* Fill in other variables: */
  vertex->tope = tope;
  vertex->index = index;
  vertex->nridges = 0;
  memcpy(vertex->position, pos, tope->dim * sizeof(double));

  return hvertex;
}



static void vertex_free(Tope* tope, u32 vertex)
{
  allocator_free(
      &tope->alc, vertex,
      sizeof(Vertex) + (tope->dim - 1) * sizeof(double)
  );
}



static void vertex_remove(Tope* tope, u32 hvertex)
{
  Vertex* vertex = allocator_mem(&tope->alc, hvertex);
  if (vertex->next != UINT32_MAX) {
    Vertex* next = allocator_mem(&tope->alc, vertex->next);
    next->prev = vertex->prev;
  }
  if (vertex->prev != UINT32_MAX) {
    Vertex* prev = allocator_mem(&tope->alc, vertex->prev);
    prev->next = vertex->next;
  }
  else {
    tope->firstvertex = vertex->next;
  }
  --tope->nverts;

  vertex_free(tope, hvertex);
}



/* Create and initialize a ridge struct with vertices: */
static u32 ridge_create(Tope* tope, u32* vertices)
{
  Allocator* alc = &tope->alc;

  /* Create ridge: */
  u32 hridge =
      allocator_alloc(alc, sizeof(Ridge) + (tope->dim - 2) * sizeof(u32));
  Ridge* ridge = allocator_mem(alc, hridge);
  ridge->next = UINT32_MAX;
  ridge->prev = UINT32_MAX;
  ridge->hvdn = UINT32_MAX;
  ridge->facets[0] = UINT32_MAX;
  ridge->facets[1] = UINT32_MAX;
  memcpy(ridge->vertices, vertices, (tope->dim - 1) * sizeof(u32));

  /* Associate vertices with ridge: */
  for (int i = 0; i < tope->dim - 1; ++i) {
    Vertex* vertex = allocator_mem(alc, vertices[i]);
    ++vertex->nridges;
  }

  /* Append ridge to tope list: */
  ridge->prev = tope->lastridge;
  ridge->next = UINT32_MAX;
  tope->lastridge = hridge;
  if (ridge->prev != UINT32_MAX) {
    Ridge* prev = allocator_mem(alc, ridge->prev);
    prev->next = hridge;
  }
  else {
    tope->firstridge = hridge;
  }
  ++tope->nridges;

  return hridge;
}



static void ridge_remove(Tope* tope, u32 hridge)
{
  Allocator* alc = &tope->alc;
  Ridge* ridge = allocator_mem(alc, hridge);
  if (ridge->prev != UINT32_MAX) {
    Ridge* prev = allocator_mem(alc, ridge->prev);
    prev->next = ridge->next;
  }
  else {
    assert(tope->firstridge == hridge);
    tope->firstridge = ridge->next;
  }
  if (ridge->next != UINT32_MAX) {
    Ridge* next = allocator_mem(alc, ridge->next);
    next->prev = ridge->prev;
  }
  else {
    assert(tope->lastridge == hridge);
    tope->lastridge = ridge->prev;
  }
  --tope->nridges;

  if (ridge->hvdn != UINT32_MAX) {
    allocator_free(alc, ridge->hvdn, (tope->dim + 1) * sizeof(double));
  }
  allocator_free(
      alc, hridge,
      sizeof(Ridge) - sizeof(u32) + (tope->dim - 1) * sizeof(u32)
  );
}



/* Create and initialize a facet struct: */
static u32 facet_new(Tope* tope)
{
  Allocator* alc = &tope->alc;
  int d = tope->dim;

  u32 hfacet = allocator_alloc(
      alc, sizeof(Facet) - sizeof(double) +
               d * (2 * sizeof(double) + 2 * sizeof(u32))
  );
  ++tope->nfacets;
  Facet* facet = allocator_mem(alc, hfacet);

  /* Append to list: */
  facet->prev = tope->lastfacet;
  if (facet->prev != UINT32_MAX) {
    Facet* prev = allocator_mem(alc, facet->prev);
    prev->next = hfacet;
  }
  else {
    tope->firstfacet = hfacet;
  }
  tope->lastfacet = hfacet;
  facet->next = UINT32_MAX;

  facet->tope = tope;
  facet->volume = 0.0;
  facet->dist = 0.0;
  facet->outsidehead = NULL;
  facet->outsidetail = NULL;
  facet->visible = false;
  vec_reset(d, facet->centroid);
  double* normal = facet->centroid + d;
  vec_reset(d, normal);
  u32* ridges = (u32*)(normal + d);
  memset(ridges, 0xff, d * sizeof(u32));
  u32* vertices = ridges + d;
  memset(vertices, 0xff, d * sizeof(u32));

  return hfacet;
}



static void facet_free(Tope* tope, u32 facet)
{
  allocator_free(
      &tope->alc, facet,
      sizeof(Facet) - sizeof(double) +
          tope->dim * (2 * sizeof(double) + 2 * sizeof(u32))
  );
}



static void facet_remove(Tope* tope, u32 hfacet)
{
  Allocator* alc = &tope->alc;
  Facet* facet = allocator_mem(alc, hfacet);

  if (facet->next != UINT32_MAX) {
    Facet* next = allocator_mem(alc, facet->next);
    next->prev = facet->prev;
  }
  else {
    assert(tope->lastfacet == hfacet);
    tope->lastfacet = facet->prev;
  }

  if (facet->prev != UINT32_MAX) {
    Facet* prev = allocator_mem(alc, facet->prev);
    prev->next = facet->next;
  }
  else {
    assert(tope->firstfacet == hfacet);
    tope->firstfacet = facet->next;
  }

  facet->next = UINT32_MAX;
  facet->prev = UINT32_MAX;
  --tope->nfacets;
}



static void facet_addoutside(Tope* tope, u32 hfacet, Point* point)
{
  Allocator* alc = &tope->alc;
  Facet* facet = allocator_mem(alc, hfacet);
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
  if (facet->prev != UINT32_MAX) {
    Facet* prev = allocator_mem(alc, facet->prev);
    /* Unlink from facet list: */
    if (facet->next != UINT32_MAX) {
      Facet* next = allocator_mem(alc, facet->next);
      next->prev = facet->prev;
    }
    else {
      assert(tope->lastfacet == hfacet);
      tope->lastfacet = facet->prev;
    }
    prev->next = facet->next;

    Facet* firstfacet = allocator_mem(alc, tope->firstfacet);
    if (firstfacet->outsidehead != NULL &&
        firstfacet->outsidehead->height > facet->outsidehead->height) {
      /* Move to second: */
      facet->next = firstfacet->next;
      if (facet->next != UINT32_MAX) {
        Facet* next = allocator_mem(alc, facet->next);
        next->prev = hfacet;
      }
      facet->prev = tope->firstfacet;
      firstfacet->next = hfacet;
    }
    else {
      /* Move to front: */
      facet->next = tope->firstfacet;
      firstfacet->prev = hfacet;
      tope->firstfacet = hfacet;
      facet->prev = UINT32_MAX;
    }
  }
}



/* Integrate facet into tope (compute volume, create ridges, etc): */
static void addfacet(Tope* tope, u32 hfacet, u32* hridges)
{
  int d = tope->dim;
  Allocator* alc = &tope->alc;
  Facet* facet = allocator_mem(alc, hfacet);

  /* Facet centroid: */
  vec_reset(d, facet->centroid);
  u32* vertices = (u32*)(facet->centroid + 2 * d) + d;
  for (int i = 0; i < d; ++i) {
    Vertex* vertex = allocator_mem(alc, vertices[i]);
    vec_add(d, facet->centroid, vertex->position);
  }
  vec_scale(d, facet->centroid, 1.0 / (double)d);

  /* Vertex span: */
  double* dirs = alloca((d - 1) * d * sizeof(double));
  Vertex* vertex = allocator_mem(alc, vertices[0]);
  double* x0 = vertex->position;
  for (int i = 0; i < d - 1; ++i) {
    /* Vertex position difference: */
    vertex = allocator_mem(alc, vertices[i + 1]);
    memcpy(dirs + i * d, vertex->position, d * sizeof(double));
    vec_sub(d, dirs + i * d, x0);
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
  double* normal = facet->centroid + d;
  memcpy(normal, facet->centroid, d * sizeof(double));
  vec_sub(d, normal, tope->center);
  for (int i = 0; i < d - 1; ++i) {
    double ip = vec_dot(d, normal, dirs + i * d);
    vec_adds(d, normal, dirs + i * d, -ip);
  }
  vec_normalize(d, normal);

  /* Plane distance: */
  facet->dist = vec_dot(d, normal, facet->centroid);

  /* Assign ridges: */
  u32* ridges = (u32*)(facet->centroid + 2 * d);
  for (int i = 0; i < d; ++i) {
    /* Get ridge: */
    Ridge* ridge = allocator_mem(alc, hridges[i]);

    /* Remember facet. */
    if (ridge->facets[0] == UINT32_MAX) {
      ridge->facets[0] = hfacet;
    }
    else {
      assert(ridge->facets[1] == UINT32_MAX);
      ridge->facets[1] = hfacet;
    }

    ridges[i] = hridges[i];
  }
}



/* Create initial simplex, add outside vertices, etc: */
static void initialsimplex(Tope* tope, int npoints, Point* points)
{
  int i;
  int j;
  int k;
  Allocator* alc = &tope->alc;

  /* For readability: */
  int d = tope->dim;

  /* Allocation: */
  double* span = malloc((npoints - 1) * d * sizeof(double));

  /* Point furthest from unit box center: */
  int maxindex = -1;
  double maxdist = 0.0;
  for (i = 0; i < npoints; ++i) {
    double dist = 0.0;
    for (j = 0; j < d; ++j) {
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
  for (i = 1; i < npoints; ++i) {
    memcpy(&span[(i - 1) * d], points[i].pos, d * sizeof(double));
    vec_sub(d, &span[(i - 1) * d], points[0].pos);
  }

  /* Permutation vector: */
  int* p = malloc(npoints * sizeof(int));
  for (i = 0; i < npoints; ++i) {
    p[i] = i;
  }

  /* Due to pivoting, the permutation vector will represent the point order
     which leads to the greatest volume. */
  for (i = 0; i < d; ++i) {
    /* Max norm: */
    double maxnormsq = 0.0;
    int pivot = -1;
    for (j = i; j < npoints - 1; ++j) {
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
    for (j = i + 1; j < npoints - 1; ++j) {
      double fac = vec_dot(d, span + i * d, span + j * d) / maxnormsq;
      vec_adds(d, span + j * d, span + i * d, -fac);
    }
  }

  /* Initialize tope center to initial simplex centroid: */
  vec_reset(tope->dim, tope->center);
  for (i = 0; i < tope->dim + 1; ++i) {
    vec_add(tope->dim, tope->center, points[p[i]].pos);
  }
  vec_scale(tope->dim, tope->center, 1.0 / (double)(tope->dim + 1));

  /* Create initial simplex vertices: */
  u32* hvertices = alloca((tope->dim + 1) * sizeof(u32));
  for (i = 0; i < tope->dim + 1; ++i) {
    hvertices[i] = vertex_new(tope, points[p[i]].pos, points[p[i]].index);
  }

  u32* hfacetridges = alloca(tope->dim * sizeof(u32));
  u32* hridgeverts = alloca((tope->dim - 1) * sizeof(u32));

  hashmap_clear(&tope->newridges);

  /* Create facets: */
  for (i = 0; i < tope->dim + 1; ++i) {
    /* Create facet: */
    u32 hfacet = facet_new(tope);
    Facet* facet = allocator_mem(alc, hfacet);
    u32* vertices = (u32*)(facet->centroid + 2 * d) + d;

    /* Add all vertices except for one: */
    for (j = 0; j < i; ++j) {
      vertices[j] = hvertices[j];
    }
    for (j = i + 1; j < tope->dim + 1; ++j) {
      vertices[j - 1] = hvertices[j];
    }

    /* Add ridges: */
    for (j = 0; j < d; ++j) {
      /* Ridge verts: */
      for (k = 0; k < j; ++k) {
        hridgeverts[k] = vertices[k];
      }
      for (k = j + 1; k < d; ++k) {
        hridgeverts[k - 1] = vertices[k];
      }

      u32* phnewridge = hashmap_get(&tope->newridges, d, hridgeverts, alc);
      if (*phnewridge == UINT32_MAX) {
        /* new ridge apparently */
        *phnewridge = ridge_create(tope, hridgeverts);
      }
      hfacetridges[j] = *phnewridge;
    }

    /* Vertex array filled, compute and integrate the facet: */
    addfacet(tope, hfacet, hfacetridges);
  }

  /* Create outside sets (assign each remaining vertex to a facet which it
   * 'sees'): */
  for (i = tope->dim + 1; i < npoints; ++i) {
    /* Determine maximum distance: */
    u32 hfacet = tope->firstfacet;
    Facet* facet = allocator_mem(alc, hfacet);
    while (hfacet != UINT32_MAX) {
      /* Distance to facet: */
      double* normal = facet->centroid + d;
      double h = vec_dot(tope->dim, normal, points[p[i]].pos) - facet->dist;

      /* If above facet, add it to outside set and move to next point. */
      if (h > EPS) {
        points[p[i]].height = h;
        facet_addoutside(tope, hfacet, points + p[i]);
        break;
      }

      hfacet = facet->next;
      facet = allocator_mem(alc, hfacet);
    }
  }

  /* Clean up: */
  free(p);
  free(span);
}



/* Add point to tope: */
static void addpoint(Tope* tope, u32 hfacet, Point* apex)
{
  int d = tope->dim;
  Allocator* alc = &tope->alc;
  int i;
  int ivertex;
  int iridge;
  double nh;

  facet_remove(tope, hfacet);
  Facet* facet = allocator_mem(alc, hfacet);
  facet->visible = true;

  /* Add facet to visible list: */
  u32 visiblelist = hfacet;

  /* Initialize outside points list: */
  Point* outsidepoints = NULL;

  while (visiblelist != UINT32_MAX) {
    /* Pop a value: */
    hfacet = visiblelist;
    facet = allocator_mem(alc, hfacet);
    visiblelist = facet->next;

    /* Loop over ridges to visit neighbours: */
    u32* ridges = (u32*)(facet->centroid + 2 * d);
    for (i = 0; i < d; ++i) {
      /* Ridge i: */
      u32 hridge = ridges[i];
      Ridge* ridge = allocator_mem(alc, hridge);

      /* Retrieve neighbour: */
      u32 hneighbour = UINT32_MAX;
      if (hfacet == ridge->facets[0]) {
        /* Current facet is 0, neighbour is 1: */
        ridge->facets[0] = UINT32_MAX; /* forget reference to this facet */
        hneighbour = ridge->facets[1];
      }
      else {
        assert(hfacet == ridge->facets[1]);
        /* Current facet is 1, neighbour is 0: */
        ridge->facets[1] = UINT32_MAX; /* forget reference to this facet */
        hneighbour = ridge->facets[0];
      }

      if (hneighbour == UINT32_MAX) {
        /* This neighbour was already removed, remove ridge. */

        /* Remove reference to this ridge from adjacent vertices: */
        for (ivertex = 0; ivertex < tope->dim - 1; ++ivertex) {
          u32 hvertex = ridge->vertices[ivertex];
          Vertex* vertex = allocator_mem(alc, hvertex);
          --vertex->nridges;
          if (vertex->nridges == 0) {
            /* Remove the vertex as well: */
            vertex_remove(tope, hvertex);
          }
        }

        /* Remove the ridge from the tope: */
        ridge_remove(tope, hridge);
      }
      else {
        Facet* neighbour = allocator_mem(alc, hneighbour);
        if (!neighbour->visible) {
          /* Neighbour is untested. */

          /* Height of apex above neighbour: */
          double* normal = neighbour->centroid + d;
          double h = vec_dot(d, apex->pos, normal) - neighbour->dist;

          /* Decide if neighbour is visible: */
          if (h > EPS) {
            /* Remove neighbour from tope list: */
            facet_remove(tope, hneighbour);

            /* Prepend neighbour to visible list: */
            neighbour->next = visiblelist;
            visiblelist = hneighbour;
            neighbour->visible = true;
          }
          else {
            /* Neighbour not visible. Ridge belongs to horizon: */
            if (tope->horizonridges_len == tope->horizonridges_cap) {
              tope->horizonridges_cap =
                  3 * tope->horizonridges_cap / 2 + 1;
              tope->horizonridges = realloc(
                  tope->horizonridges,
                  tope->horizonridges_cap * sizeof(Ridge*)
              );
            }
            tope->horizonridges[tope->horizonridges_len++] = hridge;
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
    facet_free(tope, hfacet);
  }

  /* Add apex vertex to tope: */
  u32 hvertex = vertex_new(tope, apex->pos, apex->index);
  u32* hridgeverts = alloca((tope->dim - 1) * sizeof(u32));
  u32* hfacetridges = alloca(tope->dim * sizeof(u32));

  hashmap_clear(&tope->newridges);

  /* Form new facets: */
  while (tope->horizonridges_len > 0) {
    /* New facet: */
    hfacet = facet_new(tope);
    facet = allocator_mem(alc, hfacet);

    /* Pop a horizon ridge: */
    u32 hhorizonridge = tope->horizonridges[--tope->horizonridges_len];
    Ridge* horizonridge = allocator_mem(alc, hhorizonridge);

    /* Facet vertices are apex + horizon vertices: */
    u32* vertices = (u32*)(facet->centroid + 2 * d) + d;
    vertices[0] = hvertex;
    for (ivertex = 1; ivertex < d; ++ivertex) {
      vertices[ivertex] = horizonridge->vertices[ivertex - 1];
    }

    hfacetridges[0] = hhorizonridge;
    for (iridge = 1; iridge < d; ++iridge) {
      /* Add all facet vertices except the one at index 'iridge': */
      for (ivertex = 0; ivertex < iridge; ++ivertex) {
        hridgeverts[ivertex] = vertices[ivertex];
      }
      for (ivertex = iridge + 1; ivertex < d; ++ivertex) {
        hridgeverts[ivertex - 1] = vertices[ivertex];
      }

      /* Get or create new ridge: */
      u32* phnewridge = hashmap_get(&tope->newridges, d, hridgeverts, alc);
      if (*phnewridge == UINT32_MAX) {
        /* Create new ridge: */
        *phnewridge = ridge_create(tope, hridgeverts);
        facet = allocator_mem(alc, hfacet);
        vertices = (u32*)(facet->centroid + 2 * d) + d;
        horizonridge = allocator_mem(alc, hhorizonridge);
      }
      hfacetridges[iridge] = *phnewridge;
    }

    addfacet(tope, hfacet, hfacetridges);

    /* Add new facet to newfacets array: */
    if (tope->newfacets_len == tope->newfacets_cap) {
      tope->newfacets_cap = 3 * tope->newfacets_cap / 2 + 1;
      tope->newfacets = realloc(
          tope->newfacets, tope->newfacets_cap * sizeof(Facet*)
      );
    }
    tope->newfacets[tope->newfacets_len++] = hfacet;
  }

  /* Assign outside verts: */
  hfacet = tope->newfacets[0];
  while (outsidepoints) {
    /* Remember next in list: */
    Point* next = outsidepoints->next;
    facet = allocator_mem(alc, hfacet);

    /* Height of outside point above facet: */
    double* normal = facet->centroid + d;
    double h = vec_dot(tope->dim, outsidepoints->pos, normal) - facet->dist;

    /* Visit neighbours, skipping horizon ridge: */
    u32* ridges = (u32*)(facet->centroid + 2 * d);
    u32 prev = UINT32_MAX;
    for (iridge = 1; iridge < d; ++iridge) {
      u32 hridge = ridges[iridge];
      Ridge* ridge = allocator_mem(alc, hridge);
      u32 hneighbour = UINT32_MAX;
      if (ridge->facets[0] == hfacet) {
        hneighbour = ridge->facets[1];
      }
      else {
        assert(ridge->facets[1] == hfacet);
        hneighbour = ridge->facets[0];
      }
      if (hneighbour == prev) {
        continue;
      }

      /* Neighbour height: */
      Facet* neighbour = allocator_mem(alc, hneighbour);
      double* normal = neighbour->centroid + d;
      nh = vec_dot(tope->dim, outsidepoints->pos, normal) - neighbour->dist;
      if (nh > h) {
        prev = hfacet;
        hfacet = hneighbour;
        facet = allocator_mem(alc, hfacet);
        ridges = (u32*)(facet->centroid + 2 * d);
        h = nh;
        iridge = 0; /* reset iteration (new base facet) */
      }
    }

    if (h > EPS) {
      /* Outside. */
      outsidepoints->height = h;
      facet_addoutside(tope, hfacet, outsidepoints);
    }

    /* Next in list: */
    outsidepoints = next;
  }

  assert(tope->horizonridges_len == 0);
  tope->newfacets_len = 0;
}



static void build(Tope* tope, int npoints, Point* points)
{
  Allocator* alc = &tope->alc;
  assert(npoints >= tope->dim + 1);

  initialsimplex(tope, npoints, points);

  /* Now, all vertices are either on the initial simplex OR
     added to the outside set of one of its facets OR
     inside the initial simplex. */
  Facet* firstfacet = allocator_mem(alc, tope->firstfacet);
  while (firstfacet->outsidehead) {
    /* Pop head: */
    Point* point = firstfacet->outsidehead;
    firstfacet->outsidehead = point->next;

    /* Add vertex to tope: */
    addpoint(tope, tope->firstfacet, point);
    firstfacet = allocator_mem(alc, tope->firstfacet);
  }
}



static Tope* tope_new()
{
  Tope* tope = calloc(1, sizeof(Tope));

  allocator_init(&tope->alc);
  hashmap_init(&tope->newridges);
  tope->firstfacet = UINT32_MAX;
  tope->lastfacet = UINT32_MAX;
  tope->firstridge = UINT32_MAX;
  tope->lastridge = UINT32_MAX;
  tope->firstvertex = UINT32_MAX;

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

  for (u32 hfacet = rectope->firstfacet, ifacet = 0; hfacet != UINT32_MAX;
       hfacet = ((Facet*)allocator_mem(&rectope->alc, hfacet))->next,
           ++ifacet) {
    /* Unscaled normal and distance: */
    Facet* facet = allocator_mem(&rectope->alc, hfacet);
    double* normal = facet->centroid + d;
    double dist = facet->dist;
    for (int i = 0; i < d; ++i) {
      points[ifacet * d + i] = normal[i] / rectope->scales[i];
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
  tope->shift = malloc(tope->dim * sizeof(double));
  tope->scales = malloc(tope->dim * sizeof(double));
  tope->center = malloc(tope->dim * sizeof(double));

  /* Bounding box: */
  minindices = alloca(tope->dim * sizeof(int));
  maxindices = alloca(tope->dim * sizeof(int));
  minima = alloca(tope->dim * sizeof(double));
  maxima = alloca(tope->dim * sizeof(double));
  boundingbox(
      npoints, tope->dim, orgpoints, minindices, maxindices, minima, maxima
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
  int ipoint;
  int i;
  int j;
  int k;
  int* minindices;
  int* maxindices;
  double* positions;
  Point* points;

  Tope* tope = tope_new();
  Allocator* alc = &tope->alc;

  /* Reciprocal of dimension: */
  double recd = 1.0 / (double)dim;

  /* Initialize member variables: */
  tope->dim = dim + 1;
  tope->isdelaunay = true;
  tope->shift = malloc(tope->dim * sizeof(double));
  tope->scales = malloc(tope->dim * sizeof(double));
  tope->center = malloc(tope->dim * sizeof(double));

  /* Bounding box: */
  minindices = alloca(dim * sizeof(int));
  maxindices = alloca(dim * sizeof(int));
  boundingbox(
      npoints, dim, orgpoints, minindices, maxindices, tope->shift,
      tope->scales
  );
  tope->shift[dim] = 0.0;

  /* Scales: */
  vec_sub(dim, tope->scales, tope->shift);
  tope->scales[dim] = 1.0;

  /* Points array: */
  points = malloc((npoints + 1) * sizeof(Point));
  positions = malloc((npoints + 1) * (dim + 1) * sizeof(double));
  for (ipoint = 0; ipoint < npoints + 1; ++ipoint) {
    points[ipoint].next = NULL;
    points[ipoint].index = ipoint;
    points[ipoint].height = 0.0;
    points[ipoint].pos = positions + ipoint * (dim + 1);

    if (ipoint < npoints) {
      /* Transformed point: */
      for (i = 0; i < dim; ++i) {
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
  for (ipoint = 0; ipoint < npoints; ++ipoint) {
    vec_add(dim, points[npoints].pos, points[ipoint].pos);
  }
  vec_scale(dim, points[npoints].pos, 1.0 / (double)npoints);
  points[npoints].pos[dim] = 2.0;

  /* Construct convex tope of paraboloid: */
  build(tope, npoints + 1, points);

  /* Clean up: */
  free(positions);
  free(points);

  /* Remove upper delaunay surfaces: */
  u32* hridges = malloc(tope->nridges * sizeof(u32));
  u32 hridge = tope->firstridge;
  for (i = 0; i < tope->nridges;
       ++i, hridge = ((Ridge*)allocator_mem(alc, hridge))->next) {
    hridges[i] = hridge;
  }

  for (i = 0; i < tope->nridges; ++i) {
    hridge = hridges[i];
    Ridge* ridge = allocator_mem(alc, hridge);
    for (j = 0; j < 2; ++j) {
      u32 hfacet = ridge->facets[j];
      if (hfacet == UINT32_MAX) {
        continue;
      }

      Facet* facet = allocator_mem(alc, hfacet);
      double* normal = facet->centroid + tope->dim;
      if (normal[tope->dim - 1] < 0.0) {
        continue;
      }

      /* Upper delaunay. */

      /* Disassociate facet from ridges: */
      ridge->facets[j] = UINT32_MAX;

      u32* ridges = (u32*)(facet->centroid + 2 * tope->dim);
      for (k = 0; k < tope->dim; ++k) {
        u32 hridgek = ridges[k];
        Ridge* ridgek = allocator_mem(alc, hridgek);
        if (ridgek->facets[0] == hfacet) {
          ridgek->facets[0] = UINT32_MAX;
        }
        else if (ridgek->facets[1] == hfacet) {
          ridgek->facets[1] = UINT32_MAX;
        }
      }

      /* Remove facet: */
      facet_remove(tope, hfacet);
      facet_free(tope, hfacet);
    }

    if (ridge->facets[0] == UINT32_MAX && ridge->facets[1] == UINT32_MAX) {
      /* Both adjacent facets are gone. */

      /* Disassociate from adjacent vertices: */
      for (j = 0; j < tope->dim - 1; ++j) {
        Vertex* vertex = allocator_mem(alc, ridge->vertices[j]);
        --vertex->nridges;
        if (vertex->nridges == 0) {
          /* Remove vertex as well: */
          vertex_remove(tope, ridge->vertices[j]);
        }
      }

      /* Remove ridge: */
      ridge_remove(tope, hridge);
      hridges[i] = UINT32_MAX; /* forget reference */
    }
  }
  free(hridges);

  return tope;
}



void tope_addvertex(Tope* tope, double const* point)
{
  int idim;
  Allocator* alc = &tope->alc;

  /* Allocations: */
  double* pos = alloca(tope->dim * sizeof(double));

  /* Initialize apex: */
  Point apex;
  apex.next = NULL;
  apex.index = tope->nverts;
  apex.height = 0.0;
  apex.pos = pos;

  /* Transformed position: */
  memcpy(pos, point, tope->dim * sizeof(double));
  vec_sub(tope->dim, pos, tope->shift);
  for (idim = 0; idim < tope->dim; ++idim) {
    if (tope->scales[idim] > 0.0) {
      pos[idim] /= tope->scales[idim];
    }
  }

  /* Find facet that can see the point: */
  u32 hmaxfacet = UINT32_MAX;
  u32 hfacet = tope->firstfacet;
  while (hfacet != UINT32_MAX) {
    Facet* facet = allocator_mem(alc, hfacet);
    double* normal = facet->centroid + tope->dim;
    double h = vec_dot(tope->dim, pos, normal) - facet->dist;
    if (h > EPS) {
      hmaxfacet = hfacet;
      break;
    }
    hfacet = facet->next;
  }

  if (hmaxfacet != UINT32_MAX) {
    /* Add vertex to tope: */
    addpoint(tope, hmaxfacet, &apex);
  }
}



void tope_print(Tope* tope)
{
  Allocator* alc = &tope->alc;

  int d = tope->dim;
  double* sv = alloca(d * sizeof(double));
  double* normal = alloca(d * sizeof(double));
  double* centroid = alloca(d * sizeof(double));
  double* position = alloca(d * sizeof(double));

  printf("%d facets\n", tope->nfacets);
  printf("%d ridges\n", tope->nridges);
  printf("%d vertices\n", tope->nverts);

  vec_reset(d, sv);

  int i = 0;
  u32 hfacet = tope->firstfacet;
  while (hfacet != UINT32_MAX) {
    Facet* facet = allocator_mem(alc, hfacet);
    tope_facet_getnormal(facet, normal);
    double volume = tope_facet_getvolume(facet);
    vec_adds(d, sv, normal, volume);
    tope_facet_getcentroid(facet, centroid);

    printf("facet %d\n", i + 1);
    printf("  volume: %g\n", volume);
    printf("  centroid: ");
    vec_print(d, centroid);
    printf("\n");
    printf("  normal: ");
    vec_print(d, normal);
    printf("\n");

    printf("  vertices:\n");
    u32* vertices = (u32*)(facet->centroid + 2 * d) + d;
    for (int j = 0; j < d; ++j) {
      u32 hvertex = vertices[j];
      Vertex* vertex = allocator_mem(alc, hvertex);
      tope_vertex_getposition(vertex, position);
      vec_print(d, position);
    }
    printf("\n");

    ++i;
    hfacet = facet->next;
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
  int i;
  int ivertex;
  int iridge;
  double* xiprime;
  double* verts;
  double hmin;
  double h;
  Ridge* minridge;
  double totalweight;
  double sign;
  Ridge* ridge;
  Allocator* alc = &tope->alc;

  assert(tope->isdelaunay);

  /* Some space for various tasks: */
  int d = tope->dim - 1;
  xiprime = alloca(d * sizeof(double));
  double* centroid = alloca(d * sizeof(double));
  verts = alloca(d * d * sizeof(double));

  /* Transform xi to local coordinates: */
  for (i = 0; i < d; ++i) {
    xiprime[i] = (xi[i] - tope->shift[i]) / tope->scales[i];
  }

  /* Start with first facet: */
  assert(tope->firstfacet != UINT32_MAX);
  u32 hcurrentfacet = tope->firstfacet;

  while (1) {
    Facet* currentfacet = allocator_mem(alc, hcurrentfacet);
    u32* ridges = (u32*)(currentfacet->centroid + 2 * tope->dim);
    u32* vertices = ridges + tope->dim;

    /* In the current facet, find the ridge where xi is highest above: */
    totalweight = 0.0;
    hmin = HUGE_VAL;
    minridge = NULL;
    for (iridge = 0; iridge < tope->dim; ++iridge) {
      /* Retrieve ridge: */
      u32 hridge = ridges[iridge];
      Ridge* ridge = allocator_mem(alc, hridge);
      u32 hvertex = vertices[iridge];
      Vertex* vertex = allocator_mem(alc, hvertex);

      /* On demand construction of ridge volume, distance and normal
       * (tope->dim - 1): */
      double* vdn = NULL;
      if (ridge->hvdn != UINT32_MAX) {
        vdn = allocator_mem(alc, ridge->hvdn);
      }
      else {
        ridge->hvdn = allocator_alloc(alc, (2 + d) * sizeof(double));
        vdn = allocator_mem(alc, ridge->hvdn);

        /* Matrix of vertex coordinates: */
        for (ivertex = 0; ivertex < d; ++ivertex) {
          Vertex* vertex = allocator_mem(alc, ridge->vertices[ivertex]);
          memcpy(&verts[ivertex * d], vertex->position, d * sizeof(double));
        }

        /* Analyse ridge simplex: */
        analysesimplex(d, d, verts, vdn + 0, centroid);

        /* Construct normal: */
        memcpy(vdn + 2, centroid, d * sizeof(double));
        vec_sub(d, vdn + 2, currentfacet->centroid);
        for (i = 0; i < d - 1; ++i) {
          double fac = vec_dot(d, vdn + 2, &verts[i * d]);
          vec_adds(d, vdn + 2, &verts[i * d], -fac);
        }
        vec_normalize(d, vdn + 2);

        /* Ridge plane distance: */
        vdn[1] = vec_dot(d, vdn + 2, centroid);
      }

      /* Sign of normal: */
      sign = 1.0;
      if (vec_dot(d, currentfacet->centroid, vdn + 2) < vdn[1]) {
        /* outward pointing normal */
        sign = -1.0;
      }

      /* Height of interpolation point above ridge: */
      h = sign * (vec_dot(d, xiprime, vdn + 2) - vdn[1]);

      /* Weight for this ridge: */
      weights[iridge] = h * vdn[0];
      totalweight += weights[iridge];

      /* Index for this ridge: */
      indices[iridge] = vertex->index;

      /* Keep track of minimum height (if not boundary ridge): */
      if (h < hmin && ridge->facets[0] != UINT32_MAX &&
          ridge->facets[1] != UINT32_MAX) {
        hmin = h;
        minridge = ridge;
      }
    }

    /* If hmin insignificant, stop: */
    if (hmin > -EPS) {
      break;
    }

    /* Else, goto next facet (cross ridge): */
    ridge = minridge;
    if (ridge->facets[0] == hcurrentfacet) {
      hcurrentfacet = ridge->facets[1];
    }
    else {
      assert(ridge->facets[1] == hcurrentfacet);
      hcurrentfacet = ridge->facets[0];
    }
  }

  if (totalweight > 0.0) {
    vec_scale(tope->dim, weights, 1.0 / totalweight);
  }
}



Facet* tope_firstfacet(Tope* tope)
{
  return allocator_mem(&tope->alc, tope->firstfacet);
}



Vertex* tope_firstvertex(Tope* tope)
{
  return allocator_mem(&tope->alc, tope->firstvertex);
}



Facet* tope_facet_nextfacet(Facet* facet)
{
  return allocator_mem(&facet->tope->alc, facet->next);
}



void tope_facet_getnormal(Facet* facet, double* normal)
{
  int d = facet->tope->dim;
  memcpy(normal, facet->centroid + d, d * sizeof(double));
  for (int i = 0; i < d; ++i) {
    normal[i] /= facet->tope->scales[i];
  }
  vec_normalize(d, normal);
}



void tope_facet_getcentroid(Facet* facet, double* centroid)
{
  int i;

  memcpy(centroid, facet->centroid, facet->tope->dim * sizeof(double));
  for (i = 0; i < facet->tope->dim; ++i) {
    centroid[i] *= facet->tope->scales[i];
  }
  vec_add(facet->tope->dim, centroid, facet->tope->shift);
}



double tope_facet_getoffset(Facet* facet)
{
  double* centroid = alloca(facet->tope->dim * sizeof(double));
  double* normal = alloca(facet->tope->dim * sizeof(double));
  tope_facet_getcentroid(facet, centroid);
  tope_facet_getnormal(facet, normal);
  return vec_dot(facet->tope->dim, centroid, normal);
}



double tope_facet_getvolume(Facet* facet)
{
  int d = facet->tope->dim;
  double* s = facet->tope->scales;
  double* n = facet->centroid + d;
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



int tope_facet_getnumvertices(Facet* facet)
{
  return facet->tope->dim;
}



Vertex* tope_vertex_nextvertex(Vertex* vertex)
{
  return allocator_mem(&vertex->tope->alc, vertex->next);
}



void tope_vertex_getposition(Vertex* vertex, double* position)
{
  double* scales = vertex->tope->scales;
  double* shift = vertex->tope->shift;
  int i = vertex->tope->dim;
  i -= vertex->tope->isdelaunay ? 1 : 0;
  while (i--) {
    position[i] = vertex->position[i] * scales[i] + shift[i];
  }
}
