#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "hashmap.h"
#include "math.h"
#include "util.h"

#include <polytoop.h>

#define EPS 1.0e-9

typedef polytoop_Vertex Vertex;
typedef polytoop_Facet Facet;



/* Create and initialize a vertex struct with a position and id: */
static u32 vertex_new(Polytoop* polytoop, double const* pos, int index)
{
  Allocator* alc = &polytoop->alc;

  /* Allocate vertex: */
  u32 hvertex = allocator_alloc(
      alc, sizeof(Vertex) - sizeof(double) + polytoop->dim * sizeof(double)
  );
  Vertex* vertex = allocator_mem(alc, hvertex);

  /* Prepend to list: */
  vertex->next = polytoop->firstvertex;
  vertex->prev = UINT32_MAX;
  if (vertex->next != UINT32_MAX) {
    Vertex* next = allocator_mem(alc, vertex->next);
    next->prev = hvertex;
  }
  polytoop->firstvertex = hvertex;
  ++polytoop->nverts;

  /* Fill in other variables: */
  vertex->polytoop = polytoop;
  vertex->index = index;
  vertex->nridges = 0;
  memcpy(vertex->position, pos, polytoop->dim * sizeof(double));

  return hvertex;
}



static void vertex_free(Polytoop* polytoop, u32 vertex)
{
  allocator_free(
      &polytoop->alc, vertex,
      sizeof(Vertex) + (polytoop->dim - 1) * sizeof(double)
  );
}



static void vertex_remove(Polytoop* polytoop, u32 hvertex)
{
  Vertex* vertex = allocator_mem(&polytoop->alc, hvertex);
  if (vertex->next != UINT32_MAX) {
    Vertex* next = allocator_mem(&polytoop->alc, vertex->next);
    next->prev = vertex->prev;
  }
  if (vertex->prev != UINT32_MAX) {
    Vertex* prev = allocator_mem(&polytoop->alc, vertex->prev);
    prev->next = vertex->next;
  }
  else {
    polytoop->firstvertex = vertex->next;
  }
  --polytoop->nverts;

  vertex_free(polytoop, hvertex);
}



/* Create and initialize a ridge struct with vertices: */
static u32 ridge_create(Polytoop* polytoop, u32* vertices)
{
  Allocator* alc = &polytoop->alc;

  /* Create ridge: */
  u32 hridge =
      allocator_alloc(alc, sizeof(Ridge) + (polytoop->dim - 2) * sizeof(u32));
  Ridge* ridge = allocator_mem(alc, hridge);
  ridge->next = UINT32_MAX;
  ridge->prev = UINT32_MAX;
  ridge->hvdn = UINT32_MAX;
  ridge->facets[0] = UINT32_MAX;
  ridge->facets[1] = UINT32_MAX;
  memcpy(ridge->vertices, vertices, (polytoop->dim - 1) * sizeof(u32));

  /* Associate vertices with ridge: */
  for (int i = 0; i < polytoop->dim - 1; ++i) {
    Vertex* vertex = allocator_mem(alc, vertices[i]);
    ++vertex->nridges;
  }

  /* Append ridge to polytoop list: */
  ridge->prev = polytoop->lastridge;
  ridge->next = UINT32_MAX;
  polytoop->lastridge = hridge;
  if (ridge->prev != UINT32_MAX) {
    Ridge* prev = allocator_mem(alc, ridge->prev);
    prev->next = hridge;
  }
  else {
    polytoop->firstridge = hridge;
  }
  ++polytoop->nridges;

  return hridge;
}



static void ridge_remove(Polytoop* polytoop, u32 hridge)
{
  Allocator* alc = &polytoop->alc;
  Ridge* ridge = allocator_mem(alc, hridge);
  if (ridge->prev != UINT32_MAX) {
    Ridge* prev = allocator_mem(alc, ridge->prev);
    prev->next = ridge->next;
  }
  else {
    assert(polytoop->firstridge == hridge);
    polytoop->firstridge = ridge->next;
  }
  if (ridge->next != UINT32_MAX) {
    Ridge* next = allocator_mem(alc, ridge->next);
    next->prev = ridge->prev;
  }
  else {
    assert(polytoop->lastridge == hridge);
    polytoop->lastridge = ridge->prev;
  }
  --polytoop->nridges;

  if (ridge->hvdn != UINT32_MAX) {
    allocator_free(alc, ridge->hvdn, (polytoop->dim + 1) * sizeof(double));
  }
  allocator_free(
      alc, hridge,
      sizeof(Ridge) - sizeof(u32) + (polytoop->dim - 1) * sizeof(u32)
  );
}



/* Create and initialize a facet struct: */
static u32 facet_new(Polytoop* polytoop)
{
  Allocator* alc = &polytoop->alc;
  int d = polytoop->dim;

  u32 hfacet = allocator_alloc(
      alc, sizeof(Facet) - sizeof(double) +
               d * (2 * sizeof(double) + 2 * sizeof(u32))
  );
  ++polytoop->nfacets;
  Facet* facet = allocator_mem(alc, hfacet);

  /* Append to list: */
  facet->prev = polytoop->lastfacet;
  if (facet->prev != UINT32_MAX) {
    Facet* prev = allocator_mem(alc, facet->prev);
    prev->next = hfacet;
  }
  else {
    polytoop->firstfacet = hfacet;
  }
  polytoop->lastfacet = hfacet;
  facet->next = UINT32_MAX;

  facet->polytoop = polytoop;
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



static void facet_free(Polytoop* polytoop, u32 facet)
{
  allocator_free(
      &polytoop->alc, facet,
      sizeof(Facet) - sizeof(double) +
          polytoop->dim * (2 * sizeof(double) + 2 * sizeof(u32))
  );
}



static void facet_remove(Polytoop* polytoop, u32 hfacet)
{
  Allocator* alc = &polytoop->alc;
  Facet* facet = allocator_mem(alc, hfacet);

  if (facet->next != UINT32_MAX) {
    Facet* next = allocator_mem(alc, facet->next);
    next->prev = facet->prev;
  }
  else {
    assert(polytoop->lastfacet == hfacet);
    polytoop->lastfacet = facet->prev;
  }

  if (facet->prev != UINT32_MAX) {
    Facet* prev = allocator_mem(alc, facet->prev);
    prev->next = facet->next;
  }
  else {
    assert(polytoop->firstfacet == hfacet);
    polytoop->firstfacet = facet->next;
  }

  facet->next = UINT32_MAX;
  facet->prev = UINT32_MAX;
  --polytoop->nfacets;
}



static void facet_addoutside(Polytoop* polytoop, u32 hfacet, Point* point)
{
  Allocator* alc = &polytoop->alc;
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
      assert(polytoop->lastfacet == hfacet);
      polytoop->lastfacet = facet->prev;
    }
    prev->next = facet->next;

    Facet* firstfacet = allocator_mem(alc, polytoop->firstfacet);
    if (firstfacet->outsidehead != NULL &&
        firstfacet->outsidehead->height > facet->outsidehead->height) {
      /* Move to second: */
      facet->next = firstfacet->next;
      if (facet->next != UINT32_MAX) {
        Facet* next = allocator_mem(alc, facet->next);
        next->prev = hfacet;
      }
      facet->prev = polytoop->firstfacet;
      firstfacet->next = hfacet;
    }
    else {
      /* Move to front: */
      facet->next = polytoop->firstfacet;
      firstfacet->prev = hfacet;
      polytoop->firstfacet = hfacet;
      facet->prev = UINT32_MAX;
    }
  }
}



/* Integrate facet into polytoop (compute volume, create ridges, etc): */
static void addfacet(Polytoop* polytoop, u32 hfacet, u32* hridges)
{
  int d = polytoop->dim;
  Allocator* alc = &polytoop->alc;
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
  vec_sub(d, normal, polytoop->center);
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
static void initialsimplex(Polytoop* polytoop, int npoints, Point* points)
{
  int i;
  int j;
  int k;
  Allocator* alc = &polytoop->alc;

  /* For readability: */
  int d = polytoop->dim;

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
  int* p = alloca(npoints * sizeof(int));
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

  /* Initialize polytoop center to initial simplex centroid: */
  vec_reset(polytoop->dim, polytoop->center);
  for (i = 0; i < polytoop->dim + 1; ++i) {
    vec_add(polytoop->dim, polytoop->center, points[p[i]].pos);
  }
  vec_scale(polytoop->dim, polytoop->center, 1.0 / (double)(polytoop->dim + 1));

  /* Create initial simplex vertices: */
  u32* hvertices = alloca((polytoop->dim + 1) * sizeof(u32));
  for (i = 0; i < polytoop->dim + 1; ++i) {
    hvertices[i] = vertex_new(polytoop, points[p[i]].pos, points[p[i]].index);
  }

  u32* hfacetridges = alloca(polytoop->dim * sizeof(u32));
  u32* hridgeverts = alloca((polytoop->dim - 1) * sizeof(u32));

  hashmap_clear(&polytoop->newridges);

  /* Create facets: */
  for (i = 0; i < polytoop->dim + 1; ++i) {
    /* Create facet: */
    u32 hfacet = facet_new(polytoop);
    Facet* facet = allocator_mem(alc, hfacet);
    u32* vertices = (u32*)(facet->centroid + 2 * d) + d;

    /* Add all vertices except for one: */
    for (j = 0; j < i; ++j) {
      vertices[j] = hvertices[j];
    }
    for (j = i + 1; j < polytoop->dim + 1; ++j) {
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

      u32* phnewridge = hashmap_get(&polytoop->newridges, d, hridgeverts, alc);
      if (*phnewridge == UINT32_MAX) {
        /* new ridge apparently */
        *phnewridge = ridge_create(polytoop, hridgeverts);
      }
      hfacetridges[j] = *phnewridge;
    }

    /* Vertex array filled, compute and integrate the facet: */
    addfacet(polytoop, hfacet, hfacetridges);
  }

  /* Create outside sets (assign each remaining vertex to a facet which it
   * 'sees'): */
  for (i = polytoop->dim + 1; i < npoints; ++i) {
    /* Determine maximum distance: */
    u32 hfacet = polytoop->firstfacet;
    Facet* facet = allocator_mem(alc, hfacet);
    while (hfacet != UINT32_MAX) {
      /* Distance to facet: */
      double* normal = facet->centroid + d;
      double h = vec_dot(polytoop->dim, normal, points[p[i]].pos) - facet->dist;

      /* If above facet, add it to outside set and move to next point. */
      if (h > EPS) {
        points[p[i]].height = h;
        facet_addoutside(polytoop, hfacet, points + p[i]);
        break;
      }

      hfacet = facet->next;
      facet = allocator_mem(alc, hfacet);
    }
  }

  /* Clean up: */
  free(span);
}



/* Add point to polytoop: */
static void addpoint(Polytoop* polytoop, u32 hfacet, Point* apex)
{
  int d = polytoop->dim;
  Allocator* alc = &polytoop->alc;
  int i;
  int ivertex;
  int iridge;
  double nh;

  facet_remove(polytoop, hfacet);
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
        for (ivertex = 0; ivertex < polytoop->dim - 1; ++ivertex) {
          u32 hvertex = ridge->vertices[ivertex];
          Vertex* vertex = allocator_mem(alc, hvertex);
          --vertex->nridges;
          if (vertex->nridges == 0) {
            /* Remove the vertex as well: */
            vertex_remove(polytoop, hvertex);
          }
        }

        /* Remove the ridge from the polytoop: */
        ridge_remove(polytoop, hridge);
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
            /* Remove neighbour from polytoop list: */
            facet_remove(polytoop, hneighbour);

            /* Prepend neighbour to visible list: */
            neighbour->next = visiblelist;
            visiblelist = hneighbour;
            neighbour->visible = true;
          }
          else {
            /* Neighbour not visible. Ridge belongs to horizon: */
            if (polytoop->horizonridges_len == polytoop->horizonridges_cap) {
              polytoop->horizonridges_cap =
                  3 * polytoop->horizonridges_cap / 2 + 1;
              polytoop->horizonridges = realloc(
                  polytoop->horizonridges,
                  polytoop->horizonridges_cap * sizeof(Ridge*)
              );
            }
            polytoop->horizonridges[polytoop->horizonridges_len++] = hridge;
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
    facet_free(polytoop, hfacet);
  }

  /* Add apex vertex to polytoop: */
  u32 hvertex = vertex_new(polytoop, apex->pos, apex->index);
  u32* hridgeverts = alloca((polytoop->dim - 1) * sizeof(u32));
  u32* hfacetridges = alloca(polytoop->dim * sizeof(u32));

  hashmap_clear(&polytoop->newridges);

  /* Form new facets: */
  while (polytoop->horizonridges_len > 0) {
    /* New facet: */
    hfacet = facet_new(polytoop);
    facet = allocator_mem(alc, hfacet);

    /* Pop a horizon ridge: */
    u32 hhorizonridge = polytoop->horizonridges[--polytoop->horizonridges_len];
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
      u32* phnewridge = hashmap_get(&polytoop->newridges, d, hridgeverts, alc);
      if (*phnewridge == UINT32_MAX) {
        /* Create new ridge: */
        *phnewridge = ridge_create(polytoop, hridgeverts);
        facet = allocator_mem(alc, hfacet);
        vertices = (u32*)(facet->centroid + 2 * d) + d;
        horizonridge = allocator_mem(alc, hhorizonridge);
      }
      hfacetridges[iridge] = *phnewridge;
    }

    addfacet(polytoop, hfacet, hfacetridges);

    /* Add new facet to newfacets array: */
    if (polytoop->newfacets_len == polytoop->newfacets_cap) {
      polytoop->newfacets_cap = 3 * polytoop->newfacets_cap / 2 + 1;
      polytoop->newfacets = realloc(
          polytoop->newfacets, polytoop->newfacets_cap * sizeof(Facet*)
      );
    }
    polytoop->newfacets[polytoop->newfacets_len++] = hfacet;
  }

  /* Assign outside verts: */
  hfacet = polytoop->newfacets[0];
  while (outsidepoints) {
    /* Remember next in list: */
    Point* next = outsidepoints->next;
    facet = allocator_mem(alc, hfacet);

    /* Height of outside point above facet: */
    double* normal = facet->centroid + d;
    double h = vec_dot(polytoop->dim, outsidepoints->pos, normal) - facet->dist;

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
      nh = vec_dot(polytoop->dim, outsidepoints->pos, normal) - neighbour->dist;
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
      facet_addoutside(polytoop, hfacet, outsidepoints);
    }

    /* Next in list: */
    outsidepoints = next;
  }

  assert(polytoop->horizonridges_len == 0);
  polytoop->newfacets_len = 0;
}



static void build(Polytoop* polytoop, int npoints, Point* points)
{
  Allocator* alc = &polytoop->alc;
  assert(npoints >= polytoop->dim + 1);

  initialsimplex(polytoop, npoints, points);

  /* Now, all vertices are either on the initial simplex OR
     added to the outside set of one of its facets OR
     inside the initial simplex. */
  Facet* firstfacet = allocator_mem(alc, polytoop->firstfacet);
  while (firstfacet->outsidehead) {
    /* Pop head: */
    Point* point = firstfacet->outsidehead;
    firstfacet->outsidehead = point->next;

    /* Add vertex to polytoop: */
    addpoint(polytoop, polytoop->firstfacet, point);
    firstfacet = allocator_mem(alc, polytoop->firstfacet);
  }
}



static Polytoop* polytoop_new()
{
  Polytoop* polytoop = calloc(1, sizeof(Polytoop));

  allocator_init(&polytoop->alc);
  hashmap_init(&polytoop->newridges);
  polytoop->firstfacet = UINT32_MAX;
  polytoop->lastfacet = UINT32_MAX;
  polytoop->firstridge = UINT32_MAX;
  polytoop->lastridge = UINT32_MAX;
  polytoop->firstvertex = UINT32_MAX;

  return polytoop;
}



void polytoop_delete(Polytoop* polytoop)
{
  free(polytoop->newfacets);
  free(polytoop->horizonridges);
  hashmap_destroy(&polytoop->newridges);
  allocator_destroy(&polytoop->alc);
  free(polytoop->center);
  free(polytoop->scales);
  free(polytoop->shift);
  free(polytoop);
}



Polytoop* polytoop_fromplanes(
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

  /* Construct reciprocal polytoop: */
  Polytoop* recpolytoop = polytoop_frompoints(n, d, points);
  free(points);

  /* Convert reciprocal polytoop to polytoop: */
  int nfacets = polytoop_getnumfacets(recpolytoop);
  points = malloc(nfacets * d * sizeof(double));

  for (u32 hfacet = recpolytoop->firstfacet, ifacet = 0; hfacet != UINT32_MAX;
       hfacet = ((Facet*)allocator_mem(&recpolytoop->alc, hfacet))->next,
           ++ifacet) {
    /* Unscaled normal and distance: */
    Facet* facet = allocator_mem(&recpolytoop->alc, hfacet);
    double* normal = facet->centroid + d;
    double dist = facet->dist;
    for (int i = 0; i < d; ++i) {
      points[ifacet * d + i] = normal[i] / recpolytoop->scales[i];
      dist += points[ifacet * d + i] * recpolytoop->shift[i];
    }

    /* Straight space point: */
    vec_scale(d, points + ifacet * d, 1.0 / dist);
    vec_add(d, points + ifacet * d, xc);
  }

  Polytoop* polytoop = polytoop_frompoints(nfacets, d, points);

  /* Clean up: */
  polytoop_delete(recpolytoop);
  free(points);

  return polytoop;
}



Polytoop* polytoop_frompoints(int npoints, int dim, double const* orgpoints)
{
  int ipoint;
  int idim;
  int* minindices;
  int* maxindices;
  double* minima;
  double* maxima;
  double* positions;
  Point* points;

  Polytoop* polytoop = polytoop_new();

  /* Initialize member variables: */
  polytoop->dim = dim;
  polytoop->isdelaunay = false;
  polytoop->shift = malloc(polytoop->dim * sizeof(double));
  polytoop->scales = malloc(polytoop->dim * sizeof(double));
  polytoop->center = malloc(polytoop->dim * sizeof(double));

  /* Bounding box: */
  minindices = alloca(polytoop->dim * sizeof(int));
  maxindices = alloca(polytoop->dim * sizeof(int));
  minima = alloca(polytoop->dim * sizeof(double));
  maxima = alloca(polytoop->dim * sizeof(double));
  boundingbox(
      npoints, polytoop->dim, orgpoints, minindices, maxindices, minima, maxima
  );

  /* Scales: */
  memcpy(polytoop->scales, maxima, polytoop->dim * sizeof(double));
  vec_sub(polytoop->dim, polytoop->scales, minima);

  /* Translation: */
  memcpy(polytoop->shift, minima, polytoop->dim * sizeof(double));

  /* Points array: */
  points = malloc(npoints * sizeof(Point));
  positions = malloc(npoints * dim * sizeof(double));
  for (ipoint = 0; ipoint < npoints; ++ipoint) {
    points[ipoint].next = NULL;
    points[ipoint].index = ipoint;
    points[ipoint].height = 0.0;
    points[ipoint].pos = positions + ipoint * dim;
    for (idim = 0; idim < polytoop->dim; ++idim) {
      points[ipoint].pos[idim] =
          (orgpoints[ipoint * polytoop->dim + idim] - polytoop->shift[idim]) /
          polytoop->scales[idim];
    }
  }

  /* Build the polytoop: */
  build(polytoop, npoints, points);

  /* Clean up: */
  free(positions);
  free(points);

  return polytoop;
}



Polytoop* polytoop_delaunay(int npoints, int dim, double const* orgpoints)
{
  int ipoint;
  int i;
  int j;
  int k;
  int* minindices;
  int* maxindices;
  double* positions;
  Point* points;

  Polytoop* polytoop = polytoop_new();
  Allocator* alc = &polytoop->alc;

  /* Reciprocal of dimension: */
  double recd = 1.0 / (double)dim;

  /* Initialize member variables: */
  polytoop->dim = dim + 1;
  polytoop->isdelaunay = true;
  polytoop->shift = malloc(polytoop->dim * sizeof(double));
  polytoop->scales = malloc(polytoop->dim * sizeof(double));
  polytoop->center = malloc(polytoop->dim * sizeof(double));

  /* Bounding box: */
  minindices = alloca(dim * sizeof(int));
  maxindices = alloca(dim * sizeof(int));
  boundingbox(
      npoints, dim, orgpoints, minindices, maxindices, polytoop->shift,
      polytoop->scales
  );
  polytoop->shift[dim] = 0.0;

  /* Scales: */
  vec_sub(dim, polytoop->scales, polytoop->shift);
  polytoop->scales[dim] = 1.0;

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
            (orgpoints[ipoint * dim + i] - polytoop->shift[i]) /
            polytoop->scales[i];
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

  /* Construct convex polytoop of paraboloid: */
  build(polytoop, npoints + 1, points);

  /* Clean up: */
  free(positions);
  free(points);

  /* Remove upper delaunay surfaces: */
  u32* hridges = malloc(polytoop->nridges * sizeof(u32));
  u32 hridge = polytoop->firstridge;
  for (i = 0; i < polytoop->nridges;
       ++i, hridge = ((Ridge*)allocator_mem(alc, hridge))->next) {
    hridges[i] = hridge;
  }

  for (i = 0; i < polytoop->nridges; ++i) {
    hridge = hridges[i];
    Ridge* ridge = allocator_mem(alc, hridge);
    for (j = 0; j < 2; ++j) {
      u32 hfacet = ridge->facets[j];
      if (hfacet == UINT32_MAX) {
        continue;
      }

      Facet* facet = allocator_mem(alc, hfacet);
      double* normal = facet->centroid + polytoop->dim;
      if (normal[polytoop->dim - 1] < 0.0) {
        continue;
      }

      /* Upper delaunay. */

      /* Disassociate facet from ridges: */
      ridge->facets[j] = UINT32_MAX;

      u32* ridges = (u32*)(facet->centroid + 2 * polytoop->dim);
      for (k = 0; k < polytoop->dim; ++k) {
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
      facet_remove(polytoop, hfacet);
      facet_free(polytoop, hfacet);
    }

    if (ridge->facets[0] == UINT32_MAX && ridge->facets[1] == UINT32_MAX) {
      /* Both adjacent facets are gone. */

      /* Disassociate from adjacent vertices: */
      for (j = 0; j < polytoop->dim - 1; ++j) {
        Vertex* vertex = allocator_mem(alc, ridge->vertices[j]);
        --vertex->nridges;
        if (vertex->nridges == 0) {
          /* Remove vertex as well: */
          vertex_remove(polytoop, ridge->vertices[j]);
        }
      }

      /* Remove ridge: */
      ridge_remove(polytoop, hridge);
      hridges[i] = UINT32_MAX; /* forget reference */
    }
  }
  free(hridges);

  return polytoop;
}



void polytoop_addvertex(Polytoop* polytoop, double const* point)
{
  int idim;
  Allocator* alc = &polytoop->alc;

  /* Allocations: */
  double* pos = alloca(polytoop->dim * sizeof(double));

  /* Initialize apex: */
  Point apex;
  apex.next = NULL;
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
  u32 hmaxfacet = UINT32_MAX;
  u32 hfacet = polytoop->firstfacet;
  while (hfacet != UINT32_MAX) {
    Facet* facet = allocator_mem(alc, hfacet);
    double* normal = facet->centroid + polytoop->dim;
    double h = vec_dot(polytoop->dim, pos, normal) - facet->dist;
    if (h > EPS) {
      hmaxfacet = hfacet;
      break;
    }
    hfacet = facet->next;
  }

  if (hmaxfacet != UINT32_MAX) {
    /* Add vertex to polytoop: */
    addpoint(polytoop, hmaxfacet, &apex);
  }
}



void polytoop_print(Polytoop* polytoop)
{
  Allocator* alc = &polytoop->alc;

  int d = polytoop->dim;
  double* sv = alloca(d * sizeof(double));
  double* normal = alloca(d * sizeof(double));
  double* centroid = alloca(d * sizeof(double));
  double* position = alloca(d * sizeof(double));

  printf("%d facets\n", polytoop->nfacets);
  printf("%d ridges\n", polytoop->nridges);
  printf("%d vertices\n", polytoop->nverts);

  vec_reset(d, sv);

  int i = 0;
  u32 hfacet = polytoop->firstfacet;
  while (hfacet != UINT32_MAX) {
    Facet* facet = allocator_mem(alc, hfacet);
    polytoop_facet_getnormal(facet, normal);
    double volume = polytoop_facet_getvolume(facet);
    vec_adds(d, sv, normal, volume);
    polytoop_facet_getcentroid(facet, centroid);

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
      polytoop_vertex_getposition(vertex, position);
      vec_print(d, position);
    }
    printf("\n");

    ++i;
    hfacet = facet->next;
  }
  assert(i == polytoop->nfacets);
  printf("  sv:\n");
  vec_print(d, sv);
  printf("\n");
}



int polytoop_getnumfacets(Polytoop* polytoop)
{
  return polytoop->nfacets;
}



int polytoop_getnumridges(Polytoop* polytoop)
{
  return polytoop->nridges;
}



int polytoop_getnumvertices(Polytoop* polytoop)
{
  return polytoop->nverts;
}



void polytoop_interpolate(
    Polytoop* polytoop,
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
  Allocator* alc = &polytoop->alc;

  assert(polytoop->isdelaunay);

  /* Some space for various tasks: */
  int d = polytoop->dim - 1;
  xiprime = alloca(d * sizeof(double));
  double* centroid = alloca(d * sizeof(double));
  verts = alloca(d * d * sizeof(double));

  /* Transform xi to local coordinates: */
  for (i = 0; i < d; ++i) {
    xiprime[i] = (xi[i] - polytoop->shift[i]) / polytoop->scales[i];
  }

  /* Start with first facet: */
  assert(polytoop->firstfacet != UINT32_MAX);
  u32 hcurrentfacet = polytoop->firstfacet;

  while (1) {
    Facet* currentfacet = allocator_mem(alc, hcurrentfacet);
    u32* ridges = (u32*)(currentfacet->centroid + 2 * polytoop->dim);
    u32* vertices = ridges + polytoop->dim;

    /* In the current facet, find the ridge where xi is highest above: */
    totalweight = 0.0;
    hmin = HUGE_VAL;
    minridge = NULL;
    for (iridge = 0; iridge < polytoop->dim; ++iridge) {
      /* Retrieve ridge: */
      u32 hridge = ridges[iridge];
      Ridge* ridge = allocator_mem(alc, hridge);
      u32 hvertex = vertices[iridge];
      Vertex* vertex = allocator_mem(alc, hvertex);

      /* On demand construction of ridge volume, distance and normal
       * (polytoop->dim - 1): */
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
    vec_scale(polytoop->dim, weights, 1.0 / totalweight);
  }
}



Facet* polytoop_firstfacet(Polytoop* polytoop)
{
  return allocator_mem(&polytoop->alc, polytoop->firstfacet);
}



Vertex* polytoop_firstvertex(Polytoop* polytoop)
{
  return allocator_mem(&polytoop->alc, polytoop->firstvertex);
}



Facet* polytoop_facet_nextfacet(Facet* facet)
{
  return allocator_mem(&facet->polytoop->alc, facet->next);
}



void polytoop_facet_getnormal(Facet* facet, double* normal)
{
  int d = facet->polytoop->dim;
  memcpy(normal, facet->centroid + d, d * sizeof(double));
  for (int i = 0; i < d; ++i) {
    normal[i] /= facet->polytoop->scales[i];
  }
  vec_normalize(d, normal);
}



void polytoop_facet_getcentroid(Facet* facet, double* centroid)
{
  int i;

  memcpy(centroid, facet->centroid, facet->polytoop->dim * sizeof(double));
  for (i = 0; i < facet->polytoop->dim; ++i) {
    centroid[i] *= facet->polytoop->scales[i];
  }
  vec_add(facet->polytoop->dim, centroid, facet->polytoop->shift);
}



double polytoop_facet_getoffset(Facet* facet)
{
  double* centroid = alloca(facet->polytoop->dim * sizeof(double));
  double* normal = alloca(facet->polytoop->dim * sizeof(double));
  polytoop_facet_getcentroid(facet, centroid);
  polytoop_facet_getnormal(facet, normal);
  return vec_dot(facet->polytoop->dim, centroid, normal);
}



double polytoop_facet_getvolume(Facet* facet)
{
  int d = facet->polytoop->dim;
  double* s = facet->polytoop->scales;
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



int polytoop_facet_getnumvertices(Facet* facet)
{
  return facet->polytoop->dim;
}



Vertex* polytoop_vertex_nextvertex(Vertex* vertex)
{
  return allocator_mem(&vertex->polytoop->alc, vertex->next);
}



void polytoop_vertex_getposition(Vertex* vertex, double* position)
{
  double* scales = vertex->polytoop->scales;
  double* shift = vertex->polytoop->shift;
  int i = vertex->polytoop->dim;
  i -= vertex->polytoop->isdelaunay ? 1 : 0;
  while (i--) {
    position[i] = vertex->position[i] * scales[i] + shift[i];
  }
}
