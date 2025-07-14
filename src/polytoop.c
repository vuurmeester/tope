#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "array.h"
#include "hashmap.h"
#include "math.h"
#include "util.h"

#include <polytoop.h>

#define EPS 1.0e-9



/* Create and initialize a vertex struct with a position and id: */
static polytoop_Vertex* vertex_new(Polytoop* polytoop, double* pos, int index)
{
  /* Allocate vertex: */
  polytoop_Vertex* vertex =
      allocator_alloc(polytoop->allocator, sizeof(polytoop_Vertex));

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
  /* Create ridge: */
  Ridge* ridge = allocator_alloc(polytoop->allocator, sizeof(Ridge));
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
  for (int i = 0; i < polytoop->dim - 1; ++i) {
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
  polytoop_Facet* facet =
      allocator_alloc(polytoop->allocator, sizeof(polytoop_Facet));
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
  facet->dist = 0.0;
  facet->ridges = array_new(polytoop->dim, polytoop->allocator);
  facet->vertices = array_new(polytoop->dim, polytoop->allocator);
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
  /* Facet centroid: */
  vec_reset(polytoop->dim, facet->centroid);
  for (int i = 0; i < polytoop->dim; ++i) {
    vec_add(polytoop->dim, facet->centroid,
            ((polytoop_Vertex*)facet->vertices.values[i])->position);
  }
  vec_scale(polytoop->dim, facet->centroid, 1.0 / (double)polytoop->dim);

  /* Vertex span: */
  double* dirs = alloca((polytoop->dim - 1) * polytoop->dim * sizeof(double));
  for (int i = 0; i < polytoop->dim - 1; ++i) {
    /* Vertex position difference: */
    memcpy(dirs + i * polytoop->dim,
           ((polytoop_Vertex*)facet->vertices.values[i + 1])->position,
           polytoop->dim * sizeof(double));
    vec_sub(polytoop->dim, dirs + i * polytoop->dim,
            ((polytoop_Vertex*)facet->vertices.values[0])->position);
  }

  /* Basis and volume: */
  facet->volume = 1.0;
  for (int i = 0; i < polytoop->dim - 1; ++i) {
    /* Pivot largest row on top: */
    int pivot = i;
    double maxnrmsq = vec_nrmsq(polytoop->dim, dirs + i * polytoop->dim);
    for (int j = i + 1; j < polytoop->dim - 1; ++j) {
      double nrmsq = vec_nrmsq(polytoop->dim, dirs + j * polytoop->dim);
      if (nrmsq > maxnrmsq) {
        maxnrmsq = nrmsq;
        pivot = j;
      }
    }
    if (pivot > i) {
      memswp(dirs + i * polytoop->dim, dirs + pivot * polytoop->dim,
             polytoop->dim * sizeof(double));
    }

    /* Normalize: */
    double norm = sqrt(maxnrmsq);
    vec_scale(polytoop->dim, dirs + i * polytoop->dim, 1.0 / norm);

    /* Orthogonalize: */
    for (int j = i + 1; j < polytoop->dim - 1; ++j) {
      double ip = vec_dot(polytoop->dim, dirs + j * polytoop->dim,
                          dirs + i * polytoop->dim);
      vec_adds(polytoop->dim, dirs + j * polytoop->dim,
               dirs + i * polytoop->dim, -ip);
    }

    /* Accumulate area: */
    facet->volume *= norm / (double)(i + 1);
  }

  /* Compute normal: */
  memcpy(facet->normal, facet->centroid, polytoop->dim * sizeof(double));
  vec_sub(polytoop->dim, facet->normal, polytoop->center);
  for (int i = 0; i < polytoop->dim - 1; ++i) {
    double ip = vec_dot(polytoop->dim, facet->normal, dirs + i * polytoop->dim);
    vec_adds(polytoop->dim, facet->normal, dirs + i * polytoop->dim, -ip);
  }
  vec_normalize(polytoop->dim, facet->normal);

  /* Plane distance: */
  facet->dist = vec_dot(polytoop->dim, facet->normal, facet->centroid);

  /* Assign ridges: */
  for (int i = 0; i < polytoop->dim; ++i) {
    /* Get ridge: */
    Ridge* ridge = ridges[i];

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



/* Create initial simplex, add outside vertices, etc: */
static void initialsimplex(Polytoop* polytoop, int npoints, Point* points)
{
  /* For readability: */
  int d = polytoop->dim;

  /* Allocations: */
  double* span = malloc((npoints - 1) * d * sizeof(double));
  double* vec = alloca(d * sizeof(double));

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

  /* Place extreme points up front: */
  memswp(&points[0], &points[maxindex], sizeof(Point*));

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

  /* Initialize polytoop center to initial simplex centroid: */
  vec_reset(polytoop->dim, polytoop->center);
  for (int i = 0; i < polytoop->dim + 1; ++i) {
    vec_add(polytoop->dim, polytoop->center, points[p[i]].pos);
  }
  vec_scale(polytoop->dim, polytoop->center, 1.0 / (double)(polytoop->dim + 1));

  /* Create initial simplex vertices: */
  polytoop_Vertex** vertices =
      alloca((polytoop->dim + 1) * sizeof(polytoop_Vertex*));
  for (int i = 0; i < polytoop->dim + 1; ++i) {
    vertices[i] = vertex_new(polytoop, points[p[i]].pos, points[p[i]].index);
  }

  Ridge** facetridges = alloca(polytoop->dim * sizeof(Ridge*));
  polytoop_Vertex** ridgeverts =
      alloca((polytoop->dim - 1) * sizeof(polytoop_Vertex*));

  hashmap_clear(&polytoop->newridges);

  /* Create facets: */
  for (int i = 0; i < polytoop->dim + 1; ++i) {
    /* Create facet: */
    polytoop_Facet* facet = facet_new(polytoop);

    /* Add all vertices except for one: */
    for (int j = 0; j < i; ++j) {
      array_append(&facet->vertices, vertices[j], polytoop->allocator);
    }
    for (int j = i + 1; j < polytoop->dim + 1; ++j) {
      array_append(&facet->vertices, vertices[j], polytoop->allocator);
    }

    /* Add ridges: */
    for (int j = 0; j < d; ++j) {
      /* Ridge verts: */
      for (int k = 0; k < j; ++k) {
        ridgeverts[k] = facet->vertices.values[k];
      }
      for (int k = j + 1; k < d; ++k) {
        ridgeverts[k - 1] = facet->vertices.values[k];
      }

      facetridges[j] = hashmap_retrieve(polytoop->newridges, d, ridgeverts);
      if (!facetridges[j]) {
        facetridges[j] = create_ridge(polytoop, ridgeverts);
        hashmap_insert(&polytoop->newridges, d, facetridges[j]);
      }
    }

    /* Vertex array filled, compute and integrate the facet: */
    addfacet(polytoop, facet, facetridges);
  }

  /* Create outside sets (assign each remaining vertex to a facet which it
   * 'sees'): */
  for (int i = polytoop->dim + 1; i < npoints; ++i) {
    /* Determine maximum distance: */
    for (polytoop_Facet* facet = polytoop->firstfacet; facet != NULL;
         facet = facet->next) {
      /* Distance to facet: */
      double h =
          vec_dot(polytoop->dim, facet->normal, points[p[i]].pos) - facet->dist;

      /* If above facet, add it to outside set and move to next point. */
      if (h > EPS) {
        points[p[i]].height = h;
        facet_addoutside(polytoop, facet, points + p[i]);
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
  int d = polytoop->dim;

  /* Allocations: */
  Array newfacets = array_new(8, polytoop->allocator);
  Array horizonridges = array_new(8, polytoop->allocator);
  double* vec = alloca(d * sizeof(double));

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
  polytoop_Facet* visiblelist = facet;
  facet->prev = NULL;
  facet->next = NULL;
  facet->visible = 1;

  /* Initialize outside points list: */
  Point* outsidepoints = NULL;

  while (visiblelist) {
    /* Pop a value: */
    facet = visiblelist;
    visiblelist = facet->next;

    /* Loop over ridges to visit neighbours: */
    for (int iridge = 0; iridge < facet->ridges.len; ++iridge) {
      /* Ridge i: */
      Ridge* ridge = facet->ridges.values[iridge];

      /* Retrieve neighbour: */
      polytoop_Facet* neighbour = NULL;
      if (facet == ridge->facets[0]) {
        /* Current facet is 0, neighbour is 1: */
        ridge->facets[0] = NULL; /* forget reference to this facet */
        neighbour = ridge->facets[1];
      } else {
        assert(facet == ridge->facets[1]);
        /* Current facet is 1, neighbour is 0: */
        ridge->facets[1] = NULL; /* forget reference to this facet */
        neighbour = ridge->facets[0];
      }

      if (!neighbour) {
        /* This neighbour was already removed, remove ridge. */

        /* Remove reference to this ridge from adjacent vertices: */
        for (int ivertex = 0; ivertex < polytoop->dim - 1; ++ivertex) {
          polytoop_Vertex* vertex = ridge->vertices[ivertex];
          --vertex->nridges;
          if (vertex->nridges == 0) {
            /* Remove the vertex as well: */
            vertex_remove(polytoop, vertex);
          }
        }

        /* Remove the ridge from the polytoop: */
        ridge_remove(polytoop, ridge);
      } else if (!neighbour->visible) {
        /* Neighbour is untested. */

        /* Height of apex above neighbour: */
        double h = vec_dot(d, apex->pos, neighbour->normal) - neighbour->dist;

        /* Decide if neighbour is visible: */
        if (h > EPS) {
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
          neighbour->prev = NULL; /* singly linked */
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
  polytoop_Vertex* vertex = vertex_new(polytoop, apex->pos, apex->index);

  polytoop_Vertex** ridgeverts =
      alloca((polytoop->dim - 1) * sizeof(polytoop_Vertex*));
  Ridge** facetridges = alloca(polytoop->dim * sizeof(Ridge*));

  hashmap_clear(&polytoop->newridges);

  /* Form new facets: */
  while (horizonridges.len > 0) {
    /* The horizon ridge: */
    Ridge* horizonridge = array_pop(&horizonridges);

    /* New facet: */
    facet = facet_new(polytoop);

    /* Facet vertices are apex + horizon vertices: */
    facet->vertices.len = d;
    assert(facet->vertices.len <= facet->vertices.cap);
    facet->vertices.values[0] = vertex;
    for (int ivertex = 1; ivertex < d; ++ivertex) {
      facet->vertices.values[ivertex] = horizonridge->vertices[ivertex - 1];
    }

    /* 1 existing (horizon) ridge and d - 1 new ridges: */
    facetridges[0] = horizonridge;
    for (int iridge = 1; iridge < d; ++iridge) {
      /* Add all facet vertices except the one at index 'iridge': */
      for (int ivertex = 0; ivertex < iridge; ++ivertex) {
        ridgeverts[ivertex] = facet->vertices.values[ivertex];
      }
      for (int ivertex = iridge + 1; ivertex < d; ++ivertex) {
        ridgeverts[ivertex - 1] = facet->vertices.values[ivertex];
      }

      /* Get or create new ridge: */
      Ridge* newridge = hashmap_retrieve(polytoop->newridges, d, ridgeverts);
      if (!newridge) {
        /* Create new ridge: */
        newridge = create_ridge(polytoop, ridgeverts);
        hashmap_insert(&polytoop->newridges, d, newridge);
      }
      facetridges[iridge] = newridge;
    }

    addfacet(polytoop, facet, facetridges);

    if (polytoop->merge) {
      /* Retrieve neighbour: */
      polytoop_Facet* neighbour = NULL;
      if (horizonridge->facets[0] == facet) {
        neighbour = horizonridge->facets[1];
      } else {
        neighbour = horizonridge->facets[0];
      }

      /* Normal vector difference: */
      memcpy(vec, facet->normal, polytoop->dim * sizeof(double));
      vec_sub(polytoop->dim, vec, neighbour->normal);

      if (vec_nrmsq(polytoop->dim, vec) < EPS * EPS) {
        /* Parallel, facet is to be removed. */

        /* Add facet ridges to neighbour: */
        for (int iridge2 = 0; iridge2 < facet->ridges.len; ++iridge2) {
          Ridge* ridge2 = facet->ridges.values[iridge2];
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
        double newvol = facet->volume + neighbour->volume;
        vec_scale(polytoop->dim, neighbour->centroid,
                  neighbour->volume / newvol);
        vec_adds(polytoop->dim, neighbour->centroid, facet->centroid,
                 facet->volume / newvol);
        neighbour->volume = newvol;

        /* Disassociate the dividing ridge from the neighbour: */
        array_remove(&neighbour->ridges,
                     array_find(neighbour->ridges, horizonridge));

        /* Disassociate vertices from the dividing ridge: */
        for (int ivertex = 0; ivertex < polytoop->dim - 1; ++ivertex) {
          --horizonridge->vertices[ivertex]->nridges;
          assert(horizonridge->vertices[ivertex]->nridges > 0);
        }

        /* Remove the dividing ridge: */
        ridge_remove(polytoop, horizonridge);

        /* Remove the facet and replace by neighbour: */
        facet_remove(polytoop, facet);
        facet = neighbour; /* so it is added to newfacets */
      }
    }

    /* Add new facet to newfacets array: */
    array_append(&newfacets, facet, polytoop->allocator);
  }

  /* Assign outside verts: */
  facet = newfacets.values[0];
  while (outsidepoints) {
    /* Remember next in list: */
    Point* next = outsidepoints->next;

    double h =
        vec_dot(polytoop->dim, outsidepoints->pos, facet->normal) - facet->dist;

    /* Visit neighbours, skipping horizon ridge: */
    polytoop_Facet* prev = NULL;
    for (int iridge = 1; iridge < facet->ridges.len; ++iridge) {
      Ridge* ridge = facet->ridges.values[iridge];
      assert(iridge > 0 && array_find(horizonridges, ridge) == -1);
      polytoop_Facet* neighbour = NULL;
      if (ridge->facets[0] == facet) {
        neighbour = ridge->facets[1];
      } else {
        assert(ridge->facets[1] == facet);
        neighbour = ridge->facets[0];
      }
      if (neighbour == prev) {
        continue;
      }

      /* Neighbour height: */
      double nh =
          vec_dot(polytoop->dim, outsidepoints->pos, neighbour->normal) -
          neighbour->dist;
      if (nh > h) {
        iridge = 0; /* reset iteration (new base facet) */
        prev = facet;
        facet = neighbour;
        h = nh;
      }
    }

    if (h > EPS) {
      /* Outside. */
      outsidepoints->height = h;
      facet_addoutside(polytoop, facet, outsidepoints);
    }

    /* Next in list: */
    outsidepoints = next;
  }

  array_delete(&horizonridges, polytoop->allocator);
  array_delete(&newfacets, polytoop->allocator);
}



static void build(Polytoop* polytoop, int npoints, Point* points)
{
  assert(npoints >= polytoop->dim + 1);

  initialsimplex(polytoop, npoints, points);

  /* Now, all vertices are either on the initial simplex OR
     added to the outside set of one of its facets OR
     inside the initial simplex. */
  while (polytoop->firstfacet->outsidehead) {
    /* Pop head: */
    Point* point = polytoop->firstfacet->outsidehead;
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

  hashmap_init(&polytoop->newridges);

  return polytoop;
}



void polytoop_delete(Polytoop* polytoop)
{
  polytoop_clear(polytoop);
  hashmap_destroy(&polytoop->newridges);
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

  hashmap_clear(&polytoop->newridges);
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
    double dist = dists[i] - vec_dot(d, x, &normals[i * d]);
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
    /* Get vertices: */
    for (ivert = 0; ivert < d; ++ivert) {
      polytoop_vertex_getposition(polytoop_facet_getvertex(facet, ivert),
                                  &verts[ivert * d]);
    }

    /* Simplex analysis (replaces verts with span): */
    double vol;
    analysesimplex(d, d, verts, &vol, centroid);

    /* Construct normal: */
    memcpy(normal, centroid, d * sizeof(double));
    for (i = 0; i < d - 1; ++i) {
      double fac = vec_dot(d, normal, &verts[i * d]);
      vec_adds(d, normal, &verts[i * d], -fac);
    }
    double dist2 = vec_nrmsq(d, normal);

    /* Construct vertex position: */
    memcpy(&points[ifacet * d], normal, d * sizeof(double));
    vec_scale(d, &points[ifacet * d], 1.0 / dist2);
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
  int* minindices = alloca(polytoop->dim * sizeof(int));
  int* maxindices = alloca(polytoop->dim * sizeof(int));
  double* maxima = alloca(polytoop->dim * sizeof(double));
  double* minima = alloca(polytoop->dim * sizeof(double));
  boundingbox(npoints, polytoop->dim, orgpoints, minindices, maxindices, minima,
              maxima);

  /* Scales: */
  memcpy(polytoop->scales, maxima, polytoop->dim * sizeof(double));
  vec_sub(polytoop->dim, polytoop->scales, minima);

  /* Translation: */
  memcpy(polytoop->shift, minima, polytoop->dim * sizeof(double));

  /* Points array: */
  Point* points = malloc(npoints * sizeof(Point));
  double* positions = malloc(npoints * dim * sizeof(double));
  for (int ipoint = 0; ipoint < npoints; ++ipoint) {
    points[ipoint].next = NULL;
    points[ipoint].index = ipoint;
    points[ipoint].height = 0.0;
    points[ipoint].pos = positions + ipoint * dim;
    for (int idim = 0; idim < polytoop->dim; ++idim) {
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
  Point* points = malloc((npoints + 1) * sizeof(Point));
  double* positions = malloc((npoints + 1) * (dim + 1) * sizeof(double));
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
  /* Allocations: */
  double* pos = alloca(polytoop->dim * sizeof(double));

  /* Initialize apex: */
  Point apex = {
      .next = NULL, .index = polytoop->nverts, .height = 0.0, .pos = pos};

  /* Transformed position: */
  memcpy(pos, point, polytoop->dim * sizeof(double));
  vec_sub(polytoop->dim, pos, polytoop->shift);
  for (int idim = 0; idim < polytoop->dim; ++idim) {
    if (polytoop->scales[idim] > 0.0) {
      pos[idim] /= polytoop->scales[idim];
    }
  }

  /* Find facet that can see the point: */
  polytoop_Facet* maxfacet = NULL;
  for (polytoop_Facet* facet = polytoop->firstfacet; facet != NULL;
       facet = facet->next) {
    double h = vec_dot(polytoop->dim, pos, facet->normal) - facet->dist;
    if (h > EPS) {
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
  printf("%d facets\n", polytoop->nfacets);
  printf("%d ridges\n", polytoop->nridges);
  printf("%d vertices\n", polytoop->nverts);

  int d = polytoop->dim;

  double* sv = alloca(d * sizeof(double));
  vec_reset(d, sv);

  int i;
  polytoop_Facet* facet;
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
    for (int j = 0; j < facet->vertices.len; ++j) {
      polytoop_Vertex* vertex = facet->vertices.values[j];
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
  assert(polytoop->isdelaunay);

  /* Some space for various tasks: */
  int dim = polytoop->dim - 1;
  double* xiprime = alloca(dim * sizeof(double));
  double* verts = alloca(dim * dim * sizeof(double));

  /* Transform xi to local coordinates: */
  for (int i = 0; i < dim; ++i) {
    xiprime[i] = (xi[i] - polytoop->shift[i]) / polytoop->scales[i];
  }

  /* Start with first facet: */
  assert(polytoop->firstfacet != NULL);
  polytoop_Facet* currentfacet = polytoop->firstfacet;

  double totalweight;
  while (1) {
    /* In the current facet, find the ridge where xi is highest above: */
    totalweight = 0.0;
    double hmin = HUGE_VAL;
    int minindex = -1;
    for (int iridge = 0; iridge < polytoop->dim; ++iridge) {
      /* Retrieve ridge: */
      Ridge* ridge = currentfacet->ridges.values[iridge];

      /* On demand construction of ridge normal and centroid (d - 1): */
      if (ridge->normal == NULL) {
        /* Matrix of vertex coordinates: */
        for (int ivertex = 0; ivertex < dim; ++ivertex) {
          memcpy(&verts[ivertex * dim], ridge->vertices[ivertex]->position,
                 dim * sizeof(double));
        }

        /* Analyse ridge simplex: */
        ridge->centroid =
            allocator_alloc(polytoop->allocator, dim * sizeof(double));
        analysesimplex(dim, dim, verts, &ridge->volume, ridge->centroid);

        /* Construct normal: */
        ridge->normal =
            allocator_alloc(polytoop->allocator, dim * sizeof(double));
        memcpy(ridge->normal, ridge->centroid, dim * sizeof(double));
        vec_sub(dim, ridge->normal, currentfacet->centroid);
        for (int i = 0; i < dim - 1; ++i) {
          double fac = vec_dot(dim, ridge->normal, &verts[i * dim]);
          vec_adds(dim, ridge->normal, &verts[i * dim], -fac);
        }
        vec_normalize(dim, ridge->normal);
      }

      /* Ridge plane distance: */
      double dist = vec_dot(dim, ridge->normal, ridge->centroid);

      /* Sign of normal: */
      double sign = 1.0;
      if (vec_dot(dim, currentfacet->centroid, ridge->normal) < dist) {
        /* outward pointing normal */
        sign = -1.0;
      }

      /* Height of interpolation point above ridge: */
      double h = sign * (vec_dot(dim, xiprime, ridge->normal) - dist);

      /* Weight for this ridge: */
      weights[iridge] = h * ridge->volume;
      totalweight += weights[iridge];

      /* Index for this ridge: */
      polytoop_Vertex* vertex = currentfacet->vertices.values[iridge];
      indices[iridge] = vertex->index;
      assert(find(dim, ridge->vertices, vertex) == -1);

      /* Keep track of minimum height (if not boundary ridge): */
      if (h < hmin && ridge->facets[0] != NULL && ridge->facets[1] != NULL) {
        hmin = h;
        minindex = iridge;
      }
    }

    /* If hmin insignificant, stop: */
    if (hmin > -EPS) {
      break;
    }

    /* Else, goto next facet (cross ridge): */
    Ridge* ridge = currentfacet->ridges.values[minindex];
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
  assert(i < facet->ridges.len);
  Ridge* ridge = facet->ridges.values[i];
  if (ridge->facets[0] == facet) {
    return ridge->facets[1];
  }
  assert(ridge->facets[1] == facet);
  return ridge->facets[0];
}



void polytoop_facet_getnormal(polytoop_Facet* facet, double* normal)
{
  memcpy(normal, facet->normal, facet->polytoop->dim * sizeof(double));
  for (int i = 0; i < facet->polytoop->dim; ++i) {
    normal[i] /= facet->polytoop->scales[i];
  }
  vec_normalize(facet->polytoop->dim, normal);
}



void polytoop_facet_getcentroid(polytoop_Facet* facet, double* centroid)
{
  memcpy(centroid, facet->centroid, facet->polytoop->dim * sizeof(double));
  for (int i = 0; i < facet->polytoop->dim; ++i) {
    centroid[i] *= facet->polytoop->scales[i];
  }
  vec_add(facet->polytoop->dim, centroid, facet->polytoop->shift);
}



double polytoop_facet_getoffset(polytoop_Facet* facet)
{
  double* centroid = alloca(facet->polytoop->dim * sizeof(double));
  double* normal = alloca(facet->polytoop->dim * sizeof(double));
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
  int dim = vertex->polytoop->dim;
  if (vertex->polytoop->isdelaunay) {
    --dim;
  }

  memcpy(position, vertex->position, dim * sizeof(double));
  for (int i = 0; i < dim; ++i) {
    position[i] *= vertex->polytoop->scales[i];
  }
  vec_add(dim, position, vertex->polytoop->shift);
}
