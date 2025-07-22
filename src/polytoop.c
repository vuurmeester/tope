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
static Vertex* vertex_new(Polytoop* polytoop, double* pos, int index)
{
  /* Allocate vertex: */
  Vertex* vertex = allocator_alloc(
                     &polytoop->alc, sizeof(Vertex) + (polytoop->dim - 1) * sizeof(double));

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
  memcpy(vertex->position, pos, polytoop->dim * sizeof(double));
  vertex->nridges = 0;

  return vertex;
}



static void vertex_free(Polytoop* polytoop, Vertex* vertex)
{
  allocator_free(&polytoop->alc, vertex,
                 sizeof(Vertex) + (polytoop->dim - 1) * sizeof(double));
}



static void vertex_remove(Polytoop* polytoop, Vertex* vertex)
{
  Vertex* next = vertex->next;
  Vertex* prev = vertex->prev;
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
static Ridge* create_ridge(Polytoop* polytoop, Vertex** vertices)
{
  /* Create ridge: */
  Ridge* ridge = allocator_alloc(
                   &polytoop->alc, sizeof(Ridge) + (polytoop->dim - 2) * sizeof(Vertex*));
  ridge->next = NULL;
  ridge->prev = NULL;
  ridge->volume = 0.0;
  ridge->centroid = NULL;
  ridge->normal = NULL;
  ridge->facets[0] = NULL;
  ridge->facets[1] = NULL;
  memcpy(ridge->vertices, vertices, (polytoop->dim - 1) * sizeof(Vertex*));

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



static void ridge_free(Polytoop* polytoop, Ridge* ridge)
{
  if (polytoop->isdelaunay) {
    allocator_free(&polytoop->alc, ridge->normal,
                   (polytoop->dim - 1) * sizeof(double));
    allocator_free(&polytoop->alc, ridge->centroid,
                   (polytoop->dim - 1) * sizeof(double));
  }
  allocator_free(&polytoop->alc, ridge,
                 sizeof(Ridge) + (polytoop->dim - 2) * sizeof(Vertex*));
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



/* Create and initialize a facet struct: */
static Facet* facet_new(Polytoop* polytoop)
{
  int d = polytoop->dim;

  Facet* facet = allocator_alloc(
                   &polytoop->alc, sizeof(Facet) + d * (2 * sizeof(double) + sizeof(Ridge*) +
                                                        sizeof(Vertex*)));
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
  facet->centroid = (double*)(facet + 1);
  vec_reset(d, facet->centroid);
  facet->normal = facet->centroid + d;
  vec_reset(d, facet->normal);
  facet->dist = 0.0;
  facet->ridges = (Ridge**)(facet->normal + d);
  memset(facet->ridges, 0, d * sizeof(Ridge*));
  facet->vertices = (Vertex**)(facet->ridges + d);
  memset(facet->vertices, 0, d * sizeof(Vertex*));
  facet->outsidehead = NULL;
  facet->outsidetail = NULL;
  facet->visible = 0;

  return facet;
}



static void facet_free(Polytoop* polytoop, Facet* facet)
{
  allocator_free(&polytoop->alc, facet,
                 sizeof(Facet) +
                 polytoop->dim * (2 * sizeof(double) + sizeof(Ridge*) +
                                  sizeof(Vertex*)));
}



static void facet_remove(Polytoop* polytoop, Facet* facet)
{
  Facet* next = facet->next;
  Facet* prev = facet->prev;

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



static void facet_addoutside(Polytoop* polytoop, Facet* facet, Point* point)
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



/* Integrate facet into polytoop (compute volume, create ridges, etc): */
static void addfacet(Polytoop* polytoop, Facet* facet, Ridge** ridges)
{
  int d = polytoop->dim;

  /* Facet centroid: */
  vec_reset(d, facet->centroid);
  for (int i = 0; i < d; ++i) {
    vec_add(d, facet->centroid, facet->vertices[i]->position);
  }
  vec_scale(d, facet->centroid, 1.0 / (double)d);

  /* Vertex span: */
  double* dirs = alloca((d - 1) * d * sizeof(double));
  double* x0 = facet->vertices[0]->position;
  for (int i = 0; i < d - 1; ++i) {
    /* Vertex position difference: */
    Vertex* vertex = facet->vertices[i + 1];
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
    double norm = sqrt(maxnrmsq);
    vec_scale(d, dirs + i * d, 1.0 / norm);

    /* Orthogonalize: */
    for (int j = i + 1; j < d - 1; ++j) {
      double ip = vec_dot(d, dirs + j * d, dirs + i * d);
      vec_adds(d, dirs + j * d, dirs + i * d, -ip);
    }

    /* Accumulate area: */
    facet->volume *= norm / (double)(i + 1);
  }

  /* Compute normal: */
  memcpy(facet->normal, facet->centroid, d * sizeof(double));
  vec_sub(d, facet->normal, polytoop->center);
  for (int i = 0; i < d - 1; ++i) {
    double ip = vec_dot(d, facet->normal, dirs + i * d);
    vec_adds(d, facet->normal, dirs + i * d, -ip);
  }
  vec_normalize(d, facet->normal);

  /* Plane distance: */
  facet->dist = vec_dot(d, facet->normal, facet->centroid);

  /* Assign ridges: */
  for (int i = 0; i < d; ++i) {
    /* Get ridge: */
    Ridge* ridge = ridges[i];

    /* Remember facet. */
    if (ridge->facets[0] == NULL) {
      ridge->facets[0] = facet;
    } else {
      assert(ridge->facets[1] == NULL);
      ridge->facets[1] = facet;
    }

    facet->ridges[i] = ridge;
  }
}



/* Create initial simplex, add outside vertices, etc: */
static void initialsimplex(Polytoop* polytoop, int npoints, Point* points)
{
  /* For readability: */
  int d = polytoop->dim;

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
  int* p = alloca(npoints * sizeof(int));
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
  Vertex** vertices = alloca((polytoop->dim + 1) * sizeof(Vertex*));
  for (int i = 0; i < polytoop->dim + 1; ++i) {
    vertices[i] = vertex_new(polytoop, points[p[i]].pos, points[p[i]].index);
  }

  Ridge** facetridges = alloca(polytoop->dim * sizeof(Ridge*));
  Vertex** ridgeverts = alloca((polytoop->dim - 1) * sizeof(Vertex*));

  hashmap_clear(&polytoop->newridges);

  /* Create facets: */
  for (int i = 0; i < polytoop->dim + 1; ++i) {
    /* Create facet: */
    Facet* facet = facet_new(polytoop);

    /* Add all vertices except for one: */
    for (int j = 0; j < i; ++j) {
      facet->vertices[j] = vertices[j];
    }
    for (int j = i + 1; j < polytoop->dim + 1; ++j) {
      facet->vertices[j - 1] = vertices[j];
    }

    /* Add ridges: */
    for (int j = 0; j < d; ++j) {
      /* Ridge verts: */
      for (int k = 0; k < j; ++k) {
        ridgeverts[k] = facet->vertices[k];
      }
      for (int k = j + 1; k < d; ++k) {
        ridgeverts[k - 1] = facet->vertices[k];
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
    for (Facet* facet = polytoop->firstfacet; facet != NULL;
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
  free(span);
}



/* Add point to polytoop: */
static void addpoint(Polytoop* polytoop, Facet* facet, Point* apex)
{
  int d = polytoop->dim;

  /* Allocations: */
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
  Facet* visiblelist = facet;
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
    for (int i = 0; i < d; ++i) {
      /* Ridge i: */
      Ridge* ridge = facet->ridges[i];

      /* Retrieve neighbour: */
      Facet* neighbour = NULL;
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
          Vertex* vertex = ridge->vertices[ivertex];
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
          if (polytoop->horizonridges_len == polytoop->horizonridges_cap) {
            polytoop->horizonridges_cap =
              3 * polytoop->horizonridges_cap / 2 + 1;
            polytoop->horizonridges =
              realloc(polytoop->horizonridges,
                      polytoop->horizonridges_cap * sizeof(Ridge*));
          }
          polytoop->horizonridges[polytoop->horizonridges_len++] = ridge;
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
  Vertex* vertex = vertex_new(polytoop, apex->pos, apex->index);

  Vertex** ridgeverts = alloca((polytoop->dim - 1) * sizeof(Vertex*));
  Ridge** facetridges = alloca(polytoop->dim * sizeof(Ridge*));

  hashmap_clear(&polytoop->newridges);

  /* Form new facets: */
  while (polytoop->horizonridges_len > 0) {
    /* Pop a horizon ridge: */
    Ridge* horizonridge =
      polytoop->horizonridges[--polytoop->horizonridges_len];

    /* New facet: */
    facet = facet_new(polytoop);

    /* Facet vertices are apex + horizon vertices: */
    facet->vertices[0] = vertex;
    for (int ivertex = 1; ivertex < d; ++ivertex) {
      facet->vertices[ivertex] = horizonridge->vertices[ivertex - 1];
    }

    facetridges[0] = horizonridge;
    for (int iridge = 1; iridge < d; ++iridge) {
      /* Add all facet vertices except the one at index 'iridge': */
      for (int ivertex = 0; ivertex < iridge; ++ivertex) {
        ridgeverts[ivertex] = facet->vertices[ivertex];
      }
      for (int ivertex = iridge + 1; ivertex < d; ++ivertex) {
        ridgeverts[ivertex - 1] = facet->vertices[ivertex];
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

    /* Add new facet to newfacets array: */
    if (polytoop->newfacets_len == polytoop->newfacets_cap) {
      polytoop->newfacets_cap = 3 * polytoop->newfacets_cap / 2 + 1;
      polytoop->newfacets = realloc(polytoop->newfacets,
                                    polytoop->newfacets_cap * sizeof(Facet*));
    }
    polytoop->newfacets[polytoop->newfacets_len++] = facet;
  }

  /* Assign outside verts: */
  facet = polytoop->newfacets[0];
  while (outsidepoints) {
    /* Remember next in list: */
    Point* next = outsidepoints->next;

    double h =
      vec_dot(polytoop->dim, outsidepoints->pos, facet->normal) - facet->dist;

    /* Visit neighbours, skipping horizon ridge: */
    Facet* prev = NULL;
    for (int iridge = 1; iridge < d; ++iridge) {
      Ridge* ridge = facet->ridges[iridge];
      Facet* neighbour = NULL;
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
        prev = facet;
        facet = neighbour;
        h = nh;
        iridge = 0; /* reset iteration (new base facet) */
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

  assert(polytoop->horizonridges_len == 0);
  polytoop->newfacets_len = 0;
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



static Polytoop* polytoop_new()
{
  Polytoop* polytoop = calloc(1, sizeof(Polytoop));

  allocator_init(&polytoop->alc);
  hashmap_init(&polytoop->newridges);

  return polytoop;
}



void polytoop_delete(Polytoop* polytoop)
{
  free(polytoop->newfacets);
  free(polytoop->horizonridges);
  hashmap_destroy(&polytoop->newridges);
  allocator_clear(&polytoop->alc);
  free(polytoop);
}



Polytoop* polytoop_fromplanes(int n, int d, double* normals, double* dists)
{
  /* Allocations: */
  double* mata = malloc(n * (d + 1) * sizeof(double));
  double* b = alloca(n * sizeof(double));
  double* c = alloca((d + 1) * sizeof(double));
  double* x = alloca((d + 1) * sizeof(double));

  /* Copy data to augmented matrix: */
  for (int i = 0; i < n; ++i) {
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
  free(mata);
  if (x[d] <= 0.0) {
    /* Infeasible. */
    return NULL;
  }

  /* Construct reciprocal points: */
  double* points = malloc(n * d * sizeof(double));
  for (int i = 0; i < n; ++i) {
    memcpy(&points[i * d], &normals[i * d], d * sizeof(double));
    double dist = dists[i] - vec_dot(d, x, &normals[i * d]);
    assert(dist > 0.0);
    vec_scale(d, &points[i * d], 1.0 / dist);
  }

  /* Construct reciprocal polytoop: */
  Polytoop* recpolytoop = polytoop_frompoints(n, d, points);
  free(points);

  /* Convert reciprocal polytoop to polytoop: */
  int nfacets = polytoop_getnumfacets(recpolytoop);
  points = malloc(nfacets * d * sizeof(double));

  int ifacet;
  Facet* facet;
  for (facet = recpolytoop->firstfacet, ifacet = 0; facet != NULL;
       facet = facet->next, ++ifacet) {
    /* Unscaled normal and distance: */
    double dist = facet->dist;
    for (int i = 0; i < d; ++i) {
      points[ifacet * d + i] = facet->normal[i] / recpolytoop->scales[i];
      dist += points[ifacet * d + i] * recpolytoop->shift[i];
    }

    /* Straight space point: */
    vec_scale(d, points + ifacet * d, 1.0 / dist);
    vec_add(d, points + ifacet * d, x);
  }

  Polytoop* polytoop = polytoop_frompoints(nfacets, d, points);

  /* Clean up: */
  polytoop_delete(recpolytoop);
  free(points);

  return polytoop;
}



Polytoop* polytoop_frompoints(int npoints, int dim, double* orgpoints)
{
  Polytoop* polytoop = polytoop_new();

  /* Initialize member variables: */
  polytoop->dim = dim;
  polytoop->isdelaunay = 0;
  polytoop->shift =
    allocator_alloc(&polytoop->alc, polytoop->dim * sizeof(double));
  polytoop->scales =
    allocator_alloc(&polytoop->alc, polytoop->dim * sizeof(double));
  polytoop->center =
    allocator_alloc(&polytoop->alc, polytoop->dim * sizeof(double));

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

  return polytoop;
}



Polytoop* polytoop_delaunay(int npoints, int dim, double* orgpoints)
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
  Facet* facet;

  Polytoop* polytoop = polytoop_new();

  /* Reciprocal of dimension: */
  recd = 1.0 / (double)dim;

  /* Initialize member variables: */
  polytoop->dim = dim + 1;
  polytoop->isdelaunay = 1;
  polytoop->shift =
    allocator_alloc(&polytoop->alc, polytoop->dim * sizeof(double));
  polytoop->scales =
    allocator_alloc(&polytoop->alc, polytoop->dim * sizeof(double));
  polytoop->center =
    allocator_alloc(&polytoop->alc, polytoop->dim * sizeof(double));

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

        for (int k = 0; k < polytoop->dim; ++k) {
          Ridge* ridgek = facet->ridges[k];
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

  return polytoop;
}



void polytoop_addvertex(Polytoop* polytoop, double* point)
{
  /* Allocations: */
  double* pos = alloca(polytoop->dim * sizeof(double));

  /* Initialize apex: */
  Point apex = {
    .next = NULL, .index = polytoop->nverts, .height = 0.0, .pos = pos
  };

  /* Transformed position: */
  memcpy(pos, point, polytoop->dim * sizeof(double));
  vec_sub(polytoop->dim, pos, polytoop->shift);
  for (int idim = 0; idim < polytoop->dim; ++idim) {
    if (polytoop->scales[idim] > 0.0) {
      pos[idim] /= polytoop->scales[idim];
    }
  }

  /* Find facet that can see the point: */
  Facet* maxfacet = NULL;
  for (Facet* facet = polytoop->firstfacet; facet != NULL;
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
  Facet* facet;
  for (facet = polytoop->firstfacet, i = 0; facet != NULL;
       facet = facet->next, ++i) {

    vec_adds(d, sv, facet->normal, facet->volume);

    printf("facet %d\n", i + 1);
    printf("  volume: %g\n", facet->volume);
    printf("  centroid: ");
    vec_print(d, facet->centroid);
    printf("\n");
    printf("  normal: ");
    vec_print(d, facet->normal);
    printf("\n");

    printf("  vertices:\n");
    for (int i = 0; i < d; ++i) {
      Vertex* vertex = facet->vertices[i];
      vec_print(d, vertex->position);
    }
    printf("\n");
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
  Facet* currentfacet = polytoop->firstfacet;

  double totalweight;
  while (1) {
    /* In the current facet, find the ridge where xi is highest above: */
    totalweight = 0.0;
    double hmin = HUGE_VAL;
    Ridge* minridge = NULL;
    for (int iridge = 0; iridge < polytoop->dim; ++iridge) {
      /* Retrieve ridge: */
      Ridge* ridge = currentfacet->ridges[iridge];
      Vertex* vertex = currentfacet->vertices[iridge];

      /* On demand construction of ridge normal and centroid (d - 1): */
      if (ridge->normal == NULL) {
        /* Matrix of vertex coordinates: */
        for (int ivertex = 0; ivertex < dim; ++ivertex) {
          memcpy(&verts[ivertex * dim], ridge->vertices[ivertex]->position,
                 dim * sizeof(double));
        }

        /* Analyse ridge simplex: */
        ridge->centroid = allocator_alloc(&polytoop->alc, dim * sizeof(double));
        analysesimplex(dim, dim, verts, &ridge->volume, ridge->centroid);

        /* Construct normal: */
        ridge->normal = allocator_alloc(&polytoop->alc, dim * sizeof(double));
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
    Ridge* ridge = minridge;
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



Facet* polytoop_firstfacet(Polytoop* polytoop)
{
  return polytoop->firstfacet;
}



Vertex* polytoop_firstvertex(Polytoop* polytoop)
{
  return polytoop->firstvertex;
}



Facet* polytoop_facet_nextfacet(Facet* facet)
{
  return facet->next;
}



void polytoop_facet_getnormal(Facet* facet, double* normal)
{
  memcpy(normal, facet->normal, facet->polytoop->dim * sizeof(double));
  for (int i = 0; i < facet->polytoop->dim; ++i) {
    normal[i] /= facet->polytoop->scales[i];
  }
  vec_normalize(facet->polytoop->dim, normal);
}



void polytoop_facet_getcentroid(Facet* facet, double* centroid)
{
  memcpy(centroid, facet->centroid, facet->polytoop->dim * sizeof(double));
  for (int i = 0; i < facet->polytoop->dim; ++i) {
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
  return facet->volume;
}



int polytoop_facet_getnumvertices(Facet* facet)
{
  return facet->polytoop->dim;
}



Vertex* polytoop_vertex_nextvertex(Vertex* vertex)
{
  return vertex->next;
}



void polytoop_vertex_getposition(Vertex* vertex, double* position)
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
