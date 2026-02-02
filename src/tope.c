#include <assert.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "hashmap.h"
#include "math.h"
#include "util.h"
#include "list.h"

#include <tope.h>

#define EPS 1e-9



static double height(
  int d,
  double const* x,  // point
  double const* c,  // ref
  double const* n   // normal
)
{
  /* Perpendicular distance of x from c on plane with normal n. */
  double sum = 0.0;
  for (int i = 0; i < d; ++i) {
    sum += (x[i] - c[i]) * n[i];
  }
  return sum;
}



/* Create and initialize a vertex struct with a position and id: */
static Vertex* vertex_create(Tope* tope, double const* pos, int index)
{
  ++tope->nverts;
  int d = tope->dim;
  Allocator* alc = &tope->alc;

  /* Allocate vertex: */
  Vertex* vertex = allocator_alloc(alc, sizeof(Vertex) + (d - 1) * sizeof(double));

  /* Fill in other variables: */
  vertex->index = index;
  memcpy(vertex->position, pos, d * sizeof(double));

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



/* Create and initialize a simplicial ridge with d - 1 vertices: */
static Ridge* ridge_create(Tope* tope, Vertex** verts)
{
  ++tope->nridges;
  Allocator* alc = &tope->alc;
  int d = tope->dim;

  /* Allocate ridge: */
  Ridge* ridge = allocator_alloc(alc, sizeof(Ridge) + (d - 2) * sizeof(Vertex*));

  // Assign and associate verts:
  for (int i = 0; i < d - 1; ++i) {
    ++verts[i]->nridges;
    ridge->verts[i] = verts[i];
  }

  return ridge;
}



static void ridge_remove(Tope* tope, Ridge* ridge)
{
  --tope->nridges;
  Allocator* alc = &tope->alc;
  int d = tope->dim;

  for (int i = 0; i < d - 1; ++i) {
    // Disassociate vertex:
    Vertex* vertex = ridge->verts[i];
    --vertex->nridges;
    assert(vertex->nridges >= 0);
    if (vertex->nridges == 0) {
      /* Remove the vertex as well: */
      vertex_remove(tope, vertex);
    }
  }

  // Clean up ridge struct:
  allocator_free(alc, ridge, sizeof(Ridge) + (d - 2) * sizeof(Vertex*));
}



/* Create and initialize a simplicial facet with d vertices: */
static Facet* facet_create(Tope* tope, Vertex** fverts, List* rli)
{
  ++tope->nfacets;
  Allocator* alc = &tope->alc;
  int d = tope->dim;

  assert(list_len(rli) == d);

  Facet* facet    = allocator_alloc(alc,     sizeof(Facet ));
  facet->centroid = allocator_alloc(alc, d * sizeof(double));
  facet->normal   = allocator_alloc(alc, d * sizeof(double));

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

  /* Add d ridges: */
  facet->ridges = rli;
  for (; rli; rli = rli->next) {
    /* Remember facet. */
    Ridge* ridge = rli->val;
    if (ridge->facets[0] == NULL) {
      assert(ridge->facets[0] == NULL);
      ridge->facets[0] = facet;
    }
    else {
      assert(ridge->facets[1] == NULL);
      ridge->facets[1] = facet;
    }
  }

  /* Vertex coordinates: */
  double* basis = alloca(d * d * sizeof(double));
  for (int i = 0; i < d; ++i) {
    memcpy(basis + i * d, fverts[i]->position, d * sizeof(double));
  }

  /* Extract basis, size, centroid: */
  analyzesimplex(d, d, basis, &facet->size, facet->centroid);

  /* Compute outward pointing normal: */
  memcpy(facet->normal, facet->centroid, d * sizeof(double));
  vec_sub(d, facet->normal, tope->center);
  for (int i = 0; i < d - 1; ++i) {
    double ip = vec_dot(d, facet->normal, basis + i * d);
    vec_adds(d, facet->normal, basis + i * d, -ip);
  }
  vec_normalize(d, facet->normal);

  return facet;
}



static void facet_free(Tope* tope, Facet* facet)
{
  int d = tope->dim;
  Allocator* alc = &tope->alc;
  while (facet->ridges) {
    Ridge* ridge = facet->ridges->val;
    if (ridge->facets[0] == facet) {
      assert(ridge->facets[1] != facet);
      ridge->facets[0] = ridge->facets[1];
      ridge->facets[1] = NULL;
    }
    else if (ridge->facets[1] == facet) {
      ridge->facets[1] = NULL;
    }

    if (ridge->facets[0] == NULL) {
      assert(ridge->facets[1] == NULL);
      ridge_remove(tope, ridge);
    }

    List* next = facet->ridges->next;
    allocator_free(alc, facet->ridges, sizeof(List));
    facet->ridges = next;
  }
  allocator_free(alc, facet->normal, d * sizeof(double));
  allocator_free(alc, facet->centroid, d * sizeof(double));
  allocator_free(alc, facet, sizeof(Facet));
}



static void facet_unlink(Tope* tope, Facet* facet)
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



/* Create initial simplex, add outside points, etc: */
static void initialsimplex(Tope* tope, int npoints, Point* points)
{
  int d = tope->dim;
  Allocator* alc = &tope->alc;

  /* Allocations: */
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

  /* Due to pivoting, the permutation vector will represent the point order
     which leads to the greatest size. */
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
      memswp(points + i + 1, points + pivot + 1, sizeof(Point));
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
  for (int i = 0; i <= d; ++i) {
    vec_add(d, tope->center, points[i].pos);
  }
  vec_scale(d, tope->center, 1.0 / (double)(tope->dim + 1));

  /* Simplex vertices: */
  Vertex** sverts = alloca((d + 1) * sizeof(Vertex*));
  for (int i = 0; i < d + 1; ++i) {
    sverts[i] = vertex_create(tope, points[i].pos, points[i].index);
  }

  int nridges = 0;
  Ridge** newridges = alloca(d * (d + 1) / 2 * sizeof(Ridge*));
  Vertex** fverts = alloca(d * sizeof(Vertex*));
  Vertex** rverts = alloca((d - 1) * sizeof(Vertex*));

  /* Create facets: */
  for (int i = 0; i < d + 1; ++i) {
    /* Facet vertices: */
    int nfv = 0;
    for (int j = 0; j < d + 1; ++j) {
      if (j == i) {
        continue;
      }
      fverts[nfv++] = sverts[j];
    }
    assert(nfv == d);

    /* Create ridges: */
    List* rli = NULL;
    List** prli = &rli;
    for (int j = 0; j < d; ++j) {
      /* Set d - 1 ridge verts: */
      int nrv = 0;
      for (int k = 0; k < d; ++k) {
        if (k == j) {
          continue;
        }
        rverts[nrv++] = fverts[k];
      }
      assert(nrv == d - 1);

      Ridge* ridge = NULL;
      for (int iridge = 0; iridge < nridges; ++iridge) {
        if (memcmp(rverts, newridges[iridge]->verts, (d - 1) * sizeof(Vertex*)) == 0) {
          ridge = newridges[iridge];
          break;
        }
      }
      if (ridge == NULL) {
        /* New ridge: */
        assert(nridges < d * (d + 1) / 2);
        ridge = ridge_create(tope, rverts);
        newridges[nridges++] = ridge;
      }
      prli = list_append(prli, ridge, alc);
    }

    /* Vertex array filled, compute and integrate the facet: */
    facet_create(tope, fverts, rli);
  }

  /* Create outside sets (assign each remaining point to a facet it 'sees'): */
  for (int i = tope->dim + 1; i < npoints; ++i) {
    /* Determine maximum distance: */
    double hmax = -HUGE_VAL;
    Facet* maxfacet = NULL;
    for (Facet* facet = tope->firstfacet; facet != NULL; facet = facet->next) {
      /* Distance to facet: */
      double h = height(d, points[i].pos, facet->centroid, facet->normal);
      if (h > hmax) {
        hmax = h;
        maxfacet = facet;
      }
    }

    /* If above facet, add it to outside set and move to next point. */
    if (hmax > EPS) {
      points[i].height = hmax;
      facet_addoutside(tope, maxfacet, points + i);
    }
  }

  /* Clean up: */
  free(span);
}



/* Add point to tope: */
static void addpoint(Tope* tope, Facet* facet, Point* apex)
{
  int d = tope->dim;
  Allocator* alc = &tope->alc;

  /* Add facet to visible list: */
  facet_unlink(tope, facet);
  facet->size = -1;  /* mark visible */
  Facet* visiblelist = facet;

  /* More lists: */
  Point* outsidepoints = NULL;
  List* horizonridges = NULL;
  List* newfacets = NULL;

  while (visiblelist != NULL) {
    /* Pop a value: */
    facet = visiblelist;
    visiblelist = facet->next;

    /* Loop over ridges to visit neighbours: */
    for (List* lst = facet->ridges; lst != NULL; lst = lst->next) {
      /* Ridge: */
      Ridge* ridge = lst->val;

      /* Retrieve neighbour: */
      Facet* neighbour = NULL;
      if (facet == ridge->facets[0]) {
        /* Current facet is 0, neighbour is 1: */
        ridge->facets[0] = ridge->facets[1]; /* forget reference to this facet */
        ridge->facets[1] = NULL;
      }
      else {
        /* Current facet is 1, neighbour is 0: */
        assert(facet == ridge->facets[1]);
        ridge->facets[1] = NULL; /* forget reference to this facet */
      }
      neighbour = ridge->facets[0];

      if (neighbour == NULL || neighbour->size < 0) {
        /* Removed or visited. */
        continue;
      }

      /* Neighbour is untested. */

      /* Height of apex above neighbour: */
      double h = height(d, apex->pos, neighbour->centroid, neighbour->normal);

      /* Decide if neighbour is visible: */
      if (h > EPS) {
        /* Remove neighbour from tope list: */
        facet_unlink(tope, neighbour);

        /* Prepend neighbour to visible list: */
        neighbour->next = visiblelist;
        visiblelist = neighbour;
        neighbour->size = -1;  /* mark visible */
      }
      else {
        /* Neighbour not visible. Ridge belongs to horizon: */
        list_prepend(&horizonridges, ridge, alc);
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

  /* New ridges hashmap: */
  HashMap newridges;
  hashmap_init(&newridges);

  Vertex** fverts = alloca(d * sizeof(Vertex*));
  Vertex** rverts = alloca((d - 1) * sizeof(Vertex*));

  /* Form new facets: */
  while (horizonridges) {

    /* Pop horizon ridge: */
    Ridge* horizonridge = list_pop(&horizonridges, alc);

    fverts[0] = vertex;
    for (int i = 1; i < d; ++i) {
      fverts[i] = horizonridge->verts[i - 1];
    }
    
    List* rli = NULL;
    List** prli = list_append(&rli, horizonridge, alc);
    for (int j = 1; j < d; ++j) {
      int nrv = 0;
      for (int k = 0; k < d; ++k) {
        if (k == j) {
          continue;
        }
        rverts[nrv++] = fverts[k];
      }
      Ridge** pridge = hashmap_get(&newridges, d - 1, rverts);
      if (*pridge == NULL) {
        *pridge = ridge_create(tope, rverts);
      }
      prli = list_append(prli, *pridge, alc);
    }

    Facet* facet = facet_create(tope, fverts, rli);
    list_prepend(&newfacets, facet, alc);
  }  // end loop over horizonridges

  hashmap_destroy(&newridges);

  /* Assign outside verts: */
  while (outsidepoints) {
    /* Pop: */
    Point* point = outsidepoints;
    outsidepoints = outsidepoints->next;
    point->next = NULL;

    /* Max height: */
    double hmax = -HUGE_VAL;
    Facet* maxfacet = NULL;
    for (List* fli = newfacets; fli; fli = fli->next) {
      Facet* facet = fli->val;

      /* Height of outside point above facet: */
      double h = height(tope->dim, point->pos, facet->centroid, facet->normal);
      if (h > hmax) {
        hmax = h;
        maxfacet = facet;
      }
    }

    if (hmax > EPS) {
      /* Outside. */
      point->height = hmax;
      facet_addoutside(tope, maxfacet, point);
    }
  }
  list_free(&newfacets, alc);
}



static void build(Tope* tope, int npoints, Point* points)
{
  assert(npoints >= tope->dim + 1);
  tope->maxindex = npoints - 1;
  initialsimplex(tope, npoints, points);

  /* Now, all vertices are either on the initial simplex OR
     added to the outside set of one of its facets OR
     inside the initial simplex. */
  while (tope->firstfacet->outsidehead) {
    /* Pop head: */
    Point* point = tope->firstfacet->outsidehead;
    tope->firstfacet->outsidehead = point->next;

    /* Add vertex to tope: */
    addpoint(tope, tope->firstfacet, point);
  }
}



static Tope* tope_new()
{
  Tope* tope = calloc(1, sizeof(Tope));
  allocator_init(&tope->alc);
  return tope;
}



void tope_delete(Tope* tope)
{
#ifdef USE_MALLOC
  while (tope->firstfacet) {
    Facet* next = tope->firstfacet->next;
    facet_free(tope, tope->firstfacet);
    tope->firstfacet = next;
  }
#endif
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
  tope_merge(rectope);

  /* Convert reciprocal tope to tope: */
  int nfacets = tope_getnumfacets(rectope);
  points = malloc(nfacets * d * sizeof(double));
  u32 ifacet = 0;
  for (Facet* facet = rectope->firstfacet;
       facet != NULL;
       facet = facet->next, ++ifacet) {
    /* Unscaled normal and distance: */
    double dist = vec_dot(d, facet->normal, facet->centroid);
    for (int i = 0; i < d; ++i) {
      points[ifacet * d + i] = facet->normal[i] / rectope->scales[i];
      dist += points[ifacet * d + i] * rectope->shift[i];
    }

    /* Straight space point: */
    vec_scale(d, points + ifacet * d, 1.0 / dist);
    vec_add(d, points + ifacet * d, xc);
  }

  Tope* tope = tope_frompoints(nfacets, d, points);
  free(points);
  tope_merge(tope);

  /* Clean up: */
  tope_delete(rectope);

  return tope;
}



Tope* tope_frompoints(int npoints, int dim, double const* orgpoints)
{
  Tope* tope = tope_new();

  /* Initialize member variables: */
  tope->dim = dim;
  tope->isdelaunay = false;
  tope->shift  = malloc(tope->dim * sizeof(double));
  tope->scales = malloc(tope->dim * sizeof(double));
  tope->center = malloc(tope->dim * sizeof(double));

  /* Bounding box: */
  int* minindices = alloca(tope->dim * sizeof(int));
  int* maxindices = alloca(tope->dim * sizeof(int));
  double* minima  = alloca(tope->dim * sizeof(double));
  double* maxima  = alloca(tope->dim * sizeof(double));
  boundingbox(
    npoints,
    tope->dim,
    orgpoints,
    minindices,
    maxindices,
    minima,
    maxima
  );

  /* Scale, offset: */
  memcpy(tope->scales, maxima, tope->dim * sizeof(double));
  vec_sub(tope->dim, tope->scales, minima);
  memcpy(tope->shift, minima, tope->dim * sizeof(double));

  /* Points array: */
  Point* points = malloc(npoints * sizeof(Point));
  double* positions = malloc(npoints * dim * sizeof(double));
  for (int ipoint = 0; ipoint < npoints; ++ipoint) {
    points[ipoint].next = NULL;
    points[ipoint].index = ipoint;
    points[ipoint].height = 0.0;
    points[ipoint].pos = positions + ipoint * dim;
    for (int idim = 0; idim < tope->dim; ++idim) {
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

  /* Remove upper Delaunay facets: */
  Facet* facet = tope->firstfacet;
  while (facet != NULL) {
    Facet* next = facet->next;
    if (facet->normal[dim] > -EPS) {
      // Upper delaunay
      facet_unlink(tope, facet);
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
    .index = ++tope->maxindex,
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
  Facet* facet = tope->firstfacet;
  for ( ; facet != NULL; facet = facet->next) {
    double h = height(tope->dim, pos, facet->centroid, facet->normal);
    if (h > EPS) {
      break;
    }
  }

  if (facet != NULL) {
    /* Add vertex to tope: */
    addpoint(tope, facet, &apex);
  }
}



void tope_print(Tope* tope)
{
  int d = tope->dim;

  double* sv       = alloca(d * sizeof(double));
  double* normal   = alloca(d * sizeof(double));
  double* centroid = alloca(d * sizeof(double));

  printf("%d facets\n", tope->nfacets);
  printf("%d ridges\n", tope->nridges);
  printf("%d vertices\n", tope->nverts);

  vec_reset(d, sv);

  int i = 0;
  for (Facet* facet = tope->firstfacet; facet != NULL; ++i, facet = facet->next) {
    tope_facet_getnormal(tope, facet, normal);
    double size = tope_facet_getsize(tope, facet);
    vec_adds(d, sv, normal, size);
    tope_facet_getcentroid(tope, facet, centroid);

    printf("facet %d\n", i + 1);
    printf("  size: %g\n", size);
    printf("  centroid: ");
    vec_print(d, centroid);
    printf("  normal: ");
    vec_print(d, normal);
    printf("  nridges = %d\n", list_len(facet->ridges));

    //int j = 0;
    //for (List* rli = facet->ridges; rli; rli = rli->next, ++j) {
    //  Ridge* ridge = rli->val;
    //  printf("  ridge %d\n", j + 1);
    //  printf("    size: %g\n", ridge->size);
    //  printf("    centroid: "); vec_print(d, ridge->centroid);
    //  printf("    normal1: "); vec_print(d, ridge->normal1);
    //  printf("    normal2: "); vec_print(d, ridge->normal2);
    //  printf("    nverts = %d\n", d - 1);
    //}
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
  assert(tope->isdelaunay);

  int d = tope->dim - 1;

  /* Some space for various tasks: */
  double* xiprime  = alloca(d * sizeof(double));
  double* normal   = alloca(d * sizeof(double));
  double* centroid = alloca(d * sizeof(double));
  double* basis    = alloca(d * d * sizeof(double));
  Vertex** fverts  = alloca((d + 1) * sizeof(Vertex));

  /* Transform xi to local coordinates: */
  for (int i = 0; i < d; ++i) {
    xiprime[i] = (xi[i] - tope->shift[i]) / tope->scales[i];
  }

  /* Start with first facet: */
  assert(tope->firstfacet != NULL);
  Facet* currentfacet = tope->firstfacet;

  double totalweight = 0.0;
  while (1) {
    /* In the current facet, find the ridge where xi is highest above: */
    totalweight = 0.0;
    double hmin = HUGE_VAL;
    Ridge* minridge = NULL;
    List* rli = currentfacet->ridges;
    Ridge* ridge = rli->val;
    for (int ivert = 0; ivert < d; ++ivert) {
      fverts[ivert + 1] = ridge->verts[ivert];
    }
    fverts[0] = ((Ridge*)rli->next->val)->verts[0];  // first vertex of second ridge

    for (int iridge = 0; iridge < tope->dim; ++iridge, rli = rli->next) {
      /* Retrieve ridge: */
      ridge = rli->val;

      for (int i = 0; i < d; ++i) {
        memcpy(basis + i * d, ridge->verts[i]->position, d * sizeof(double));
      }

      double size;
      analyzesimplex(d, d, basis, &size, centroid);

      /* Construct normal: */
      memcpy(normal, fverts[iridge]->position, d * sizeof(double));
      vec_sub(d, normal, centroid);
      for (int i = 0; i < d - 1; ++i) {
        double ip = vec_dot(d, normal, basis + i * d);
        vec_adds(d, normal, basis + i * d, -ip);
      }
      vec_normalize(d, normal);

      /* Height of xi above ridge: */
      double h = height(d, xiprime, centroid, normal);

      /* Weight for this ridge: */
      weights[iridge] = h * size;
      totalweight += weights[iridge];

      /* Index for this ridge: */
      indices[iridge] = fverts[iridge]->index;

      /* Keep track of minimum height (if not boundary ridge): */
      if (h < hmin && ridge->facets[1] != NULL) {
        assert(ridge->facets[0] != NULL);
        hmin = h;
        minridge = ridge;
      }
    }
    assert(rli == NULL);

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



double tope_facet_getsize(Tope* tope, Facet* facet)
{
  int d = tope->dim;
  double* s = tope->scales;
  double* n = facet->normal;
  double size = facet->size;

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

  size *= prod * sqrt(scale);

  return size;
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
  while (i--) {
    position[i] = vertex->position[i] * scales[i] + shift[i];
  }
}



u64 tope_bytes_used(Tope* tope)
{
  return tope->alc.used;
}



void tope_merge(Tope* tope)
{
  int d = tope->dim;
  Allocator* alc = &tope->alc;
  double* diff = alloca(d * sizeof(double));

  // Loop over facets:
  for (Facet* facet = tope->firstfacet; facet; facet = facet->next) {
    // Loop over ridges:
    List** rli_end = &facet->ridges;
    for (; *rli_end; rli_end = &(*rli_end)->next);  // scream to end of list
    for (List** prli = &facet->ridges; *prli;) {
      Ridge* divridge = (*prli)->val;  // dividing ridge
      Facet* neighbour = divridge->facets[0] == facet ? divridge->facets[1] : divridge->facets[0];
      if (facet == neighbour || neighbour == NULL) {
        // Remove dangling ridge:
        List* next = (*prli)->next;
        allocator_free(alc, *prli, sizeof(List));
        ridge_remove(tope, divridge);
        *prli = next;
        continue;
      }

      assert(neighbour && neighbour->size > 0);
      memcpy(diff, facet->normal, d * sizeof(double));
      vec_sub(d, diff, neighbour->normal);
      if (vec_nrmsq(d, diff) > EPS * EPS) {
        prli = &(*prli)->next;
      }
      else {
        // Parallel. Merge with facet:
        vec_scale(d, facet->centroid, facet->size);
        vec_adds(d, facet->centroid, neighbour->centroid, neighbour->size);
        facet->size += neighbour->size;
        vec_scale(d, facet->centroid, 1.0 / facet->size);
        while (neighbour->ridges) {
          // Pop value:
          List* nrli = neighbour->ridges;
          neighbour->ridges = neighbour->ridges->next;
          Ridge* nridge = nrli->val;

          // If ridge is dividing ridge, remove reference:
          if (nridge == divridge) {
            allocator_free(alc, nrli, sizeof(List));
            continue;
          }
          if (nridge->facets[0] == facet) {
            nridge->facets[1] = NULL;
            allocator_free(alc, nrli, sizeof(List));
            continue;
          }
          if (nridge->facets[1] == facet) {
            nridge->facets[0] = nridge->facets[1];
            nridge->facets[1] = NULL;
            allocator_free(alc, nrli, sizeof(List));
            continue;
          }

          // Append to facet ridges:
          *rli_end = nrli;
          rli_end = &(*rli_end)->next;
          *rli_end = NULL;

          // Assign ridge to facet:
          if (nridge->facets[0] == neighbour) {
            nridge->facets[0] = facet;
          }
          else {
            assert(nridge->facets[1] == neighbour);
            nridge->facets[1] = facet;
          }
        }
        assert(*rli_end == NULL);

        // Remove neighbour facet:
        facet_unlink(tope, neighbour);
        assert(neighbour->ridges == NULL);
        facet_free(tope, neighbour);

        // Remove dividing ridge:
        List* next = (*prli)->next;
        allocator_free(alc, *prli, sizeof(List));
        ridge_remove(tope, divridge);
        *prli = next;
      }
    }
  }
}
