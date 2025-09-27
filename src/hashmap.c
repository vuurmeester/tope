#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "hashmap.h"
#include "util.h"

/* Initial capacity: */
#define MIN_CAP 32



static uint32_t hashvertset(int d, polytoop_Vertex** verts)
{
  uint32_t hash;

  d -= 2; /* d - 1 vertices, so d - 2 is the last index */
  hash = verts[d]->index;
  while (d--) {
    hash = 31 * hash ^ verts[d]->index;
  }
  return hash;
}



static int compvertsets(
    int d,
    polytoop_Vertex** vertset1,
    polytoop_Vertex** vertset2
)
{
  --d; /* d - 1 vertices per ridge */
  while (d--) {
    if (vertset1[d]->index != vertset2[d]->index) {
      return -1;
    }
  }
  return 0;
}



static void expand(HashMap* hashmap)
{
  uint32_t newcap;
  uint32_t i;
  uint32_t newmask;
  uint32_t newindex;
  uint32_t* newhashes;
  Ridge** newridges;

  newcap = hashmap->cap * 2;
  newridges = calloc(newcap, sizeof(Ridge*) + sizeof(uint32_t));
  newhashes = (uint32_t*)(newridges + newcap);
  newmask = newcap - 1;

  for (i = 0; i < hashmap->cap; ++i) {
    if (hashmap->ridges[i] == NULL) {
      continue;
    }
    newindex = hashmap->hashes[i] & newmask;
    while (newridges[newindex] != NULL) {
      newindex = (newindex + 1) & newmask;
    }
    newridges[newindex] = hashmap->ridges[i];
    newhashes[newindex] = hashmap->hashes[i];
  }

  free(hashmap->ridges);
  hashmap->ridges = newridges;
  hashmap->hashes = newhashes;
  hashmap->cap = newcap;
}



void hashmap_init(HashMap* hashmap)
{
  hashmap->cap = MIN_CAP;
  hashmap->len = 0;
  hashmap->ridges = calloc(MIN_CAP, sizeof(Ridge*) + sizeof(uint32_t));
  hashmap->hashes = (uint32_t*)(hashmap->ridges + MIN_CAP);
}



void hashmap_destroy(HashMap* hashmap)
{
  free(hashmap->ridges);
  memset(hashmap, 0, sizeof(HashMap));
}



void hashmap_insert(HashMap* hashmap, int d, Ridge* ridge)
{
  uint32_t hash;
  uint32_t index;

  if (5 * hashmap->len > 4 * hashmap->cap) {
    /* More than 4/5 filled. */
    expand(hashmap);
  }

  hash = hashvertset(d, ridge->vertices);
  index = hash & (hashmap->cap - 1);

  while (hashmap->ridges[index] != NULL) {
    /* Don't insert stuff that is already in here: */
    assert(
        compvertsets(d, hashmap->ridges[index]->vertices, ridge->vertices) != 0
    );
    index = (index + 1) & (hashmap->cap - 1); /* next in cluster */
  }

  hashmap->ridges[index] = ridge;
  hashmap->hashes[index] = hash;

  ++hashmap->len;
}



void hashmap_clear(HashMap* hashmap)
{
  memset(hashmap->ridges, 0, hashmap->cap * sizeof(Ridge*));
  memset(hashmap->hashes, 0, hashmap->cap * sizeof(uint32_t));
  hashmap->len = 0;
}



Ridge* hashmap_retrieve(HashMap hashmap, int d, polytoop_Vertex** verts)
{
  uint32_t hash;
  uint32_t index;

  hash = hashvertset(d, verts);
  index = hash & (hashmap.cap - 1);

  while (hashmap.ridges[index] != NULL) {
    if (hashmap.hashes[index] == hash &&
        compvertsets(d, verts, hashmap.ridges[index]->vertices) == 0) {
      return hashmap.ridges[index];
    }
    index = (index + 1) & (hashmap.cap - 1); /* next in cluster */
  }

  return NULL;
}
