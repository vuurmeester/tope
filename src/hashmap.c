#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "hashmap.h"
#include "util.h"

/* Initial capacity: */
#define MIN_CAP 32



static uint32_t hashvertset(int d, u32* verts)
{
  uint32_t hash;

  d -= 2; /* d - 1 vertices, so d - 2 is the last index */
  hash = verts[d];
  while (d--) {
    hash = 31 * hash ^ verts[d];
  }
  return hash;
}



static int compvertsets(int d, u32* vertset1, u32* vertset2)
{
  --d; /* d - 1 vertices per ridge */
  while (d--) {
    if (vertset1[d] != vertset2[d]) {
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
  u32* newridges;

  newcap = hashmap->cap * 2;
  newridges = calloc(2 * newcap, sizeof(u32));
  memset(newridges, 0xff, newcap * sizeof(u32));
  newhashes = (u32*)(newridges + newcap);
  newmask = newcap - 1;

  for (i = 0; i < hashmap->cap; ++i) {
    if (hashmap->ridges[i] == (u32)-1) {
      continue;
    }
    newindex = hashmap->hashes[i] & newmask;
    while (newridges[newindex] != (u32)-1) {
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
  hashmap->ridges = malloc(2 * MIN_CAP * sizeof(u32));
  memset(hashmap->ridges, 0xff, MIN_CAP * sizeof(u32));
  hashmap->hashes = (u32*)(hashmap->ridges + MIN_CAP);
}



void hashmap_destroy(HashMap* hashmap)
{
  free(hashmap->ridges);
  memset(hashmap, 0, sizeof(HashMap));
}



void hashmap_insert(HashMap* hashmap, int d, u32* verts, u32 ridge)
{
  uint32_t hash;
  uint32_t index;

  if (5 * hashmap->len > 4 * hashmap->cap) {
    /* More than 4/5 filled. */
    expand(hashmap);
  }

  hash = hashvertset(d, verts);
  index = hash & (hashmap->cap - 1);

  while (hashmap->ridges[index] != (u32)-1) {
    /* Don't insert stuff that is already in here: */
    assert(hashmap->ridges[index] != ridge);
    index = (index + 1) & (hashmap->cap - 1); /* next in cluster */
  }

  hashmap->ridges[index] = ridge;
  hashmap->hashes[index] = hash;

  ++hashmap->len;
}



void hashmap_clear(HashMap* hashmap)
{
  memset(hashmap->ridges, 0xff, hashmap->cap * sizeof(u32));
  hashmap->len = 0;
}



u32 hashmap_retrieve(HashMap hashmap, int d, u32* verts, Allocator* alc)
{
  uint32_t hash;
  uint32_t index;

  hash = hashvertset(d, verts);
  index = hash & (hashmap.cap - 1);

  while (hashmap.ridges[index] != (u32)-1) {
    if (hashmap.hashes[index] == hash) {
      Ridge* ridge = allocator_mem(alc, hashmap.ridges[index]);
      if (compvertsets(d, verts, ridge->vertices) == 0) {
        return hashmap.ridges[index];
      }
    }
    index = (index + 1) & (hashmap.cap - 1); /* next in cluster */
  }

  return (u32)-1;
}
