#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "allocator.h"
#include "hashmap.h"
#include "util.h"

/* Initial capacity: */
#define MIN_CAP 32



static u32 hashvertset(int d, tope_Vertex** verts)
{
  u32 hash = 0x811c9dc5;  /* fnv-1a */
  d -= 1; /* d - 1 vertices */
  while (d--) {
    hash = hash ^ (u32)((u64)verts[d] >> 3);
    hash *= 0x01000193;
  }
  return hash;
}



static bool vertsets_equal(int d, Vertex** vertset1, Vertex** vertset2)
{
  --d; /* d - 1 vertices per ridge */
  while (d--) {
    if (vertset1[d] != vertset2[d]) {
      return false;
    }
  }
  return true; /* same */
}



static void initarrays(Ridge*** ridges, u32** hashes, u32 cap)
{
  *ridges = malloc(cap * (sizeof(Ridge*) + sizeof(u32)));
  memset(*ridges, 0x00, cap * sizeof(Ridge*));
  *hashes = (u32*)(*ridges + cap);
}



static void expand(HashMap* hashmap)
{
  u32 newcap = hashmap->cap * 2;
  u32 newmask = newcap - 1;
  Ridge** newridges;
  u32* newhashes;
  initarrays(&newridges, &newhashes, newcap);

  for (u32 i = 0; i < hashmap->cap; ++i) {
    if (hashmap->ridges[i] == NULL) {
      continue;
    }
    u32 newindex = hashmap->hashes[i] & newmask;
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
  initarrays(&hashmap->ridges, &hashmap->hashes, MIN_CAP);
}



void hashmap_destroy(HashMap* hashmap)
{
  free(hashmap->ridges);
  memset(hashmap, 0, sizeof(HashMap));
}



void hashmap_clear(HashMap* hashmap)
{
  memset(hashmap->ridges, 0x00, hashmap->cap * sizeof(Ridge*));
  hashmap->len = 0;
}



Ridge** hashmap_get(HashMap* hashmap, int d, Vertex** verts, Allocator* alc)
{
  if (4 * hashmap->len > 3 * hashmap->cap) {
    /* More than 75% filled. */
    expand(hashmap);
  }

  u32 hash = hashvertset(d, verts);
  u32 index = hash & (hashmap->cap - 1);

  while (hashmap->ridges[index] != NULL) {
    if (hashmap->hashes[index] == hash) {
      Ridge* ridge = hashmap->ridges[index];
      if (vertsets_equal(d, verts, ridge->vertices)) {
        return hashmap->ridges + index;
      }
    }
    index = (index + 1) & (hashmap->cap - 1); /* next in cluster */
  }

  /* Not found, create new entry and return reference: */
  hashmap->hashes[index] = hash;
  ++hashmap->len;
  return hashmap->ridges + index;
}
