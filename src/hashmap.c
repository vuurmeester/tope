#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "allocator.h"
#include "hashmap.h"
#include "util.h"

/* Initial capacity: */
#define MIN_CAP 32



static u32 hashvertset(int nverts, Vertex** verts)
{
  u32 hash = 0x811c9dc5;  /* fnv-1a */
  while (nverts--) {
    Vertex* vert = verts[nverts];
    hash = hash ^ vert->index;
    hash *= 0x01000193;
  }
  return hash;
}



static void initarrays(Ridge*** ridges, u32** hashes, u32 cap)
{
  /* Allocate ridge and hash arrays together: */
  int numbytes = cap * (sizeof(Ridge*) + sizeof(u32));
  *ridges = malloc(numbytes);
  memset(*ridges, 0x00, numbytes);
  *hashes = (u32*)(*ridges + cap);  /* hash array is past ridge array */
}



static void expand(HashMap* hashmap)
{
  u32 newcap = hashmap->cap * 2;  /* double capacity */
  u32 newmask = newcap - 1;

  /* Allocate new arrays: */
  Ridge** newridges;
  u32* newhashes;
  initarrays(&newridges, &newhashes, newcap);

  /* Put existing entries in new arrays: */
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

  /* Release old arrays and assign new ones: */
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



Ridge** hashmap_get(HashMap* hashmap, int nverts, Vertex** vertset)
{
  if (4 * hashmap->len >= 3 * hashmap->cap) {
    /* More than 75% filled. */
    expand(hashmap);
  }

  u32 hash = hashvertset(nverts, vertset);
  u32 index = hash & (hashmap->cap - 1);

  while (hashmap->ridges[index] != NULL) {
    if (hashmap->hashes[index] == hash &&
        memcmp(vertset, hashmap->ridges[index]->verts, nverts * sizeof(Vertex*)) == 0) {
      /* same hash and same vertset */
      return hashmap->ridges + index;
    }
    index = (index + 1) & (hashmap->cap - 1);  /* next in cluster */
  }

  /* Not found, create new entry and return reference (so that it can be filled): */
  hashmap->hashes[index] = hash;
  ++hashmap->len;
  return hashmap->ridges + index;
}
