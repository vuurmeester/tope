#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "hashmap.h"
#include "util.h"

/* Initial capacity: */
#define MIN_CAP 32



static u32 hashvertset(int d, u32* verts)
{
  d -= 2; /* d - 1 vertices, so d - 2 is the last index */
  u32 hash = verts[d];
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
  return 0; /* same */
}



static void initarrays(u32** ridges, u32** hashes, u32 cap)
{
  *ridges = malloc(cap * 2 * sizeof(u32));
  memset(*ridges, 0xff, cap * sizeof(u32));
  *hashes = *ridges + cap;
}



static void expand(HashMap* hashmap)
{
  u32 newcap = hashmap->cap * 2;
  u32 newmask = newcap - 1;
  u32* newridges;
  u32* newhashes;
  initarrays(&newridges, &newhashes, newcap);

  for (u32 i = 0; i < hashmap->cap; ++i) {
    if (hashmap->ridges[i] == UINT32_MAX) {
      continue;
    }
    u32 newindex = hashmap->hashes[i] & newmask;
    while (newridges[newindex] != UINT32_MAX) {
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
  memset(hashmap->ridges, 0xff, hashmap->cap * sizeof(u32));
  hashmap->len = 0;
}



u32* hashmap_get(HashMap* hashmap, int d, u32* verts, Allocator* alc)
{
  if (4 * hashmap->len > 3 * hashmap->cap) {
    /* More than 75% filled. */
    expand(hashmap);
  }

  u32 hash = hashvertset(d, verts);
  u32 index = hash & (hashmap->cap - 1);

  while (hashmap->ridges[index] != UINT32_MAX) {
    if (hashmap->hashes[index] == hash) {
      Ridge* ridge = allocator_mem(alc, hashmap->ridges[index]);
      if (compvertsets(d, verts, ridge->vertices) == 0) {
        return hashmap->ridges + index;
      }
    }
    index = (index + 1) & (hashmap->cap - 1); /* next in cluster */
  }

  hashmap->hashes[index] = hash;

  ++hashmap->len;

  return hashmap->ridges + index;
}
