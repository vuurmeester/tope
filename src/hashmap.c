#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "hashmap.h"
#include "util.h"

/* Initial capacity: */
#define MIN_CAP 32



static unsigned hashvertset(int d, polytoop_Vertex** verts)
{
  d -= 2; /* d - 1 vertices, so d - 2 is the last index */
  unsigned hash = verts[d]->index;
  while (d--) {
    hash = 31 * hash ^ verts[d]->index;
  }
  return hash;
}



static int compvertsets(int d, polytoop_Vertex** vertset1,
                        polytoop_Vertex** vertset2)
{
  --d; /* d - 1 vertices per ridge */
  while (d--) {
    if (vertset1[d]->index != vertset2[d]->index) {
      return -1;
    }
  }
  return 0;
}



static void resize(HashMap* hashmap, int newcap)
{
  Ridge** newridges = calloc(newcap, sizeof(Ridge*) + sizeof(unsigned));
  unsigned* newhashes = (unsigned*)(newridges + newcap);
  unsigned newmask = newcap - 1;
  for (int i = 0; i < hashmap->cap; ++i) {
    if (!hashmap->ridges[i]) {
      continue;
    }
    unsigned newindex = hashmap->hashes[i] & newmask;
    while (newridges[newindex]) {
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
  hashmap->ridges = calloc(MIN_CAP, sizeof(Ridge*) + sizeof(unsigned));
  hashmap->hashes = (unsigned*)(hashmap->ridges + MIN_CAP);
}



void hashmap_destroy(HashMap* hashmap)
{
  free(hashmap->ridges);
  memset(hashmap, 0, sizeof(HashMap));
}



void hashmap_insert(HashMap* hashmap, int d, Ridge* ridge)
{
  unsigned hash = hashvertset(d, ridge->vertices);
  int index = hash & (hashmap->cap - 1);
  while (hashmap->ridges[index]) {
    /* Don't insert stuff that is already in here: */
    assert(compvertsets(d, hashmap->ridges[index]->vertices, ridge->vertices));
    index = (index + 1) & (hashmap->cap - 1); /* next in cluster */
  }

  hashmap->ridges[index] = ridge;
  hashmap->hashes[index] = hash;

  ++hashmap->len;

  if (3 * hashmap->len > 2 * hashmap->cap) {
    resize(hashmap, hashmap->cap << 1);
  }
}



void hashmap_clear(HashMap* hashmap)
{
  memset(hashmap->ridges, 0, hashmap->cap * sizeof(Ridge*));
  memset(hashmap->hashes, 0, hashmap->cap * sizeof(unsigned));
  hashmap->len = 0;
}



Ridge* hashmap_retrieve(HashMap hashmap, int d, polytoop_Vertex** verts)
{
  unsigned hash = hashvertset(d, verts);
  int index = hash & (hashmap.cap - 1);

  while (hashmap.ridges[index]) {
    if (hashmap.hashes[index] == hash &&
        compvertsets(d, verts, hashmap.ridges[index]->vertices) == 0) {
      return hashmap.ridges[index];
    }
    index = (index + 1) & (hashmap.cap - 1); /* next in cluster */
  }

  return NULL;
}
