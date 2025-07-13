#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "hashmap.h"

/* Initial capacity: */
#define MIN_CAP 8



static unsigned hashvertset(int d, polytoop_Vertex** verts)
{
  int i;
  int nverts;
  unsigned ret;

  /* Vertices: */
  nverts = d - 1;

  /* Sum vertex id's: */
  ret = 0;
  for (i = 0; i < nverts; ++i) {
    ret += verts[i]->index;
  }

  return ret;
}



/* Compare sets of d - 1 sorted vertices (for ordering ridges): */
static int compvertsets(int d, polytoop_Vertex** vertset1,
                        polytoop_Vertex** vertset2)
{
  int i;

  for (i = 0; i < d - 1; ++i) {
    if (vertset1[i]->index < vertset2[i]->index) {
      return -1;
    }
    if (vertset1[i]->index > vertset2[i]->index) {
      return 1;
    }
  }

  return 0;
}



static void resize(HashMap* hashmap, int newcap)
{
  int i;
  int index;

  Ridge** newridges = calloc(newcap, sizeof(Ridge*) + sizeof(unsigned));
  unsigned* newhashes = (unsigned*)(newridges + newcap);
  unsigned newmask = newcap - 1;
  for (i = 0; i < hashmap->cap; ++i) {
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



HashMap* hashmap_new()
{
  int i;
  HashMap* hashmap;

  hashmap = malloc(sizeof(HashMap));
  hashmap->cap = MIN_CAP;
  hashmap->len = 0;
  hashmap->ridges = calloc(MIN_CAP, sizeof(Ridge*) + sizeof(unsigned));
  hashmap->hashes = (unsigned*)(hashmap->ridges + MIN_CAP);

  return hashmap;
}



void hashmap_delete(HashMap* hashmap)
{
  free(hashmap->ridges);
  free(hashmap);
}



void hashmap_insert(HashMap* hashmap, int d, Ridge* ridge)
{
  unsigned hash;
  int index;

  hash = hashvertset(d, ridge->vertices);
  index = hash & (hashmap->cap - 1);
  while (hashmap->ridges[index]) {
    assert(hashmap->hashes[index] != hash ||
           compvertsets(d, hashmap->ridges[index]->vertices, ridge->vertices) !=
               0);
    index = (index + 1) & (hashmap->cap - 1);
  }

  /* New entry: */
  hashmap->ridges[index] = ridge;
  hashmap->hashes[index] = hash;

  /* Increase size: */
  ++hashmap->len;

  /* Double capacity if loadfactor > 2/3 */
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



Ridge* hashmap_retrieve(HashMap* hashmap, int d, polytoop_Vertex** verts)
{
  unsigned hash;
  int index;

  /* Generate hash: */
  hash = hashvertset(d, verts);

  /* Get proper chain: */
  index = hash & (hashmap->cap - 1);

  while (hashmap->ridges[index]) {
    if (hashmap->hashes[index] == hash &&
        compvertsets(d, verts, hashmap->ridges[index]->vertices) == 0) {
      return hashmap->ridges[index];
    }
    index = (index + 1) & (hashmap->cap - 1);
  }

  return NULL;
}
