#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "hashmap.h"

/* Minimum capacity: */
#define MIN_CAP 8



static void resize(HashMap* hashmap, int newcap)
{
  int i;
  int index;
  HashEntry* entry;
  HashEntry** pentry;

  HashEntry* newentries = calloc(newcap, sizeof(HashEntry));
  for (i = 0; i < hashmap->cap; ++i) {
    entry = &hashmap->entries[i];
    if (!entry->ridge) {
      continue;
    }
    unsigned newindex = entry->hash & (newcap - 1);
    newentries[newindex] = hashmap->entries[i];
  }

  free(hashmap->entries);
  hashmap->entries = newentries;
  hashmap->cap = newcap;
}



HashMap* hashmap_new(unsigned (*hashfunc)(void* key, void* data),
                     int (*compar)(void* key1, void* key2, void* data))
{
  int i;
  HashMap* hashmap;

  hashmap = malloc(sizeof(HashMap));
  hashmap->cap = MIN_CAP;
  hashmap->len = 0;
  hashmap->data = NULL;
  hashmap->hashfunc = hashfunc;
  hashmap->compar = compar;
  hashmap->entries = calloc(MIN_CAP, sizeof(HashEntry));

  return hashmap;
}



void hashmap_delete(HashMap* hashmap)
{
  free(hashmap->entries);
  free(hashmap);
}



void hashmap_setdata(HashMap* hashmap, void* data) { hashmap->data = data; }



void hashmap_insert(HashMap* hashmap, Ridge* ridge)
{
  unsigned hash;
  int index;
  HashEntry* entry;

  hash = hashmap->hashfunc(ridge->vertices, hashmap->data);
  index = hash & (hashmap->cap - 1);
  entry = &hashmap->entries[index];
  while (entry->ridge) {
    assert(entry->hash != hash ||
           hashmap->compar(entry->ridge->vertices, ridge->vertices,
                           hashmap->data) != 0);
    index = (index + 1) & (hashmap->cap - 1);
    entry = &hashmap->entries[index];
  }

  /* New entry: */
  entry->hash = hash;
  entry->ridge = ridge;

  /* Increase size: */
  ++hashmap->len;

  /* Double capacity if loadfactor > 2/3 */
  if (3 * hashmap->len > 2 * hashmap->cap) {
    resize(hashmap, hashmap->cap << 1);
  }
}



void hashmap_clear(HashMap* hashmap)
{
  memset(hashmap->entries, 0, hashmap->cap * sizeof(HashEntry));
  hashmap->len = 0;
}



Ridge* hashmap_retrieve(HashMap* hashmap, polytoop_Vertex** verts)
{
  unsigned hash;
  int index;
  HashEntry* entry;

  if (hashmap->entries != NULL) {
    /* Generate hash: */
    hash = hashmap->hashfunc(verts, hashmap->data);

    /* Get proper chain: */
    index = hash & (hashmap->cap - 1);
    entry = &hashmap->entries[index];

    while (entry->ridge) {
      if (entry->hash == hash &&
          hashmap->compar(verts, entry->ridge->vertices, hashmap->data) == 0) {
        return entry->ridge;
      }
      index = (index + 1) & (hashmap->cap - 1);
      entry = &hashmap->entries[index];
    }
  }

  return NULL;
}
