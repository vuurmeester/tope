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

  hashmap->cap = newcap;
  hashmap->mask = newcap - 1;
  hashmap->entries = realloc(hashmap->entries, newcap * sizeof(HashEntry*));
  for (i = 0; i < newcap; ++i) {
    hashmap->entries[i] = NULL;
  }

  /* Find a spot for all the entries in the new array: */
  for (entry = hashmap->first; entry != NULL; entry = entry->nextinlist) {
    index = entry->hash & hashmap->mask; /* new index */
    entry->nextinchain = NULL;
    pentry = &hashmap->entries[index];
    while (*pentry) {
      pentry = &(*pentry)->nextinchain;
    }
    *pentry = entry;
  }
}



HashMap* hashmap_new(unsigned (*hashfunc)(void* key, void* data),
                     int (*compar)(void* key1, void* key2, void* data))
{
  int i;
  HashMap* hashmap;

  hashmap = malloc(sizeof(HashMap));
  hashmap->cap = MIN_CAP;
  hashmap->mask = MIN_CAP - 1;
  hashmap->len = 0;
  hashmap->data = NULL;
  hashmap->hashfunc = hashfunc;
  hashmap->compar = compar;
  hashmap->entries = malloc(MIN_CAP * sizeof(HashEntry*));
  for (i = 0; i < MIN_CAP; ++i) {
    hashmap->entries[i] = NULL;
  }
  hashmap->allocator = allocator_new();
  hashmap->first = NULL;
  hashmap->last = NULL;

  return hashmap;
}



void hashmap_delete(HashMap* hashmap)
{
  /* Empty the hashmap: */
  hashmap_clear(hashmap);

  /* Release allocator: */
  allocator_delete(hashmap->allocator);

  /* Release chain pointers: */
  free(hashmap->entries);

  /* Release object itself: */
  free(hashmap);
}



void hashmap_setdata(HashMap* hashmap, void* data) { hashmap->data = data; }



HashEntry* hashmap_insert(HashMap* hashmap, void* key, void* value)
{
  unsigned hash;
  int index;
  HashEntry** pentry;
  HashEntry* entry;

  hash = hashmap->hashfunc(key, hashmap->data);
  index = hash & hashmap->mask;
  pentry = &hashmap->entries[index];
  while (*pentry) {
    if ((*pentry)->hash == hash &&
        hashmap->compar(key, (*pentry)->key, hashmap->data) == 0) {
      /* Found existing key/value pair. Replace them: */
      (*pentry)->key = key;
      (*pentry)->value = value;

      /* Done. */
      return *pentry;
    }

    pentry = &(*pentry)->nextinchain;
  }

  /* New entry: */
  *pentry = allocator_alloc(hashmap->allocator, sizeof(HashEntry));
  (*pentry)->hash = hash;
  (*pentry)->key = key;
  (*pentry)->value = value;
  (*pentry)->nextinchain = NULL;

  /* Append to linked list: */
  (*pentry)->nextinlist = NULL;
  if (!hashmap->first) {
    assert(hashmap->last == NULL && hashmap->len == 0);
    hashmap->first = *pentry;
    (*pentry)->previnlist = NULL;
  } else {
    assert(hashmap->last->nextinlist == NULL);
    hashmap->last->nextinlist = *pentry;
    (*pentry)->previnlist = hashmap->last;
  }
  hashmap->last = *pentry;

  /* Increase size: */
  ++hashmap->len;

  entry = *pentry;

  /* Double capacity if loadfactor > 2/3 */
  if (3 * hashmap->len > 2 * hashmap->cap) {
    resize(hashmap, hashmap->cap << 1);
  }

  return entry;
}



void hashmap_clear(HashMap* hashmap)
{
  int i;
  HashEntry* entry;

  /* Release all entries: */
#ifndef NDEBUG
  entry = hashmap->first;
  while (entry) {
    HashEntry* next;
    next = entry->nextinlist;
    allocator_free(hashmap->allocator, entry, sizeof(HashEntry));
    entry = next;
  }
#endif
  allocator_clear(hashmap->allocator);

  /* No more elements: */
  hashmap->cap = MIN_CAP;
  hashmap->mask = MIN_CAP - 1;
  hashmap->len = 0;
  hashmap->entries = realloc(hashmap->entries, MIN_CAP * sizeof(HashEntry*));
  for (i = 0; i < MIN_CAP; ++i) {
    hashmap->entries[i] = NULL;
  }
  hashmap->first = NULL;
  hashmap->last = NULL;
}



void* hashmap_retrieve(HashMap* hashmap, void* key)
{
  unsigned hash;
  int index;
  HashEntry* entry;

  if (hashmap->entries != NULL) {
    /* Generate hash: */
    hash = hashmap->hashfunc(key, hashmap->data);

    /* Get proper chain: */
    index = hash & hashmap->mask;
    entry = hashmap->entries[index];

    while (entry != NULL) {
      if (entry->hash == hash &&
          hashmap->compar(key, entry->key, hashmap->data) == 0) {
        return entry->value;
      }
      entry = entry->nextinchain;
    }
  }

  return NULL;
}
