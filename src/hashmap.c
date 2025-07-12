#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"
#include "hashmap.h"

/* Minimum capacity: */
#define MIN_CAP 8

struct _HashEntry {
  unsigned hash;
  void* key;
  void* value;
  HashEntry* nextinchain;
  HashEntry* nextinlist;
  HashEntry* previnlist;
};

struct _HashMap {
  int capacity;                       /* current capacity */
  unsigned mask;                      /* for efficient modulo */
  int size;                           /* number of elements */
  void* data;                         /* userdata */
  unsigned (*hashfunc)(void*, void*); /* hash function */
  int (*compar)(void*, void*, void*); /* key comparison function */
  void (*freekey)(void*);             /* key deletion function */
  void (*freeval)(void*);             /* value deletion function */
  void (*freedata)(void*);            /* user data deletion function */
  HashEntry** entries;                /* chain pointers */
  Allocator* allocator;
  HashEntry* first; /* first entry in linked list */
  HashEntry* last;  /* last entry in linked list */
};



static void resize(HashMap* hashmap, int newcap)
{
  int i;
  int index;
  HashEntry* entry;
  HashEntry** pentry;

  hashmap->capacity = newcap;
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
  hashmap->capacity = MIN_CAP;
  hashmap->mask = MIN_CAP - 1;
  hashmap->size = 0;
  hashmap->data = NULL;
  hashmap->hashfunc = hashfunc;
  hashmap->compar = compar;
  hashmap->freekey = NULL;
  hashmap->freeval = NULL;
  hashmap->freedata = NULL;
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



void hashmap_setfreekey(HashMap* hashmap, void (*freekey)(void*))
{
  hashmap->freekey = freekey;
}



void hashmap_setfreeval(HashMap* hashmap, void (*freeval)(void*))
{
  hashmap->freeval = freeval;
}



void hashmap_setfreedata(HashMap* hashmap, void (*freedata)(void* data))
{
  hashmap->freedata = freedata;
}



int hashmap_getsize(HashMap* hashmap) { return hashmap->size; }



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
      if (hashmap->freekey && key != (*pentry)->key) {
        hashmap->freekey((*pentry)->key);
      }
      if (hashmap->freeval && value != (*pentry)->value) {
        hashmap->freeval((*pentry)->value);
      }

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
    assert(hashmap->last == NULL && hashmap->size == 0);
    hashmap->first = *pentry;
    (*pentry)->previnlist = NULL;
  } else {
    assert(hashmap->last->nextinlist == NULL);
    hashmap->last->nextinlist = *pentry;
    (*pentry)->previnlist = hashmap->last;
  }
  hashmap->last = *pentry;

  /* Increase size: */
  ++hashmap->size;

  entry = *pentry;

  /* Double capacity if loadfactor > 80% */
  if (10 * hashmap->size > 8 * hashmap->capacity) {
    resize(hashmap, hashmap->capacity << 1);
  }

  return entry;
}



static void remove(HashMap* hashmap, HashEntry** pentry)
{
  HashEntry* nextinchain;
  HashEntry* next;
  HashEntry* prev;

  /* Remove entry: */
  nextinchain = (*pentry)->nextinchain;
  next = (*pentry)->nextinlist;
  prev = (*pentry)->previnlist;

  /* Deallocate: */
  if (hashmap->freekey) {
    hashmap->freekey((*pentry)->key);
  }
  if (hashmap->freeval) {
    hashmap->freeval((*pentry)->value);
  }
  allocator_free(hashmap->allocator, *pentry, sizeof(HashEntry));

  /* Update chain: */
  *pentry = nextinchain;

  /* Update linked list: */
  if (next) {
    next->previnlist = prev;
  } else {
    hashmap->last = prev;
  }
  if (prev) {
    prev->nextinlist = next;
  } else {
    hashmap->first = next;
  }

  /* Decrement size: */
  --hashmap->size;

  /* Halve capacity if loadfactor < 30% */
  if (10 * hashmap->size < 3 * hashmap->capacity &&
      hashmap->capacity > MIN_CAP) {
    resize(hashmap, hashmap->capacity >> 1);
  }
}



void hashmap_remove_entry(HashMap* hashmap, HashEntry* entry)
{
  unsigned index;
  HashEntry** pentry;

  /* Get chain: */
  index = entry->hash & hashmap->mask;
  pentry = &hashmap->entries[index];

  /* Find key in chain: */
  while (*pentry != entry) {
    pentry = &(*pentry)->nextinchain;
  }

  remove(hashmap, pentry);
}



void hashmap_remove(HashMap* hashmap, void* key)
{
  /* Generate hash: */
  unsigned hash = hashmap->hashfunc(key, hashmap->data);

  /* Get chain: */
  unsigned index = hash & hashmap->mask;
  HashEntry** pentry = &hashmap->entries[index];

  /* Find key in chain: */
  while (*pentry) {
    if ((*pentry)->hash == hash &&
        hashmap->compar(key, (*pentry)->key, hashmap->data) == 0) {
      /* Found it. Remove entry: */
      remove(hashmap, pentry);

      /* Done. */
      break;
    }

    /* Try next in chain: */
    pentry = &(*pentry)->nextinchain;
  }
}



void hashmap_clear(HashMap* hashmap)
{
  int i;
  HashEntry* entry;

  /* Call free functions: */
  if (hashmap->freekey) {
    for (entry = hashmap->first; entry != NULL; entry = entry->nextinlist) {
      hashmap->freekey(entry->key);
    }
  }
  if (hashmap->freeval) {
    for (entry = hashmap->first; entry != NULL; entry = entry->nextinlist) {
      hashmap->freeval(entry->value);
    }
  }

  /* Free user data: */
  if (hashmap->freedata) {
    hashmap->freedata(hashmap->data);
  }

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
  hashmap->capacity = MIN_CAP;
  hashmap->mask = MIN_CAP - 1;
  hashmap->size = 0;
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



void hashmap_getarrays(HashMap* hashmap, void** keys, void** values)
{
  int i;
  HashEntry* entry;

  if (keys) {
    i = 0;
    for (entry = hashmap->first; entry != NULL; entry = entry->nextinlist) {
      keys[i++] = entry->key;
    }
    assert(i == hashmap->size);
  }
  if (values) {
    i = 0;
    for (entry = hashmap->first; entry != NULL; entry = entry->nextinlist) {
      values[i++] = entry->value;
    }
    assert(i == hashmap->size);
  }
}



HashEntry* hashmap_first(HashMap* hashmap) { return hashmap->first; }



HashEntry* hashmap_last(HashMap* hashmap) { return hashmap->last; }



HashEntry* hashentry_next(HashEntry* entry) { return entry->nextinlist; }



HashEntry* hashentry_prev(HashEntry* entry) { return entry->previnlist; }



void* hashentry_getkey(HashEntry* entry) { return entry->key; }



void* hashentry_getvalue(HashEntry* entry) { return entry->value; }
