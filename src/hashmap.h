#ifndef __hashmap_h__
#define __hashmap_h__

#include "allocator.h"

typedef struct _HashMap HashMap;
typedef struct _HashEntry HashEntry;

/** Create hashmap. */
HashMap* hashmap_new(
  unsigned(*hashfunc)(void* key, void* data),
  int (*compar)(void* key1, void* key2, void* data)
);

/** Delete hashmap. */
void hashmap_delete(HashMap* hashmap);

/** Set user data for comparison function. */
void hashmap_setdata(HashMap* hashmap, void* data);

/** Number of key/value pairs in hash map. */
int hashmap_getsize(HashMap* hashmap);

/** Insert key/value pair. */
HashEntry* hashmap_insert(HashMap* hashmap, void* key, void* value);

/** Remove entry. */
void hashmap_remove_entry(HashMap* hashmap, HashEntry* entry);

/** Remove key/value pair. */
void hashmap_remove(HashMap* hashmap, void* key);

/** Remove all key/value pairs. */
void hashmap_clear(HashMap* hashmap);

/** Retrieve value by key. */
void* hashmap_retrieve(HashMap* hashmap, void* key);

/** Key/value enumeration. */
void hashmap_getarrays(HashMap* hashmap, void** keys, void** values);

/** First entry. */
HashEntry* hashmap_first(HashMap* hashmap);

/** Last entry. */
HashEntry* hashmap_last(HashMap* hashmap);



/** Next entry. */
HashEntry* hashentry_next(HashEntry* entry);

/** Previous entry. */
HashEntry* hashentry_prev(HashEntry* entry);

/** Key. */
void* hashentry_getkey(HashEntry* entry);

/** Value. */
void* hashentry_getvalue(HashEntry* entry);

#endif  /* __hashmap_h__ */
