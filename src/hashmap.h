#include "allocator.h"

#include "types.h"

/** Create hashmap. */
HashMap* hashmap_new(
  unsigned(*hashfunc)(void* key, void* data),
  int (*compar)(void* key1, void* key2, void* data)
);

/** Delete hashmap. */
void hashmap_delete(HashMap* hashmap);

/** Set user data for comparison function. */
void hashmap_setdata(HashMap* hashmap, void* data);

/** Insert key/value pair. */
HashEntry* hashmap_insert(HashMap* hashmap, void* key, void* value);

/** Remove all key/value pairs. */
void hashmap_clear(HashMap* hashmap);

/** Retrieve value by key. */
void* hashmap_retrieve(HashMap* hashmap, void* key);
