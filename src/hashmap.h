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
void hashmap_insert(HashMap* hashmap, Ridge* ridge);

/** Remove all key/value pairs. */
void hashmap_clear(HashMap* hashmap);

/** Retrieve value by key. */
Ridge* hashmap_retrieve(HashMap* hashmap, polytoop_Vertex** verts);
