#pragma once

#include "types.h"

/** Initialize hashmap. */
void hashmap_init(HashMap* hashmap);

/** Release hashmap resources. */
void hashmap_destroy(HashMap* hashmap);

/** Insert ridge. */
void hashmap_insert(HashMap* hashmap, int d, u32* verts, u32 ridge);

/** Clear hashmap for reuse. */
void hashmap_clear(HashMap* hashmap);

/** Retrieve ridge by its vertices. */
u32 hashmap_retrieve(HashMap hashmap, int d, u32* verts, Allocator* alc);
