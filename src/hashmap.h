#pragma once

#include "types.h"

/** Initialize hashmap. */
void hashmap_init(HashMap* hashmap);

/** Release hashmap. */
void hashmap_destroy(HashMap* hashmap);

/** Insert ridge. */
void hashmap_insert(HashMap* hashmap, int d, Ridge* ridge);

/** Clear hashmap for reuse. */
void hashmap_clear(HashMap* hashmap);

/** Retrieve ridge by its vertices. */
Ridge* hashmap_retrieve(HashMap hashmap, int d, polytoop_Vertex** verts);
