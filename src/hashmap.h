#pragma once

#include "types.h"

/** Create hashmap. */
HashMap* hashmap_new(void);

/** Delete hashmap. */
void hashmap_delete(HashMap* hashmap);

/** Insert ridge. */
void hashmap_insert(HashMap* hashmap, int d, Ridge* ridge);

/** Clear hashmap for reuse. */
void hashmap_clear(HashMap* hashmap);

/** Retrieve ridge by its vertices. */
Ridge* hashmap_retrieve(HashMap* hashmap, int d, polytoop_Vertex** verts);
