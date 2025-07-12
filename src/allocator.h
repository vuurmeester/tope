#pragma once

typedef struct _Allocator Allocator;

/** New allocator object. */
Allocator* allocator_new(void);

/** Allocate single unit. */
void* allocator_alloc(Allocator* allocator, int numbytes);

/** Release unit. */
void allocator_free(Allocator* allocator, void* mem, int numbytes);

/** Resize memory block. */
void* allocator_realloc(Allocator* allocator, void* oldmem, int oldsize, int newsize);

/** Release all units. */
void allocator_clear(Allocator* allocator);

/** Release the allocator object, and all memory still associated with it. */
void allocator_delete(Allocator* allocator);
