#pragma once

#define NPOOLS 12 /* 8, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512 */

typedef struct _Block Block;

typedef struct _Allocator {
  int blocksize;
  void* freeptrs[NPOOLS];
  Block* firstblock;
  void* freeptr;
} Allocator;

/** New allocator object. */
void allocator_init(Allocator* alc);

/** Allocate single unit. */
void* allocator_alloc(Allocator* alc, int numbytes);

/** Release memory. */
void allocator_free(Allocator* alc, void* mem, int numbytes);

/** Resize memory block. */
void* allocator_realloc(Allocator* alc, void* oldmem, int oldsize, int newsize);

/** Release the allocator object, and all memory still associated with it. */
void allocator_clear(Allocator* alc);
