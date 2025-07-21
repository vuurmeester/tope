#pragma once

#include <stdint.h>

#define NPOOLS 8

typedef struct _Block Block;

typedef struct _Allocator {
  int blocksize;
  void* freeptrs[NPOOLS];
  uint16_t sizes[NPOOLS];
  uint8_t indices[64];  // 8, 16, ..., 512
  Block* firstblock;
  void* freeptr;
} Allocator;

/** New allocator object. */
void allocator_init(Allocator* alc);

/** Allocate single unit. */
void* allocator_alloc(Allocator* alc, int numbytes);

/** Release memory. */
void allocator_free(Allocator* alc, void* mem, int numbytes);

/** Release the allocator object, and all memory still associated with it. */
void allocator_clear(Allocator* alc);
