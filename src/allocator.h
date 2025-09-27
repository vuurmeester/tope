#pragma once

#include <stdint.h>

#define ALLOCATOR_MAXSIZE 512
#define ALLOCATOR_MAXPOOLS 8

typedef struct _Block Block;

typedef struct _Allocator {
  uint32_t blocksize;
  uint32_t npools;
  Block* freeps[ALLOCATOR_MAXPOOLS];
  uint8_t indices[ALLOCATOR_MAXSIZE / sizeof(Block*)];
  Block* curblock;
  Block* curblock_freep;
} Allocator;

/** New allocator object. */
void allocator_init(Allocator* alc);

/** Allocate number of bytes <= ALLOCATOR_MAXSIZE. */
void* allocator_alloc(Allocator* alc, uint16_t numbytes);

/** Release memory. */
void allocator_free(Allocator* alc, void* mem, uint16_t numbytes);

/** Release the allocator object, and all memory associated with it. */
void allocator_destroy(Allocator* alc);
