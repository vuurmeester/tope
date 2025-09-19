#pragma once

#include <stdint.h>

typedef unsigned char u8;
typedef unsigned short u16;
typedef unsigned int u32;

#define ALLOCATOR_MAXSIZE 512
#define ALLOCATOR_NPOOLS 8

typedef struct _Block Block;

typedef struct _Allocator {
  u32 blocksize;
  u32 npools;
  Block* pool_freeps[ALLOCATOR_NPOOLS];
  u8 pool_sizes[ALLOCATOR_NPOOLS];
  Block* curblock;
  Block* curblock_freep;
} Allocator;

/** New allocator object. */
void allocator_init(Allocator* alc);

/** Allocate single unit. */
Block* allocator_alloc(Allocator* alc, u16 numbytes);

/** Release memory. */
void allocator_free(Allocator* alc, Block* mem, u16 numbytes);

/** Release the allocator object, and all memory still associated with it. */
void allocator_clear(Allocator* alc);
