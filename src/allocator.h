#pragma once

#include <stdint.h>

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;

#define ALLOCATOR_MAXSIZE 512
#define ALLOCATOR_MAXPOOLS 8

typedef struct _Block {
  u32 next;
  u32 _;  /* 8 byte alignment */
} Block;

typedef struct _Allocator {
  u32 blocksize;
  u32 npools;
  u32 freeps[ALLOCATOR_MAXPOOLS];
  u32 blockfreep;
  u8 indices[ALLOCATOR_MAXSIZE / sizeof(Block)];
  Block* block;
} Allocator;

/** New allocator object. */
void allocator_init(Allocator* alc);

/** Allocate number of bytes <= ALLOCATOR_MAXSIZE. */
u32 allocator_alloc(Allocator* alc, u16 numbytes);

/** The actual memory */
#define allocator_mem(alc, handle) ((void*)((alc)->block + handle))

/** Release memory. */
void allocator_free(Allocator* alc, u32 handle, u16 numbytes);

/** Release the allocator object, and all memory associated with it. */
void allocator_destroy(Allocator* alc);
