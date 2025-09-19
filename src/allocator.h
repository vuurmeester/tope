#pragma once

#include <stdint.h>

typedef unsigned char u8;
typedef unsigned short u16;
typedef unsigned int u32;

#define ALLOCATOR_MAXSIZE 512
#define ALLOCATOR_MINSIZE 8
#define ALLOCATOR_NPOOLS 8

typedef struct _Block Block;

typedef struct _Allocator {
  u32 blocksize;
  void* pool_freeptrs[ALLOCATOR_NPOOLS];
  u16 sizes[ALLOCATOR_NPOOLS];
  u8 indices[ALLOCATOR_MAXSIZE / ALLOCATOR_MINSIZE]; /* 8, 16, ..., 512 */
  Block* curblock;
  void* curblock_freeptr;
} Allocator;

/** New allocator object. */
void allocator_init(Allocator* alc);

/** Allocate single unit. */
void* allocator_alloc(Allocator* alc, u32 numbytes);

/** Release memory. */
void allocator_free(Allocator* alc, void* mem, u32 numbytes);

/** Release the allocator object, and all memory still associated with it. */
void allocator_clear(Allocator* alc);
