#pragma once

#include <stdint.h>

typedef uint8_t u8;
typedef uint16_t u16;
typedef uint32_t u32;
typedef uint64_t u64;

#define ALLOCATOR_MAXSIZE 512
#define ALLOCATOR_MAXPOOLS 8

typedef struct _Block Block;
struct _Block {
  Block* next;
};

typedef struct _Allocator {
  u32 blocksize;
  u32 npools;
  Block* freeps[ALLOCATOR_MAXPOOLS];
  u8 indices[ALLOCATOR_MAXSIZE / sizeof(Block)];
  Block* block;
  Block* blockfreep;
} Allocator;

/** Initialize allocator object. */
void allocator_init(Allocator* alc);

/** Allocate number of bytes <= ALLOCATOR_MAXSIZE. */
void* allocator_alloc(Allocator* alc, u16 numbytes);

/** Release memory. */
void allocator_free(Allocator* alc, void* mem, u16 numbytes);

/** Release memory associated with allocator. */
void allocator_destroy(Allocator* alc);
