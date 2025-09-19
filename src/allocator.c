#include <assert.h>
#include <stdint.h>
#include <stdlib.h> /* memcpy */
#include <string.h> /* malloc, realloc, free */

#include "allocator.h"

#define FIRST_BLOCKSIZE 4096

struct _Block {
  Block* next;
};



void allocator_init(Allocator* alc)
{
  memset(alc, 0, sizeof *alc);
  memset(alc->pool_sizes, 0xff, sizeof alc->pool_sizes);
}



Block* allocator_alloc(Allocator* alc, u16 numbytes)
{
  u32 numwords;
  u32 pool_index;
  Block* ret;

  assert(0 < numbytes && numbytes <= ALLOCATOR_MAXSIZE);

  numwords = ((numbytes - 1) / sizeof(Block)) + 1;
  for (pool_index = 0; alc->pool_sizes[pool_index] < numwords; ++pool_index)
    ;
  assert(pool_index < ALLOCATOR_NPOOLS);
  if (numwords < alc->pool_sizes[pool_index]) {
    u32 index;

    /* Shift arrays up: */
    for (index = alc->npools; index > pool_index; --index) {
      alc->pool_sizes[index] = alc->pool_sizes[index - 1];
      alc->pool_freeps[index] = alc->pool_freeps[index - 1];
    }

    /* Insert new size: */
    alc->pool_sizes[pool_index] = numwords;
    alc->pool_freeps[pool_index] = NULL;
    ++alc->npools;
  } else if (alc->pool_freeps[pool_index] != NULL) {
    /* Prefer recycled memory: */
    ret = alc->pool_freeps[pool_index];
    alc->pool_freeps[pool_index] = ret->next;
    return ret;
  }

  /* Check if requested memory exceeds current block: */
  if (alc->curblock_freep + numwords > alc->curblock + alc->blocksize) {
    Block* block;

    /* New block, double size: */
    alc->blocksize = alc->blocksize > 0 ? alc->blocksize << 1
                                        : (FIRST_BLOCKSIZE / sizeof(Block));
    block = malloc(alc->blocksize * sizeof(Block));

    /* Prepend new block: */
    block->next = alc->curblock;
    alc->curblock = block;

    /* Set block freeptr: */
    alc->curblock_freep = block + 1;
  }

  /* The memory to return: */
  ret = alc->curblock_freep;

  /* Next free memory in block: */
  alc->curblock_freep = ret + numwords;

  return ret;
}



void allocator_free(Allocator* alc, Block* mem, u16 numbytes)
{
  u32 numwords;
  u32 pool_index;

  /* Check input: */
  assert(0 < numbytes && numbytes <= ALLOCATOR_MAXSIZE);
  assert(mem != NULL);

  numwords = ((numbytes - 1) / sizeof(Block)) + 1;
  for (pool_index = 0; alc->pool_sizes[pool_index] != numwords; ++pool_index)
    ;
  assert(pool_index < ALLOCATOR_NPOOLS);

  /* Prepend to freelist: */
  mem->next = alc->pool_freeps[pool_index];
  alc->pool_freeps[pool_index] = mem;
}



void allocator_clear(Allocator* alc)
{
  assert(alc != NULL);
  while (alc->curblock) {
    Block* next;

    next = alc->curblock->next;
    free(alc->curblock);
    alc->curblock = next;
  }
  allocator_init(alc);
}
