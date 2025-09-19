#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"

#define FIRST_BLOCKSIZE 4096

struct _Block {
  Block* next;
};


void allocator_init(Allocator* alc)
{
  memset(alc, 0, sizeof *alc);

  /* Invalidate indices: */
  memset(alc->indices, 0xff, sizeof alc->indices);
}



void* allocator_alloc(Allocator* alc, u16 numbytes)
{
  u32 numwords;
  u32 pool_index;
  Block* ret;

  assert(0 < numbytes && numbytes <= ALLOCATOR_MAXSIZE);

  /* Round size to wordcount and find associated pool: */
  numwords = (numbytes - 1) / sizeof(Block);
  pool_index = alc->indices[numwords];

  if (pool_index == 0xff) {
    /* Append new pool: */
    assert(alc->npools < ALLOCATOR_MAXPOOLS);
    pool_index = alc->npools;
    alc->freeps[pool_index] = NULL;
    alc->indices[numwords] = pool_index;
    ++alc->npools;
  } else if (alc->freeps[pool_index] != NULL) {
    /* Prefer recycled memory: */
    ret = alc->freeps[pool_index];
    alc->freeps[pool_index] = ret->next;
    return ret;
  }

  ++numwords;

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



void allocator_free(Allocator* alc, void* mem, u16 numbytes)
{
  u32 numwords;
  u32 pool_index;

  /* Check input: */
  assert(0 < numbytes && numbytes <= ALLOCATOR_MAXSIZE);
  assert(mem != NULL);

  /* Find pool: */
  numwords = (numbytes - 1) / sizeof(Block);
  pool_index = alc->indices[numwords];
  assert(pool_index < alc->npools);

  /* Prepend to freelist: */
  ((Block*)mem)->next = alc->freeps[pool_index];
  alc->freeps[pool_index] = mem;
}



void allocator_destroy(Allocator* alc)
{
  assert(alc != NULL);
  while (alc->curblock) {
    Block* next;

    next = alc->curblock->next;
    free(alc->curblock);
    alc->curblock = next;
  }
}
