#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"

#define BLOCKSIZE_STEP 512



void allocator_init(Allocator* alc)
{
  alc->blocksize = 0;
  alc->npools = 0;
  memset(alc->freeps , 0x00, sizeof alc->freeps );
  memset(alc->indices, 0xff, sizeof alc->indices);
  alc->block = NULL;
  alc->blockfreep = NULL;
  alc->used = 0;
}



void* allocator_alloc(Allocator* alc, u16 numbytes)
{
  Block* mem;

#ifdef USE_MALLOC
  mem = malloc(numbytes);
  memset(mem, 0x00, numbytes);
  return mem;
#endif

  assert(0 < numbytes && numbytes <= ALLOCATOR_MAXSIZE);

  /* Rounded size and pool index: */
  u32 numblocks = (numbytes - 1) / sizeof(Block) + 1;
  u8 pool_index = alc->indices[numblocks - 1];

  if (pool_index == 0xff) {
    /* Initialize new pool: */
    assert(alc->npools < ALLOCATOR_MAXPOOLS);
    pool_index = alc->npools;
    alc->indices[numblocks - 1] = pool_index;
    ++alc->npools;
  }
  else if (alc->freeps[pool_index] != NULL) {
    /* Prefer recycled memory: */
    mem = alc->freeps[pool_index];
    alc->freeps[pool_index] = mem->next;
    goto _return;
  }

  /* Check if requested memory exceeds current block: */
  if (alc->blockfreep + numblocks > alc->block + alc->blocksize) {
    /* Increase block size : */
    alc->blocksize += BLOCKSIZE_STEP;
    Block* newblock = malloc(alc->blocksize * sizeof(Block));
    alc->used += alc->blocksize * sizeof(Block);
    newblock->next = alc->block;
    alc->block = newblock;
    alc->blockfreep = alc->block + 1;  /* step over 'next' pointer */
  }

  /* The memory to return: */
  mem = alc->blockfreep;

  /* Next free memory in block: */
  alc->blockfreep += numblocks;

_return:
  memset(mem, 0x00, numblocks * sizeof(Block));
  return mem;
}



void allocator_free(Allocator* alc, void* mem, u16 numbytes)
{
#ifdef USE_MALLOC
  free(mem);
  return;
#endif

  /* Check input: */
  assert(0 < numbytes && numbytes <= ALLOCATOR_MAXSIZE);
  assert(mem != NULL);

  /* Rounded size: */
  u32 numblocks = (numbytes - 1) / sizeof(Block) + 1;

#ifndef NDEBUG
  memset(mem, 0xcd, numblocks * sizeof(Block));
#endif

  /* Memory at end of block?: */
  if (alc->blockfreep == (Block*)mem + numblocks) {
    alc->blockfreep = mem;
    return;
  }

  /* Find pool: */
  u32 pool_index = alc->indices[numblocks - 1];
  assert(pool_index < alc->npools);

  /* Prepend to pool freelist: */
  ((Block*)mem)->next = alc->freeps[pool_index];
  alc->freeps[pool_index] = mem;
}



void allocator_destroy(Allocator* alc)
{
  assert(alc != NULL);
  while (alc->block != NULL) {
    Block* next = alc->block->next;
    free(alc->block);
    alc->block = next;
  }
  memset(alc, 0x00, sizeof *alc);
}
