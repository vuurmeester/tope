#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"

#define FIRST_BLOCKSIZE 4096



void allocator_init(Allocator* alc)
{
  alc->blocksize = FIRST_BLOCKSIZE;
  alc->npools = 0;
  memset(alc->freeps, 0xff, sizeof(alc->freeps));
  alc->blockfreep = 0;
  memset(alc->indices, 0xff, sizeof alc->indices);
  alc->block = malloc(alc->blocksize * sizeof(Block));
}



u32 allocator_alloc(Allocator* alc, uint16_t numbytes)
{
  u32 ret;

  assert(0 < numbytes && numbytes <= ALLOCATOR_MAXSIZE);

  /* Rounded size and pool index: */
  u32 numblocks = (numbytes - 1) / sizeof(Block) + 1;
  u32 pool_index = alc->indices[numblocks - 1];

  if (pool_index == 0xff) {
    /* Initialize new pool: */
    assert(alc->npools < ALLOCATOR_MAXPOOLS);
    pool_index = alc->npools;
    alc->indices[numblocks - 1] = pool_index;
    ++alc->npools;
  }
  else if (alc->freeps[pool_index] != UINT32_MAX) {
    /* Prefer recycled memory: */
    ret = alc->freeps[pool_index];
    alc->freeps[pool_index] = alc->block[ret].next;
    return ret;
  }

  /* Check if requested memory exceeds current block: */
  if (alc->blockfreep + numblocks > alc->blocksize) {
    /* New block, double size: */
    alc->blocksize <<= 1;
    alc->block = realloc(alc->block, alc->blocksize * sizeof(Block));
  }

  /* The memory to return: */
  ret = alc->blockfreep;

#ifndef NDEBUG
  memset(allocator_mem(alc, ret), 0x00, numblocks * sizeof(Block));
#endif

  /* Next free memory in block: */
  alc->blockfreep += numblocks;

  return ret;
}



void allocator_free(Allocator* alc, u32 handle, uint16_t numbytes)
{
  /* Check input: */
  assert(0 < numbytes && numbytes <= ALLOCATOR_MAXSIZE);
  assert(handle != UINT32_MAX);

  /* Rounded size: */
  u32 numblocks = (numbytes - 1) / sizeof(Block) + 1;

#ifndef NDEBUG
  memset(allocator_mem(alc, handle), 0xcd, numblocks * sizeof(Block));
#endif

  /* Prefer returning memory to block freepointer: */
  if (alc->blockfreep == handle + numblocks) {
    alc->blockfreep = handle;
    return;
  }

  /* Find pool: */
  u32 pool_index = alc->indices[numblocks - 1];
  assert(pool_index < alc->npools);

  /* Prepend to pool freelist: */
  alc->block[handle].next = alc->freeps[pool_index];
  alc->freeps[pool_index] = handle;
}



void allocator_destroy(Allocator* alc)
{
  assert(alc != NULL);
  free(alc->block);
}
