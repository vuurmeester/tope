#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"

#define FIRST_BLOCKSIZE 4096


void allocator_init(Allocator* alc)
{
  memset(alc, 0, sizeof *alc);

  /* Invalidate indices: */
  memset(alc->indices, 0xff, sizeof alc->indices);
}



u32 allocator_alloc(Allocator* alc, uint16_t numbytes)
{
  u32 ret;

  assert(0 < numbytes && numbytes <= ALLOCATOR_MAXSIZE);

  /* Rounded size and pool index: */
  u32 numwords = (numbytes - 1) / sizeof(Block) + 1;
  u32 pool_index = alc->indices[numwords - 1];

  if (pool_index == 0xff) {
    /* Initialize new pool: */
    assert(alc->npools < ALLOCATOR_MAXPOOLS);
    pool_index = alc->npools;
    alc->freeps[pool_index] = -1;
    alc->indices[numwords - 1] = pool_index;
    ++alc->npools;
  }
  else if (alc->freeps[pool_index] != UINT32_MAX) {
    /* Prefer recycled memory: */
    ret = alc->freeps[pool_index];
    alc->freeps[pool_index] = alc->block[ret].next;
    return ret;
  }

  /* Check if requested memory exceeds current block: */
  if (alc->blockfreep + numwords > alc->blocksize) {
    /* New block, double size: */
    alc->blocksize = alc->blocksize > 0 ? alc->blocksize << 1
                                        : (FIRST_BLOCKSIZE / sizeof(Block));
    alc->block = realloc(alc->block, alc->blocksize * sizeof(Block));
  }

  /* The memory to return: */
  ret = alc->blockfreep;

#ifndef NDEBUG
  memset(allocator_mem(alc, ret), 0x00, numwords * sizeof(Block));
#endif

  /* Next free memory in block: */
  alc->blockfreep += numwords;

  return ret;
}



void allocator_free(Allocator* alc, u32 handle, uint16_t numbytes)
{
  /* Check input: */
  assert(0 < numbytes && numbytes <= ALLOCATOR_MAXSIZE);
  assert(handle != UINT32_MAX);

  /* Rounded size: */
  u32 numwords = (numbytes - 1) / sizeof(Block) + 1;

#ifndef NDEBUG
  memset(allocator_mem(alc, handle), 0xff, numwords * sizeof(Block));
#endif

  /* Prefer returning memory to current block: */
  if (alc->blockfreep == handle + numwords) {
    alc->blockfreep = handle;
    return;
  }

  /* Find pool: */
  u32 pool_index = alc->indices[numwords - 1];
  assert(pool_index < alc->npools);

  /* Prepend to freelist: */
  alc->block[handle].next = alc->freeps[pool_index];
  alc->freeps[pool_index] = handle;
}



void allocator_destroy(Allocator* alc)
{
  assert(alc != NULL);
  free(alc->block);
}
