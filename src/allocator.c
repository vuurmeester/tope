#include <assert.h>
#include <stdlib.h>
#include <string.h>

#include "allocator.h"

#define FIRST_BLOCKSIZE 4096

struct _Block {
  u32 next;
  u32 _;
};


void allocator_init(Allocator* alc)
{
  memset(alc, 0, sizeof *alc);

  /* Invalidate indices: */
  memset(alc->indices, 0xff, sizeof alc->indices);
}



u32 allocator_alloc(Allocator* alc, uint16_t numbytes)
{
  uint32_t numwords;
  uint32_t pool_index;
  u32 ret;

  assert(0 < numbytes && numbytes <= ALLOCATOR_MAXSIZE);

  /* Rounded size and pool index: */
  numwords = (numbytes - 1) / sizeof(void*) + 1;
  pool_index = alc->indices[numwords - 1];

  if (pool_index == 0xff) {
    /* Initialize new pool: */
    assert(alc->npools < ALLOCATOR_MAXPOOLS);
    pool_index = alc->npools;
    alc->freeps[pool_index] = -1;
    alc->indices[numwords - 1] = pool_index;
    ++alc->npools;
  }
  else if (alc->freeps[pool_index] != (u32)-1) {
    /* Prefer recycled memory: */
    ret = alc->freeps[pool_index];
    alc->freeps[pool_index] = alc->curblock[ret].next;
    return ret;
  }

  /* Check if requested memory exceeds current block: */
  if (alc->curblock_freep + numwords > alc->blocksize) {
    /* New block, double size: */
    alc->blocksize = alc->blocksize > 0 ? alc->blocksize << 1
                                        : (FIRST_BLOCKSIZE / sizeof(Block));
    alc->curblock = realloc(alc->curblock, alc->blocksize * sizeof(Block));
  }

  /* The memory to return: */
  ret = alc->curblock_freep;

  /* Next free memory in block: */
  alc->curblock_freep += numwords;

  return ret;
}



void* allocator_mem(Allocator* alc, u32 handle)
{
  return &alc->curblock[handle];
}



void allocator_free(Allocator* alc, u32 handle, uint16_t numbytes)
{
  uint32_t numwords;
  uint32_t pool_index;

  /* Check input: */
  assert(0 < numbytes && numbytes <= ALLOCATOR_MAXSIZE);
  assert(handle != (u32)-1);

  /* Rounded size: */
  numwords = (numbytes - 1) / sizeof(Block) + 1;

  memset(allocator_mem(alc, handle), 0xff, numwords * sizeof(Block));

  /* Prefer returning memory to current block: */
  if (alc->curblock_freep == handle + numwords) {
    alc->curblock_freep = handle;
    return;
  }

  /* Find pool: */
  pool_index = alc->indices[numwords - 1];
  assert(pool_index < alc->npools);

  /* Prepend to freelist: */
  alc->curblock[handle].next = alc->freeps[pool_index];
  alc->freeps[pool_index] = handle;
}



void allocator_destroy(Allocator* alc)
{
  assert(alc != NULL);
  free(alc->curblock);
}
