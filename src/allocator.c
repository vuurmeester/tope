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
  memset(alc->indices, 0xff, sizeof(alc->indices));
}



void* allocator_alloc(Allocator* alc, u32 size)
{
  u32 size_8;
  u32 pool_index;
  void* ret;

  assert(0 < size && size <= ALLOCATOR_MAXSIZE);

  size_8 = (size - 1) >> 3;
  pool_index = alc->indices[size_8];
  if (pool_index == 0xff) {
    /* This size has not been encountered before. */

    /* Find available pool slot: */
    assert(alc->sizes[ALLOCATOR_NPOOLS - 1] == 0);
    for (pool_index = 0; alc->sizes[pool_index] > 0; ++pool_index)
      ;

    /* Assign pool to this size: */
    alc->indices[size_8] = pool_index;
    alc->sizes[pool_index] = (size_8 + 1) << 3;
  }

  /* Recalibrate size: */
  size = alc->sizes[pool_index];

  /* Prefer recycled memory: */
  if (alc->pool_freeptrs[pool_index] != NULL) {
    ret = alc->pool_freeptrs[pool_index];
    alc->pool_freeptrs[pool_index] = *(void**)ret;
    return ret;
  }

  /* Check if requested memory exceeds current block: */
  if ((char*)alc->curblock_freeptr + size >
      (char*)alc->curblock + alc->blocksize) {
    Block* block;

    /* New block, double size: */
    alc->blocksize = alc->blocksize > 0 ? alc->blocksize << 1 : FIRST_BLOCKSIZE;
    block = malloc(alc->blocksize);

    /* Prepend new block: */
    block->next = alc->curblock;
    alc->curblock = block;

    /* Set block freeptr: */
    alc->curblock_freeptr = block + 1;
  }

  /* The memory to return: */
  ret = alc->curblock_freeptr;

  /* Next free memory in block: */
  alc->curblock_freeptr = (char*)ret + size;

  return ret;
}



void allocator_free(Allocator* alc, void* mem, u32 size)
{
  u32 index;

  /* Check input: */
  assert(0 < size && size <= ALLOCATOR_MAXSIZE);
  assert(mem != NULL);

  /* Where to find the pool freeptr: */
  index = alc->indices[(size - 1) >> 3];
  assert(index < ALLOCATOR_NPOOLS);

  /* Prepend to freelist: */
  *(void**)mem = alc->pool_freeptrs[index];
  alc->pool_freeptrs[index] = mem;
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
