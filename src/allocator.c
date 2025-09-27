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



void* allocator_alloc(Allocator* alc, uint16_t numbytes)
{
  uint32_t numwords;
  uint32_t pool_index;
  Block* ret;

  assert(0 < numbytes && numbytes <= ALLOCATOR_MAXSIZE);

  /* Rounded size and pool index: */
  numwords = (numbytes - 1) / sizeof(Block) + 1;
  pool_index = alc->indices[numwords - 1];

  if (pool_index == 0xff) {
    /* Initialize new pool: */
    assert(alc->npools < ALLOCATOR_MAXPOOLS);
    pool_index = alc->npools;
    alc->freeps[pool_index] = NULL;
    alc->indices[numwords - 1] = pool_index;
    ++alc->npools;
  }
  else if (alc->freeps[pool_index] != NULL) {
    /* Prefer recycled memory: */
    ret = alc->freeps[pool_index];
    alc->freeps[pool_index] = ret->next;
    return ret;
  }

  /* Check if requested memory exceeds current block: */
  if (alc->curblock_freep + numwords > alc->curblock + alc->blocksize) {
    /* New block, double size: */
    alc->blocksize = alc->blocksize > 0 ? alc->blocksize << 1
                                        : (FIRST_BLOCKSIZE / sizeof(Block));
    ret = malloc(alc->blocksize * sizeof(Block));

    /* Prepend new block: */
    ret->next = alc->curblock;
    alc->curblock = ret;

    /* User memory: */
    ++ret;

    /* Set block freeptr: */
    alc->curblock_freep = ret + numwords;
  }
  else {
    /* The memory to return: */
    ret = alc->curblock_freep;

    /* Next free memory in block: */
    alc->curblock_freep += numwords;
  }

  return ret;
}



void allocator_free(Allocator* alc, void* mem, uint16_t numbytes)
{
  uint32_t numwords;
  uint32_t pool_index;

  /* Check input: */
  assert(0 < numbytes && numbytes <= ALLOCATOR_MAXSIZE);
  assert(mem != NULL);

  /* Rounded size: */
  numwords = (numbytes - 1) / sizeof(Block) + 1;

  /* Prefer returning memory to current block: */
  if (alc->curblock_freep == (Block*)mem + numwords) {
    alc->curblock_freep = mem;
    return;
  }

  /* Find pool: */
  pool_index = alc->indices[numwords - 1];
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
