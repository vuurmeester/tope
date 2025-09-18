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
  memset(alc, 0, sizeof(Allocator));
  memset(alc->indices, 0xff, sizeof(alc->indices));
  alc->freeptr = (void*)-1;
}



void* allocator_alloc(Allocator* alc, u32 size)
{
  u32 size_8;
  u32 index;
  u32 outsize;
  void* ret;
  Block* block;

  /* Must be between zero and MAXSIZE: */
  assert(0 < size && size <= MAXSIZE);

  size_8 = (size - 1) >> 3;
  index = alc->indices[size_8];
  if (index == 0xff) {
    assert(alc->sizes[NPOOLS - 1] == 0); /* at least one slot free */
    for (index = 0; alc->sizes[index] > 0; ++index)
      ;
    alc->indices[size_8] = index;
    alc->sizes[index] = (size_8 + 1) << 3;
  }
  outsize = alc->sizes[index];

  ret = NULL;
  if (alc->freeptrs[index]) {
    /* Recycle previously freed memory. */
    ret = alc->freeptrs[index];
    alc->freeptrs[index] = *(void**)ret;
    return ret;
  }

  if ((char*)alc->freeptr > (char*)alc->firstblock + alc->blocksize - outsize) {
    /* Allocate new block: */
    alc->blocksize = alc->blocksize ? alc->blocksize << 1 : FIRST_BLOCKSIZE;
    block = malloc(alc->blocksize);

    /* Prepend new block: */
    block->next = alc->firstblock;
    alc->firstblock = block;

    /* Set freeptr: */
    alc->freeptr = block + 1;
  }

  /* The memory to return: */
  ret = alc->freeptr;

  /* Next free memory: */
  alc->freeptr = (char*)ret + outsize;

  return ret;
}



void allocator_free(Allocator* alc, void* mem, u32 size)
{
  u32 index;

  /* Check sensibility: */
  assert(0 < size && size <= MAXSIZE);
  assert(mem != NULL);

  /* Where to find the freeptr: */
  index = alc->indices[(size - 1) >> 3];
  assert(index < NPOOLS);

  /* Prepend to freelist: */
  *(void**)mem = alc->freeptrs[index];
  alc->freeptrs[index] = mem;
}



void allocator_clear(Allocator* alc)
{
  Block* next;

  assert(alc != NULL);
  while (alc->firstblock) {
    next = alc->firstblock->next;
    free(alc->firstblock);
    alc->firstblock = next;
  }
  allocator_init(alc);
}
