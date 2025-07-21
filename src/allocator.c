#include <assert.h>
#include <stdlib.h> /* memcpy */
#include <string.h> /* malloc, realloc, free */
#include <stdint.h>

#include "allocator.h"

#define MAXSIZE 512
#define FIRST_BLOCKSIZE 4096

struct _Block {
  Block* next;
};



void allocator_init(Allocator* alc)
{
  memset(alc, 0, sizeof(Allocator));
  memset(alc->indices, 0xff, sizeof(alc->indices));
  alc->freeptr = (void*) -1;
}



void* allocator_alloc(Allocator* alc, int size)
{
  // Must be between zero and MAXSIZE:
  assert(0 < size && size <= MAXSIZE);

  int size_8 = (size - 1) >> 3;
  int index = alc->indices[size_8];
  if (index == 0xff) {
    assert(alc->sizes[NPOOLS - 1] == 0);  // at least one slot free
    for (index = 0; alc->sizes[index] > 0; ++index);
    alc->indices[size_8] = index;
    alc->sizes[index] = (size_8 + 1) << 3;
  }
  int outsize = alc->sizes[index];

  void* ret = NULL;
  if (alc->freeptrs[index]) {
    /* Recycle previously freed memory. */
    ret = alc->freeptrs[index];
    alc->freeptrs[index] = *(void**)ret;
    return ret;
  }

  if ((char*)alc->freeptr >
      (char*)alc->firstblock + alc->blocksize - outsize) {
    /* Allocate new block: */
    alc->blocksize = alc->blocksize ? alc->blocksize << 1 : FIRST_BLOCKSIZE;
    Block* block = malloc(alc->blocksize);

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



void allocator_free(Allocator* alc, void* mem, int size)
{
  // Must be between zero and MAXSIZE:
  assert(0 <= size && size <= MAXSIZE);

  if (!mem || size == 0) {
    return;
  }

  /* Where to find the freeptrs: */
  int index = alc->indices[(size - 1) >> 3];
  assert(index >= 0);

  /* Prepend to freelist: */
  *(void**)mem = alc->freeptrs[index];
  alc->freeptrs[index] = mem;
}



void allocator_clear(Allocator* alc)
{
  assert(alc);
  while (alc->firstblock) {
    Block* next = alc->firstblock->next;
    free(alc->firstblock);
    alc->firstblock = next;
  }
  memset(alc, 0, sizeof(Allocator));
  memset(alc->indices, 0xff, sizeof(alc->indices));
  alc->freeptr = (void*) -1;
}
