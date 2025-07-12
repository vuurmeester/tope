#include <assert.h>
#include <stdlib.h> /* memcpy */
#include <string.h> /* malloc, realloc, free */

#include "allocator.h"

#define MAXSIZE 256
#define FIRST_BLOCKSIZE (8 * MAXSIZE)
#define NPOOLS 8 /* 32, 64, 96, 128, 160, 192, 224, 256 */
#define OUTSIZE(insize, outsize, i)                                            \
  {                                                                            \
    i = (insize - 1) >> 5;                                                     \
    outsize = (i + 1) << 5;                                                    \
  }

typedef struct _Block Block;
struct _Block {
  Block* next;
};

struct _Allocator {
  int blocksize;
  void* freeptrs[NPOOLS];
  Block* firstblock;
  void* freeptr;
};



Allocator* allocator_new(void)
{
  int i;
  Allocator* allocator;

  allocator = malloc(sizeof(Allocator));
  allocator->blocksize = FIRST_BLOCKSIZE;

  for (i = 0; i < NPOOLS; ++i) {
    allocator->freeptrs[i] = NULL;
  }
  allocator->firstblock = malloc(sizeof(Block) + allocator->blocksize);
  allocator->firstblock->next = NULL;
  allocator->freeptr = (char*)allocator->firstblock + sizeof(Block);

  return allocator;
}



void* allocator_alloc(Allocator* allocator, int size)
{
  int outsize;
  int index;
  void* ret;

  if (size <= 0) {
    return NULL;
  }

  OUTSIZE(size, outsize, index);

  if (outsize > MAXSIZE) {
    return malloc(outsize);
  }

  if (allocator->freeptrs[index]) {
    /* Recycle previously free'd memory. */
    ret = allocator->freeptrs[index];
    allocator->freeptrs[index] = *(void**)ret;
    return ret;
  }

  if ((char*)allocator->freeptr > (char*)allocator->firstblock + sizeof(Block) +
                                      allocator->blocksize - outsize) {
    Block* block;

    /* Allocate new block: */
    allocator->blocksize <<= 1;
    block = malloc(sizeof(Block) + allocator->blocksize);

    /* Remember new block: */
    block->next = allocator->firstblock;
    allocator->firstblock = block;

    /* Set freeptr: */
    allocator->freeptr = (char*)block + sizeof(Block);
  }

  /* The memory to return: */
  ret = allocator->freeptr;

  /* Next free memory: */
  allocator->freeptr = (char*)ret + outsize;

  return ret;
}



void allocator_free(Allocator* allocator, void* mem, int size)
{
  int index;
  int outsize;

  if (!mem) {
    return;
  }

  if (size > MAXSIZE) {
    free(mem);
    return;
  }

  /* Where to find the freeptr: */
  OUTSIZE(size, outsize, index);

  /* Prepend to freelist: */
  *(void**)mem = allocator->freeptrs[index];
  allocator->freeptrs[index] = mem;
}



void* allocator_realloc(Allocator* allocator, void* oldmem, int oldsize,
                        int newsize)
{
  int oldindex;
  int newindex;
  int oldoutsize;
  int newoutsize;
  void* newmem;

  if (!oldmem) {
    /* New allocation */
    assert(oldsize == 0);
    return allocator_alloc(allocator, newsize);
  }

  /* Actual sizes and where to find the freeptrs: */
  OUTSIZE(oldsize, oldoutsize, oldindex);
  OUTSIZE(newsize, newoutsize, newindex);

  if (newoutsize == oldoutsize) {
    /* Ideal case, no resizing necessary: */
    return oldmem;
  }

  if (newoutsize > MAXSIZE && oldoutsize > MAXSIZE) {
    /* Both old and new buffers outside of allocator. */
    return realloc(oldmem, newoutsize);
  }

  /* Allocate new memory: */
  newmem = NULL;
  if (newsize > 0) {
    newmem = allocator_alloc(allocator, newoutsize);
  }

  /* Copy contents: */
  memcpy(newmem, oldmem, newsize < oldsize ? newsize : oldsize);

  /* Release old memory: */
  allocator_free(allocator, oldmem, oldsize);

  return newmem;
}



void allocator_clear(Allocator* allocator)
{
  int i;

  /* Clear pools: */
  for (i = 0; i < NPOOLS; ++i) {
    allocator->freeptrs[i] = NULL;
  }

  /* Delete allocated blocks: */
  while (allocator->firstblock->next) {
    Block* next;
    next = allocator->firstblock->next;
    free(allocator->firstblock);
    allocator->firstblock = next;
  }
  allocator->blocksize = FIRST_BLOCKSIZE;
  allocator->freeptr = (char*)allocator->firstblock + sizeof(Block);
}



void allocator_delete(Allocator* allocator)
{
  if (!allocator) {
    return;
  }

  /* Clear the allocator: */
  allocator_clear(allocator);

  /* Free first block: */
  free(allocator->firstblock);

  /* Delete the allocator object itself: */
  free(allocator);
}
