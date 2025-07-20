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

static uint16_t sizes[NPOOLS] = {
  8, 16, 24, 32, 48, 64, 96, 128, 192, 256, 384, 512,
};

static uint8_t indices[] = {
  0, // 8
  1, // 16
  2, // 24
  3, // 32
  4, 4, // 40, 48
  5, 5, // 56, 64
  6, 6, 6, 6, // 72, 80, 88, 96
  7, 7, 7, 7, // 104, 112, 120, 128
  8, 8, 8, 8, 8, 8, 8, 8, // 136, ..., 192
  9, 9, 9, 9, 9, 9, 9, 9, // 200, ..., 256
  10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
  11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11, 11,
};



void allocator_init(Allocator* alc)
{
  memset(alc, 0, sizeof(Allocator));
  alc->freeptr = (void*) -1;
}



void* allocator_alloc(Allocator* alc, int size)
{
  if (size <= 0) {
    return NULL;
  }

  if (size > MAXSIZE) {
    /* out of band: */
    return malloc(size);
  }

  int index = indices[(size - 1) >> 3];
  int outsize = sizes[index];

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
  if (!mem || size <= 0) {
    return;
  }

  if (size > MAXSIZE) {
    free(mem);
    return;
  }

  /* Where to find the freeptr: */
  int index = indices[(size - 1) >> 3];

  /* Prepend to freelist: */
  *(void**)mem = alc->freeptrs[index];
  alc->freeptrs[index] = mem;
}



void* allocator_realloc(Allocator* allocator, void* oldmem, int oldsize,
                        int newsize)
{
  if (!oldmem) {
    /* New allocation */
    assert(oldsize == 0);
    return allocator_alloc(allocator, newsize);
  }

  int newoutsize = newsize;
  if (newsize <= MAXSIZE) {
    int newindex = indices[(newsize - 1) >> 3];
    newoutsize = sizes[newindex];
  }

  int oldoutsize = oldsize;
  if (oldsize <= MAXSIZE) {
    int oldindex = indices[(oldsize - 1) >> 3];
    oldoutsize = sizes[oldindex];
  }

  if (newoutsize == oldoutsize) {
    /* No resizing necessary: */
    return oldmem;
  }

  /* Allocate new memory: */
  void* newmem = NULL;
  if (newsize > 0) {
    newmem = allocator_alloc(allocator, newoutsize);
    memcpy(newmem, oldmem, newsize < oldsize ? newsize : oldsize);
  }

  /* Release old memory: */
  allocator_free(allocator, oldmem, oldsize);

  return newmem;
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
  alc->freeptr = (void*) -1;
}
