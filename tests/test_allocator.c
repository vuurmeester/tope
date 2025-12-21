#include <stdlib.h>
#include <string.h>

#include "../src/allocator.h"



int main()
{
  int ntests = 100;
  for (int itest = 0; itest < ntests; ++itest) {
    Allocator alc;
    allocator_init(&alc);
    u16 sizes[ALLOCATOR_MAXPOOLS];
    for (int i = 0; i < ALLOCATOR_MAXPOOLS; ++i) {
      sizes[i] = rand() % ALLOCATOR_MAXSIZE + 1;
    }
    int nalloc = 1000;
    char** allocs = calloc(nalloc, sizeof(char*));
    u16* allocsizes = calloc(nalloc, sizeof(u16));
    for (int i = 0; i < nalloc; ++i) {
      int isize = rand() % ALLOCATOR_MAXPOOLS;
      allocsizes[i] = sizes[isize];
      allocs[i] = allocator_alloc(&alc, allocsizes[i]);
      memset(allocs[i], 'g', allocsizes[i]);

      // Some random deallocation:
      int ifree = rand() % nalloc;
      if (allocsizes[ifree] > 0) {
        for (int j = 0; j < allocsizes[ifree]; ++j) {
          if (allocs[ifree][j] != 'g') {
            return -1;
          }
        }
        allocator_free(&alc, allocs[ifree], allocsizes[ifree]);
        allocs[ifree] = NULL;
        allocsizes[ifree] = 0;
      }
    }

    free(allocsizes);
    free(allocs);

    allocator_destroy(&alc);
  }
  return 0;
}
