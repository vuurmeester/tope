#include <string.h>

#include "array.h"

#define MINCAP 2



int find(int n, void** values, void* value)
{
  for (int i = 0; i < n; ++i) {
    if (values[i] == value) {
      return i;
    }
  }
  return -1;
}



Array array_new(int cap, Allocator* alc)
{
  int newcap = MINCAP;
  while (newcap < cap) {
    newcap <<= 1;
  }
  return (Array) {
    .values = allocator_alloc(alc, newcap * sizeof(void*)),
    .cap = newcap,
    .len = 0
  };
}



void array_delete(Array* arr, Allocator* alc)
{
  allocator_free(alc, arr->values, arr->cap * sizeof(void*));
  memset(arr, 0, sizeof(Array));
}



void array_append(Array* arr, void* value, Allocator* alc)
{
  if (arr->len == arr->cap) {
    int newcap = arr->cap << 1;
    arr->values = allocator_realloc(alc, arr->values, arr->cap * sizeof(void*),
                                    newcap * sizeof(void*));
    arr->cap = newcap;
  }
  arr->values[arr->len] = value;
  ++arr->len;
}



void array_remove(Array* arr, int i)
{
  --arr->len;
  if (i < arr->len) {
    arr->values[i] = arr->values[arr->len];
  }
}



int array_find(Array arr, void* value)
{
  return find(arr.len, arr.values, value);
}



void* array_pop(Array* arr)
{
  --arr->len;
  return arr->values[arr->len];
}
