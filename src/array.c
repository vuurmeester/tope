#include "array.h"



Array array_new(void)
{
  Array arr;
  arr.len = 0;
  arr.cap = 0;
  arr.values = NULL;
  return arr;
}



void array_delete(Array* arr, Allocator* alc)
{
  allocator_free(alc, arr->values, arr->cap * sizeof(void*));
  memset(arr, 0, sizeof(Array));
}



void array_append(Array* arr, void* value, Allocator* alc)
{
  if (arr->len == arr->cap) {
    int newcap = arr->cap > 0 ? 2 * arr->cap : 4;
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
  for (int i = 0; i < arr.len; ++i) {
    if (arr.values[i] == value) {
      return i;
    }
  }
  return -1;
}



void* array_pop(Array* arr)
{
  --arr->len;
  return arr->values[arr->len];
}
