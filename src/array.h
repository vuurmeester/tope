#pragma once

#include <stdlib.h>

#include "allocator.h"

typedef struct _Array {
  int len;
  int cap;
  void** values;
} Array;



Array array_new(void);
void array_delete(Array* arr, Allocator* alc);
void array_append(Array* arr, void* value, Allocator* alc);
void array_remove(Array* arr, int i);
int array_find(Array arr, void* value);
void* array_pop(Array* arr);
