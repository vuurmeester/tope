#pragma once

#include "allocator.h"

/* Non-intrusive singly linked list. */
typedef struct _List List;
struct _List {
  List* next;  /* next in list */
  void* val;  /* generic pointer */
};

/* Append to list, return next append address O(1) or O(n). */
List** list_append(List** plist, void* val, Allocator* alc);

/* Prepend to list O(1). */
void list_prepend(List** plist, void* val, Allocator* alc);

/* Pop value from list O(1). */
void* list_pop(List** plist, Allocator* alc);

/* Free list O(n). */
void list_free(List** plist, Allocator* alc);

/* Number of items in list O(n). */
int list_len(List* list);
