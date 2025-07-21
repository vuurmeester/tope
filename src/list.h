#pragma once

#include <stdlib.h>

#include "allocator.h"

typedef struct _ListItem ListItem;

struct _ListItem {
  ListItem* prev;
  ListItem* next;
  void* value;
};

typedef struct _List {
  ListItem* first;
  ListItem* last;
  int len;
} List;



List list_new(void);
void list_delete(List* lst, Allocator* alc);
void list_append(List* lst, void* value, Allocator* alc);
void list_remove(List* lst, ListItem* item, Allocator* alc);
void list_remove_value(List* lst, void* value, Allocator* alc);
void* list_pop(List* lst, Allocator* alc);
