#include <string.h>
#include <assert.h>

#include "list.h"

#define MINCAP 2



List list_new()
{
  return (List) {
    0
  };
}



void list_delete(List* lst, Allocator* alc)
{
  ListItem* li = lst->last;
  while (li) {
    ListItem* prev = li->prev;
    allocator_free(alc, li, sizeof(ListItem));
    li = prev;
  }
  memset(lst, 0, sizeof(List));
}



void list_append(List* lst, void* value, Allocator* alc)
{
  ListItem* li = allocator_alloc(alc, sizeof(ListItem));
  li->prev = lst->last;
  li->next = NULL;
  li->value = value;
  lst->last = li;
  if (li->prev) {
    li->prev->next = li;
  } else {
    lst->first = li;
  }
  ++lst->len;
}



void list_remove(List* lst, ListItem* li, Allocator* alc)
{
  if (li->prev) {
    li->prev->next = li->next;
  } else {
    lst->first = li->next;
  }
  if (li->next) {
    li->next->prev = li->prev;
  } else {
    lst->last = li->prev;
  }
  allocator_free(alc, li, sizeof(ListItem));
  --lst->len;
}



void list_remove_value(List* lst, void* value, Allocator* alc)
{
  for (ListItem* li = lst->first; li; li = li->next) {
    if (li->value == value) {
      list_remove(lst, li, alc);
      return;  // done
    }
  }
}



void* list_pop(List* lst, Allocator* alc)
{
  assert(lst->len > 0);
  void* ret = lst->last->value;
  ListItem* prev = lst->last->prev;
  if (prev) {
    prev->next = NULL;
  } else {
    assert(lst->len == 1);
    lst->first = NULL;
  }
  allocator_free(alc, lst->last, sizeof(ListItem));
  lst->last = prev;
  --lst->len;
  return ret;
}
