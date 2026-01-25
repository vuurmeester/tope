#include "list.h"



/* Append to list, return next append address. */
List** list_append(List** plist, void* val, Allocator* alc)
{
  for (; *plist != NULL; plist = &(*plist)->next);  // scream to end of list
  *plist = allocator_alloc(alc, sizeof(List));
  (*plist)->val = val;
  return &(*plist)->next;
}



/* Prepend to list. */
void list_prepend(List** plist, void* val, Allocator* alc)
{
  List* newnode = allocator_alloc(alc, sizeof(List));
  newnode->next = *plist;
  newnode->val = val;
  *plist = newnode;
}



/* Pop value from list. */
void* list_pop(List** plist, Allocator* alc)
{
  List* next = (*plist)->next;
  void* val = (*plist)->val;
  allocator_free(alc, *plist, sizeof(List));
  *plist = next;
  return val;
}



void list_free(List** plist, Allocator* alc)
{
  while (*plist) {
    List* next = (*plist)->next;
    allocator_free(alc, *plist, sizeof(List));
    (*plist) = next;
  }
}



int list_len(List* list)
{
  int n = 0;
  for (; list; list = list->next, ++n);
  return n;
}
