#ifndef HEAP_H
#define HEAP_H

#include <stdbool.h>
#include <pthread.h>

#include "constants.h"
#include "list.h"

#define INIT_HEAP_CAPACITY 8

#define HEAP_GET_LEFT_IDX(idx) (idx * 2)
#define HEAP_GET_RIGHT_IDX(idx) (idx * 2 + 1)
#define HEAP_GET_PARENT_IDX(idx) (idx / 2)

#define HEAP_IDX_NULL 0
#define HEAP_IDX_ROOT 1

typedef struct heap {
  size_t capacity;
  size_t size;
  list_t node_list;
  list_node_t* node_array;
} heap_t;

list_node_t heap_allocate(heap_t* heap, const size_t size);
bool heap_free(heap_t* heap, list_node_t node);
bool heap_init(heap_t* heap, uint8_t* addr, const size_t size);

void FORCE_INLINE
heap_check(heap_t* heap) {
#if DEBUG_HEAP
  list_t* list = &heap->node_list;
  assert(!list->virtual_head->is_free && "Virtual head is dirty.");
  assert(!list->virtual_tail->is_free && "Virtual tail is dirty.");
  list_node_t curr = list->virtual_head->next;
  list_node_t end = list->virtual_tail;
  while (curr != end) {
    assert(curr);
    if (curr->is_free) {
      if (!(curr->idx > 0 && curr->idx < heap->size)) {
        fprintf(stderr, "Idx: %ld; Addr: %p\n", curr->idx, curr);
        assert(0 && "Out of bound!");
      }
    } else {
      assert(curr->idx == 0 && "Invalid null heap idx.");
    }
    curr = curr->next;
  }
  for (int i = 1; i < heap->size; ++i) {
    assert(heap->node_array[i] && "Node array being NULL.");
    assert(heap->node_array[i]->is_free && "Heap not clean.");
    assert(heap->node_array[i]->idx == i && "Idx mismatch.");
  }
#endif
}

#endif