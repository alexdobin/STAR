#include <assert.h>
#include <stdio.h>

#include "heap.h"
#include "list.h"
#include <stdlib.h>

/*
Since these functions are recursive, they cannot be inlined. 
Put them in the source code instead.
*/

static void heap_heapify_down(heap_t* heap, const size_t idx);
static void heap_heapify_up(heap_t* heap, const size_t idx);

void LOCAL_HELPER
heap_set_idx(heap_t* heap, const size_t idx, list_node_t node) {
  heap->node_array[idx] = node;
  node->idx = idx;
}

size_t LOCAL_HELPER
heap_get_num_elements(heap_t* heap) {
  return heap->size - 1;
}

size_t LOCAL_HELPER 
heap_compare_GT(heap_t* heap, const size_t l, const size_t r) {
  return (heap->node_array[l]->size > heap->node_array[r]->size) ? l : r;
}

void LOCAL_HELPER 
heap_swap(heap_t* heap, const size_t a, const size_t b) {
  heap->node_array[a]->idx = b;
  heap->node_array[b]->idx = a;
  list_node_t temp = heap->node_array[a];
  heap->node_array[a] = heap->node_array[b];
  heap->node_array[b] = temp;
}

void LOCAL_HELPER
heap_remove_idx(heap_t* heap, const size_t idx) {
  if (!(idx > 0 && idx < heap->size && "Heap remove invalid idx.")) {
    fprintf(stderr, "idx: %ld\n", idx);
    assert(0);
  }

  const size_t heap_elements = heap_get_num_elements(heap);
  if (heap_elements == idx) {
    heap->size -= 1;
    heap->node_array[idx]->idx = HEAP_IDX_NULL;
    return;
  }

  heap->node_array[idx]->idx = HEAP_IDX_NULL;
  heap_set_idx(heap, idx, heap->node_array[heap->size - 1]);
  heap->size -= 1;
  heap_heapify_down(heap, idx);
}

list_node_t LOCAL_HELPER
heap_get_root(heap_t* heap) {
  if (heap->size <= HEAP_IDX_ROOT) {
    return NULL;
  }
  return heap->node_array[HEAP_IDX_ROOT];
}

// Insert a new node to the heap with heap property maintained.
bool LOCAL_HELPER
heap_insert(heap_t* heap, list_node_t node) {
  if (heap->size == heap->capacity) {
    // Needs to expand.
    list_node_t* new_array = realloc(
      heap->node_array, 
      sizeof(list_node_base_t) * heap->capacity * 2
    );
    if (new_array == NULL) {
      fprintf(
        stderr, 
        "Error: failed to allocate more memory for internal heap.\n"
      );
      return false;
    }
    heap->node_array = new_array;
    heap->capacity *= 2;
  }

  heap_set_idx(heap, heap->size, node);
  node->is_free = true;
  heap->size += 1;
  heap_heapify_up(heap, node->idx);

  return true;
}

/*
Heapify downward. This happens when we replace a node by a smaller number.
*/
static void heap_heapify_down(heap_t* heap, const size_t idx) {
  const size_t left_idx = HEAP_GET_LEFT_IDX(idx);
  const size_t right_idx = HEAP_GET_RIGHT_IDX(idx);
  size_t largest_idx = idx;
  if (left_idx < heap->size) {
    largest_idx = heap_compare_GT(heap, largest_idx, left_idx);
  }
  if (right_idx < heap->size) {
    largest_idx = heap_compare_GT(heap, largest_idx, right_idx);
  }
  if (largest_idx != idx) {
    // Swap. Heapify needs to propagate downward.
    heap_swap(heap, idx, largest_idx);
    heap_heapify_down(heap, largest_idx);
  }
}

/*
Heapify upward. This happens when we replace a node by a larger number.
*/
static void heap_heapify_up(heap_t* heap, const size_t idx) {
  if (idx <= HEAP_IDX_ROOT) return;

  const size_t parent = HEAP_GET_PARENT_IDX(idx);
  if (heap_compare_GT(heap, idx, parent) != parent) {
    heap_swap(heap, idx, parent);
    heap_heapify_up(heap, parent);
  }
}

list_node_t heap_allocate(heap_t* heap, const size_t size) {
  list_node_t root = heap_get_root(heap);
  if (!root) {
    return NULL;
  }

  if (size > root->size) {
    // The required size is larger than the largest block size.
    return NULL;
  }

  // First, try to find the smallest node on the top of the heap >= size
  list_node_t victim = root;
  const size_t heap_size = heap->size;
  const size_t left_idx = HEAP_GET_LEFT_IDX(HEAP_IDX_ROOT);
  const size_t right_idx = HEAP_GET_RIGHT_IDX(HEAP_IDX_ROOT);
  list_node_t left = left_idx < heap_size ? 
    heap->node_array[left_idx] : NULL;
  list_node_t right = right_idx < heap_size ? 
    heap->node_array[right_idx] : NULL;
  
  if (left && left->size >= size) {
    victim = left;
  }
  if (right && right->size >= size && right->size < victim->size) {
    victim = right;
  }

  if (victim->size == size) {
    // Perfect match!
    heap_remove_idx(heap, victim->idx);
    victim->is_free = false;
    return victim;
  }

  // The block is larger than what's being asked. Need to split...
  list_node_t ret_node = list_create_node(
    victim->prev,
    victim,
    victim->addr,
    size,
    HEAP_IDX_NULL,
    false
  ); // Notice that this is the only place that we create new node.
  if (!ret_node) {
    // Running out of memory to allocate node...
    // But we can still return the allocated block
    heap_remove_idx(heap, victim->idx);
    victim->is_free = false;
    return victim;
  }

  victim->prev->next = ret_node;
  victim->prev = ret_node;
  victim->addr += size;
  victim->size -= size;
  heap_heapify_down(heap, victim->idx);
  return ret_node;
}

bool heap_free(heap_t* heap, list_node_t node) {
  // The idea is that we want to merge internal free blocks,
  // and try to avoid insert new node into the heap.
  assert(!node->is_free && "The current node must not be free already");
  if (node->next->is_free) {
    // The next block is free.
    if (node->prev->is_free) {
      // First check if the previous block is also free.
      // Merge and remove the previous is true.
      list_node_t previous = node->prev;
      node->size += previous->size;
      node->addr = previous->addr;
      assert(heap->node_array[previous->idx] == previous);
      heap_remove_idx(heap, previous->idx);
      node_unlink(previous); 
      free(previous);
    }
    
    // Merge the current into next.
    list_node_t next = node->next;
    next->size += node->size;
    next->addr = node->addr;
    node_unlink(node);
    free(node);
    // Restructure the heap as we increase the size of "next" node.
    heap_heapify_up(heap, next->idx);
    return true;

  } else if (node->prev->is_free) {
    // Merge the current block into the previous.
    list_node_t previous = node->prev;
    previous->size += node->size;
    node_unlink(node);
    free(node);

    // Restructure the heap as we increase the size of "next" node.
    heap_heapify_up(heap, previous->idx);
    return true;

  } else {
    // Nothing we can merge. Insert the current node into the heap.
    bool ret_val = heap_insert(heap, node);
    return ret_val;
  }
}

bool heap_init(heap_t* heap, uint8_t* addr, const size_t size) {
  list_node_t* array = 
    (list_node_t*) malloc(sizeof(list_node_t) * INIT_HEAP_CAPACITY);
  if (!array) {
    fprintf(stderr, "Error: failed to allocate memory for the heap.\n");
    return false;
  }

  array[0] = NULL;

  *heap = (heap_t) {
    .capacity = INIT_HEAP_CAPACITY,
    .size = 2,
    .node_array = array
  };

  if (!list_init(&(heap->node_list), addr, size)) {
    free(array);
    return false;
  }

  list_node_t head = heap->node_list.virtual_head->next;
  heap_set_idx(heap, 1, head);

  return true;
}