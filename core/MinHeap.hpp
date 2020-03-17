#ifndef MINHEAP_HPP
#define MINHEAP_HPP

// A C++ program to demonstrate common Binary Heap Operations

#include<iostream>
#include<climits>

#include "Dijkstra.hpp"
#include <queue>

  
// A class for Min Heap
class MinHeap
{
    public:
    Dijkstra *harr; // pointer to an array of Dijkstra objects in heap
    int capacity; // maximum possible size of min heap
    int heap_size; // Current number of elements in min heap
public:
    // Constructor
    MinHeap(int capacity);
    ~MinHeap();
  
    // to heapify a subtree with the root at given index
    void MinHeapify(int );
    int parent(int i);
  
    // to get index of left child of node at index i
    int left(int i);
  
    // to get index of right child of node at index i
    int right(int i);
  
    // to extract the root which is the minimum element
    Dijkstra extractMin();
  
    // Decreases key value of key at index i to new_val
    void decreaseKey(int i, int new_val);
  
    // Returns the minimum key (key at root) from min heap
    Dijkstra getMin();
  
    // Deletes a key stored at index i
    void deleteKey(int i);
  
    // Inserts a new key 'k'
    void insertKey(Dijkstra k);
    
    int findElement(int v);
    
    void changePrevVertex(int i, int v);
    
    int getHeapSize();
    
    // Prototype of a utility function to swap two integers
    void swap(int x, int y);
};




#endif
