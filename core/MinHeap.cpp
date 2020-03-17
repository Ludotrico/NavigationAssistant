#include "MinHeap.hpp"


// Constructor: Builds a heap from a given array a[] of given size
MinHeap::MinHeap(int cap)
{
    heap_size = 0;
    capacity = cap;
    harr = new Dijkstra[cap];
}
  
MinHeap::~MinHeap()
{
    delete [] harr;
}

int MinHeap::parent(int i)
{
    return (i-1)/2;
}
  
int MinHeap::left(int i)
{
    return (2*i + 1);
}
  
int MinHeap::right(int i)
{
    return (2*i + 2);
}
  
Dijkstra MinHeap::getMin()
{
    return harr[0];
}
      
int MinHeap::getHeapSize()
{
    return heap_size;
}


// Inserts a new key 'k'
void MinHeap::insertKey(Dijkstra d)
{
    if (heap_size == capacity)
        return;
  
    // First insert the new key at the end
    heap_size++;
    int i = heap_size - 1;
    harr[i] = d;
  
    // Fix the min heap property if it is violated
    while (i != 0 && harr[parent(i)] > harr[i])
    {
       swap(i, parent(i));
       i = parent(i);
    }
}
  
// Decreases value of key at index 'i' to new_val.  It is assumed that
// new_val is smaller than harr[i].
void MinHeap::decreaseKey(int i, int new_val)
{
    harr[i].d = new_val;
    while (i != 0 && harr[parent(i)] > harr[i])
    {
       swap(i, parent(i));
       i = parent(i);
    }
}
  
// Method to remove minimum element (or root) from min heap
Dijkstra MinHeap::extractMin()
{
    if (heap_size <= 0)
        return Dijkstra();

    if (heap_size == 1)
    {
        heap_size--;
        return harr[0];
    }
  
    // Store the minimum value, and remove it from heap
    Dijkstra root = harr[0];
    harr[0] = harr[heap_size-1];
    heap_size--;
    MinHeapify(0);
  
    return root;
}
  
  
// This function deletes key at index i. It first reduced value to minus
// infinite, then calls extractMin()
void MinHeap::deleteKey(int i)
{
    decreaseKey(i, std::numeric_limits<int>::min());
    extractMin();
}
  
// A recursive method to heapify a subtree with the root at given index
// This method assumes that the subtrees are already heapified
void MinHeap::MinHeapify(int i)
{
    int l = left(i);
    int r = right(i);
    int smallest = i;
    if (l < heap_size && harr[l] < harr[i])
        smallest = l;
    if (r < heap_size && harr[r] < harr[smallest])
        smallest = r;
    if (smallest != i)
    {
        swap(i, smallest);
        MinHeapify(smallest);
    }
}
  
// A utility function to swap two elements
void MinHeap::swap(int x, int y)
{
    Dijkstra temp = harr[x];
    harr[x] = harr[y];
    harr[y] = temp;
}

void MinHeap::changePrevVertex(int i, int v) {
    harr[i].pV = v;
}

int MinHeap::findElement(int v) {
    
    for(int i = 0; i < heap_size; i++)
        if(harr[i].v == v)
            return i;
  
    // We reach here when element is not
    // present in array
    return -1;
}
