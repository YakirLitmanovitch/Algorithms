// Mail : yakirli45@gmail.com
#include "MinHeap.hpp"
#include <climits>
#include <stdexcept>

namespace graph {

    MinHeap::MinHeap(int cap) : capacity(cap), size(0) {
        heapArray = new HeapNode[capacity]; // Array to store heap nodes
        positions = new int[capacity]; // Array to store positions of vertices in the heap
        for (int i = 0; i < capacity; i++) {
            positions[i] = -1;
        }
    }

    MinHeap::~MinHeap() {
        delete[] heapArray;
        delete[] positions;
    }

    bool MinHeap::isEmpty() const {
        return size == 0;
    }

    // Insert a new vertex with its distance into the heap
    void MinHeap::insert(int vertex, int distance) {
        if (size == capacity) {
            throw std::overflow_error("Heap is full");
        }
        heapArray[size] = {vertex, distance}; 
        positions[vertex] = size; 
        heapifyUp(size); // Maintain heap property
        size++;
    }

    // Maintain the heap property by moving the node at index up
    void MinHeap::heapifyUp(int index) {
        // Move the node at index up until the heap property is satisfied
        while (index > 0) {
            int parent = (index - 1) / 2;
            if (heapArray[index].distance < heapArray[parent].distance) {
                std::swap(heapArray[index], heapArray[parent]);
    
                // Update the positions of the vertices in the heap
                positions[heapArray[index].vertex] = index;
                positions[heapArray[parent].vertex] = parent;
    
                index = parent;
            } else {
                break;
            }
        }
    }

    // Maintain the heap property by moving the node at index down
    void MinHeap::heapifyDown(int index) {
        int smallest = index;
        int left = 2 * index + 1;
        int right = 2 * index + 2;

        if (left < size && heapArray[left].distance < heapArray[smallest].distance) {
            smallest = left;
        }
        if (right < size && heapArray[right].distance < heapArray[smallest].distance) {
            smallest = right;
        }
        if (smallest != index) {
            std::swap(heapArray[index], heapArray[smallest]);
            positions[heapArray[index].vertex] = index;
            positions[heapArray[smallest].vertex] = smallest;
            heapifyDown(smallest);
        }
    }

    // Extract the minimum node from the heap
    HeapNode MinHeap::extractMin() {
        if (isEmpty()) throw std::runtime_error("Heap is empty");

        HeapNode root = heapArray[0];
        heapArray[0] = heapArray[size - 1];
        positions[heapArray[0].vertex] = 0;
        size--;
        heapifyDown(0);
        return root;
    }

    // Decrease the key of a vertex in the heap
    void MinHeap::decreaseKey(int vertex, int newDistance) {
        int i = positions[vertex];
        heapArray[i].distance = newDistance;
        heapifyUp(i);
    }

    // Check if the heap contains a vertex
    bool MinHeap::contains(int vertex) const {
        return positions[vertex] != -1 && positions[vertex] < size;
    }

}
