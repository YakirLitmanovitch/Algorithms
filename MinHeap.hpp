#ifndef MINHEAP_HPP
#define MINHEAP_HPP

namespace graph {

    struct HeapNode {
        int vertex;
        int distance;
    };

    class MinHeap {
        public:
        MinHeap(int capacity);
        ~MinHeap();

        bool isEmpty() const;
        void insert(int vertex, int distance);
        HeapNode extractMin();
        void decreaseKey(int vertex, int newDistance);
        bool contains(int vertex) const;

    private:
        int capacity;
        int size;
        HeapNode* heapArray;
        int* positions;

        void heapifyUp(int index);
        void heapifyDown(int index);
        

    
    };

}

#endif
