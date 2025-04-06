#include "Graph.hpp"
#include <iostream>
#include <stdexcept>
namespace graph{
    Graph::Graph(int vertices) : Vertices(vertices){
        neighbors = new Node*[Vertices];
        for (int i = 0; i < Vertices; i++) {
            neighbors[i] = nullptr;
        }
    }
    Graph::~Graph(){
        for (int i = 0; i < Vertices; i++)
        {
            Node* current = neighbors[i];
            while (current){
                Node* temp = current;
                current = current->next;
                //std::cout << "Deleting node for vertex " << i << std::endl;
                delete temp;
            }
        }
        //std::cout << "Deleting neighbors array" << std::endl;
        delete[] neighbors;
    }

    void Graph::addEdge(int src, int dest, int weight ){
        if(src < 0 || src >= Vertices || dest < 0 || dest >= Vertices){
            throw std :: out_of_range("invalid vertex");
        }
        Node* newNode1 = new Node{dest, weight, neighbors[src]};
        neighbors[src] = newNode1;

        Node* newNode2 = new Node{src, weight, neighbors[dest]};
        neighbors[dest] = newNode2;
    }
    void Graph::addDirectedEdge(int src, int dest, int weight){
        if(src < 0 || src >= Vertices || dest < 0 || dest >= Vertices){
            throw std :: out_of_range("invalid vertex");
        }
        Node* newNode = new Node{dest, weight, neighbors[src]};
        neighbors[src] = newNode;
    }
    void Graph::removeEdge(int src, int dest){
        if(src < 0 || src >= Vertices || dest < 0 || dest >= Vertices){
            throw std :: out_of_range("invalid vertex");
        }
        if (!neighbors[src]) return;

        Node* cur = neighbors[src];
        Node* prev = nullptr;
        
        if (cur->vertexID == dest){
            neighbors[src] = cur->next;
            delete cur;
            removeEdge(dest, src);
            return;
        }
        
        while (cur && cur->vertexID != dest){
            prev = cur;
            cur = cur->next;
        }
        if (cur)
        {
            prev->next = cur->next;
            delete cur; 
            removeEdge(dest, src);
        }
        
        }
    void Graph::printGraph() const {
        for (int i = 0; i < Vertices; i++) {
            std::cout << "vertex " << i << ":";
            Node* current = neighbors[i];
            while (current) {
                std::cout << " -> (" << current->vertexID << ", " << current->weight << ")";
                current = current->next;
            }
            std::cout << std::endl;
        }
    }
    int* Graph::getNeighbors(int vertex, int& size) const {
        if (vertex < 0 || vertex >= Vertices) {
            throw std::out_of_range("Invalid vertex");
        }
        size = 0;
        Node* current = neighbors[vertex];
        while (current) {
            size++;
            current = current->next;
        }
        if (size == 0) {
            return nullptr; // Return nullptr when there are no neighbors
        }
        
        int* result = new int[size];
        //std::cout << "Allocating neighbors array for vertex " << vertex << " with size " << size << std::endl;
        current = neighbors[vertex];
        for (int i = 0; i < size; i++) {
            result[i] = current->vertexID;
            current = current->next;
        }
        return result;
    }
    bool Graph::edgeExists(int src, int dest) const {
        if (src < 0 || src >= Vertices || dest < 0 || dest >= Vertices) {
            throw std::out_of_range("Invalid vertex");
        }
        Node* current = neighbors[src];
        while (current) {
            if (current->vertexID == dest) {
                return true;
            }
            current = current->next;
        }
        return false;
    }
    int Graph::getEdgeWeight(int src, int dest) const {
        if (src < 0 || src >= Vertices || dest < 0 || dest >= Vertices) {
            throw std::out_of_range("Invalid vertex");
        }
        Node* current = neighbors[src];
        while (current) {
            if (current->vertexID == dest) {
                return current->weight;
            }
            current = current->next;
        }
        throw std::runtime_error("Edge does not exist");
    }
}