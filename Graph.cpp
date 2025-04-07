// Mail : yakirli45@gmail.com
#include "Graph.hpp"
#include <iostream>
#include <stdexcept>
namespace graph{
    /**
     * @brief Constructs a graph with a specified number of vertices.
     * 
     * @param vertices The number of vertices in the graph.
     * 
     * @throws std::invalid_argument If the number of vertices is less than or equal to zero.
     */
    Graph::Graph(int vertices) : Vertices(vertices){
        neighbors = new Node*[Vertices];
        for (int i = 0; i < Vertices; i++) {
            neighbors[i] = nullptr;
        }
    }
    /**
     * @brief Destructor for the Graph class.
     * 
     * This destructor deallocates the memory used for the adjacency list of the graph.
     */
    Graph::~Graph(){
        for (int i = 0; i < Vertices; i++)
        {
            Node* current = neighbors[i];
            while (current){
                Node* temp = current;
                current = current->next;
                delete temp;
            }
        }
        delete[] neighbors;
    }

    /**
     * @brief Adds an edge between two vertices in the graph with a specified weight.
     * 
     * @param src The source vertex of the edge.
     * @param dest The destination vertex of the edge.
     * @param weight The weight of the edge. Must be non-negative.
     * 
     * @throws std::invalid_argument If the weight is negative.
     * @throws std::out_of_range If either `src` or `dest` is not a valid vertex index.
     */
    void Graph::addEdge(int src, int dest, int weight ){
        if (weight < 0){
            throw std::invalid_argument("weight cannot be negative");
        }
        if(src < 0 || src >= Vertices || dest < 0 || dest >= Vertices){
            throw std :: out_of_range("invalid vertex");
        }
        Node* newNode1 = new Node{dest, weight, neighbors[src]};
        neighbors[src] = newNode1;

        Node* newNode2 = new Node{src, weight, neighbors[dest]};
        neighbors[dest] = newNode2;
    }
    /**
     * @brief Adds a directed edge to the graph from a source vertex to a destination vertex with a specified weight.
     * 
     * @param src The source vertex of the edge. Must be within the range [0, Vertices-1].
     * @param dest The destination vertex of the edge. Must be within the range [0, Vertices-1].
     * @param weight The weight of the edge. Must be non-negative.
     * 
     * @throws std::invalid_argument If the weight is negative.
     * @throws std::out_of_range If the source or destination vertex is out of the valid range.
     */
    void Graph::addDirectedEdge(int src, int dest, int weight){
        if (weight < 0){
            throw std::invalid_argument("weight cannot be negative");
        }
        if(src < 0 || src >= Vertices || dest < 0 || dest >= Vertices){
            throw std :: out_of_range("invalid vertex");
        }
        Node* newNode = new Node{dest, weight, neighbors[src]};
        neighbors[src] = newNode;
    }
    /**
     * @brief Removes an edge between two vertices in the graph.
     * 
     * @param src The source vertex of the edge to be removed.
     * @param dest The destination vertex of the edge to be removed.
     * 
     * @throws std::out_of_range If either `src` or `dest` is out of the valid 
     *         range of vertex indices (0 to Vertices - 1).
     * 
     * @note If the edge does not exist, the function does nothing.
     */
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
    /**
     * @brief Prints the adjacency list representation of the graph.
     * @note This function assumes that the graph is represented as an adjacency
     * list, where each vertex has a linked list of its neighbors.
     */
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
    /**
     * @brief Retrieves the neighbors of a given vertex in the graph.
     * 
     * This function returns an array of integers representing the neighbors
     * of the specified vertex. The size of the array is returned via the 
     * reference parameter `size`. If the vertex has no neighbors, the function
     * returns a nullptr and sets `size` to 0.
     * 
     * @param vertex The vertex for which neighbors are to be retrieved.
     * @param size A reference to an integer where the size of the neighbors array will be stored.
     * @return A dynamically allocated array of integers containing the neighbors of the vertex,
     *         or nullptr if the vertex has no neighbors.
     * 
     * @throws std::out_of_range If the specified vertex is invalid (less than 0 or greater than or equal to the number of vertices).
     * 
     * @note The caller is responsible for deallocating the returned array using `delete[]` to avoid memory leaks.
     */
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
    /**
     * @brief Checks if an edge exists between two vertices in the graph.
     * 
     * @param src The source vertex ID.
     * @param dest The destination vertex ID.
     * @return true if an edge exists from src to dest, false otherwise.
     * 
     * @throws std::out_of_range if src or dest is not a valid vertex ID.
     */
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
    /**
     * @brief Retrieves the weight of the edge between two vertices in the graph.
     * 
     * @param src The source vertex ID.
     * @param dest The destination vertex ID.
     * @return The weight of the edge between the source and destination vertices.
     * 
     * @throws std::out_of_range If either the source or destination vertex ID is invalid.
     * @throws std::runtime_error If there is no edge between the source and destination vertices.
     */
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