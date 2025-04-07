// Mail : yakirli45@gmail.com
#ifndef GRAPH_H
#define GRAPH_H
#include <iostream>

namespace graph {
    /**
     * @brief Class representing a graph using an adjacency list.
     * 
     * This class provides methods to add edges, remove edges, check for edge existence,
     * and retrieve neighbors of a vertex. It also allows printing the graph and checking
     * the number of vertices.
     */
    class Graph {
    private:
        int Vertices;
        struct Node
        {
            int vertexID, weight;
            Node* next;
        };
        Node** neighbors;

    public:
        Graph(int vertices);
        ~Graph();

        void addEdge(int src, int dest, int weight = 1);
        void addDirectedEdge(int src, int dest, int weight = 1);
        void removeEdge(int src, int dest);
        void printGraph() const;
        int getVertices() const { return Vertices; }
        int* getNeighbors(int vertex, int& size) const;
        bool edgeExists(int src, int dest) const;
        int getEdgeWeight(int src, int dest) const;

    };
}
#endif