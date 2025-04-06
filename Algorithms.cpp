#include "Algorithms.hpp"
#include "Graph.hpp"
#include "MinHeap.hpp"
#include "UnionFind.hpp"
#include <iostream>
#include <climits>


namespace graph
{
    // Edge structure for Kruskal's algorithm
    struct Edge {
        int src;
        int dest;
        int weight;
    };
    // Function to create a new edge
    Edge createEdge(int src, int dest, int weight) {
        Edge e;
        e.src = src;
        e.dest = dest;
        e.weight = weight;
        return e;
    }
    // Function to sort edges by weight
    void sortEdges(Edge* edges, int edgeCount) {
        // Using simple bubble sort 
        for (int i = 0; i < edgeCount - 1; i++) {
            for (int j = 0; j < edgeCount - i - 1; j++) {
                if (edges[j].weight > edges[j + 1].weight) {
                    Edge temp = edges[j];
                    edges[j] = edges[j + 1];
                    edges[j + 1] = temp;
                }
            }
        }
    }
    Graph Algorithms::bfs(const Graph& graph, int startVertex)
    {
        std::cout << "\nBFS starting from vertex " << startVertex << std::endl;
        if (startVertex < 0 || startVertex >= graph.getVertices()) {
            throw std::out_of_range("Invalid start vertex");
        }
        Graph bfsTree(graph.getVertices());

        bool* visited = new bool[graph.getVertices()]();
        int* queue = new int[graph.getVertices()];
        int front = 0, rear = 0;

        queue[rear++] = startVertex;
        visited[startVertex] = true;

        while (front < rear) {
            int currentVertex = queue[front++];
            std::cout << "Visited: " << currentVertex << std::endl;

            int size = 0;
            int* neighbors = graph.getNeighbors(currentVertex, size);
            if (neighbors){
                for (int i = 0; i < size; i++) {
                    int neighbor = neighbors[i];
                    std::cout << "Checking neighbor: " << neighbor << std::endl;
                    std::cout<<"In the visited list: " << visited[neighbor] << std::endl;
                    if (!visited[neighbor]) {
                        int weight = graph.getEdgeWeight(currentVertex, neighbor); 
                        bfsTree.addDirectedEdge(currentVertex, neighbor, weight); 
                        std::cout << "Adding edge to BFS tree: " << currentVertex << " -> " << neighbor << std::endl;
                        queue[rear++] = neighbor; 
                        visited[neighbor] = true;
                        }
                }
                delete[] neighbors;
            }
        }
        delete[] visited;
        delete[] queue;
        
        std::cout << "BFS completed." << std::endl;
        return bfsTree;
        
    }
    Graph Algorithms::dfs(const Graph& graph, int startVertex)
    {
        std::cout << "\nDFS starting from vertex " << startVertex << std::endl;
        if (startVertex < 0 || startVertex >= graph.getVertices()) {
            throw std::out_of_range("Invalid start vertex");
        }

        Graph dfsTree(graph.getVertices());
        bool* visited = new bool[graph.getVertices()];
        int* stack = new int[graph.getVertices()];
        int top = -1;

        for (int i = 0; i < graph.getVertices(); i++) {
            visited[i] = false;
        }

        stack[++top] = startVertex;
        visited[startVertex] = true;

        while (top >= 0) {
            int currentVertex = stack[top--];
            std::cout << "Visited: " << currentVertex << std::endl;

            int size;
            int* neighbors = graph.getNeighbors(currentVertex, size);
            for (int i = 0; i < size; i++) {
                int neighbor = neighbors[i];
                if (!visited[neighbor]) {
                    dfsTree.addDirectedEdge(currentVertex, neighbor);
                    stack[++top] = neighbor;
                    visited[neighbor] = true;
                }
            }
            delete[] neighbors;
        }

        delete[] visited;
        delete[] stack;

        std::cout << "DFS completed." << std::endl;
        return dfsTree;
    }
    Graph Algorithms::dijkstra(const Graph& graph, int startVertex) {
        std::cout << "\nDijkstra's algorithm starting from vertex " << startVertex <<":\n" << std::endl;

        if (startVertex < 0 || startVertex >= graph.getVertices()) {
            throw std::out_of_range("Invalid start vertex");
        }

        // Initialize distances and visited arrays
        int V = graph.getVertices();
        int* dist = new int[V]; // Shortest distances from startVertex
        bool* visited = new bool[V]; // Visited vertices
        int* parent = new int[V]; // Array to track the parent of each vertex
        Graph shortestPathTree(V); // Tree to store the shortest path
        MinHeap minHeap(V);

        // Initialize all arrays to infinity and visited to false
        for (int i = 0; i < V; i++) {
            dist[i] = INT_MAX;
            visited[i] = false;
            parent[i] = -1; 
        }

        dist[startVertex] = 0;
        minHeap.insert(startVertex, 0);

        
        while (!minHeap.isEmpty()) {
            HeapNode u = minHeap.extractMin(); // Get the vertex with the minimum distance
            std::cout << "Processing vertex: " << u.vertex << " with distance: " << dist[u.vertex] << std::endl;

            // If the vertex is already visited, skip it
            if (visited[u.vertex]) continue;
            visited[u.vertex] = true;

            int neighborCount;
            int* neighbors = graph.getNeighbors(u.vertex, neighborCount);

            // Relax all edges from u
            for (int i = 0; i < neighborCount; ++i) {
                int v = neighbors[i];
                if (visited[v]) continue;

                int weight = graph.getEdgeWeight(u.vertex, v);
                if (dist[u.vertex] + weight < dist[v]) {
                    // Remove the old edge from the shortest path tree
                    if (parent[v] != -1) {
                        shortestPathTree.removeEdge(parent[v], v);
                        std::cout << "Removing edge from shortest path tree: " << parent[v] << " -> " << v << std::endl;
                    }
                    // Update the distance and parent
                    dist[v] = dist[u.vertex] + weight;
                    parent[v] = u.vertex;

                    std::cout << "Updating distance of vertex " << v << " to " << dist[v] << std::endl;
                    minHeap.insert(v, dist[v]);
                    shortestPathTree.addDirectedEdge(u.vertex, v, weight);
                    std::cout << "Adding edge to shortest path tree: " << u.vertex << " -> " << v << std::endl;
                }
            }

            delete[] neighbors;
        }
        std::cout << "\nAll shortest paths results: " << std::endl;
        for (int i = 0; i < V; i++) {
            if (dist[i] == INT_MAX) {
                std::cout << "Vertex " << i << ": No path from " << startVertex << std::endl;
            } else {
                std::cout << "Vertex " << i << ": Distance = " << dist[i] << ", Parent = " << parent[i] << std::endl;
            }
        }
        delete[] dist;
        delete[] visited;
        delete[] parent;

        return shortestPathTree;
    }
    Graph Algorithms::prim(const Graph& graph) {
        std::cout << "\nPrim's algorithm starting\n" << std::endl;

        int V = graph.getVertices();
        if (V == 0) {
            throw std::invalid_argument("Graph is empty");
        }

        int startVertex = 0; // Default to starting from vertex 0
        bool* inMST = new bool[V];
        int* key = new int[V];
        int* parent = new int[V];
        Graph mst(V);
        MinHeap minHeap(V);

        for (int i = 0; i < V; i++) {
            key[i] = INT_MAX;
            inMST[i] = false;
            parent[i] = -1;
        }

        key[startVertex] = 0;
        minHeap.insert(startVertex, 0);

        while (!minHeap.isEmpty()) {
            HeapNode u = minHeap.extractMin();
            inMST[u.vertex] = true;

            int neighborCount;
            int* neighbors = graph.getNeighbors(u.vertex, neighborCount);

            if (neighbors && neighborCount > 0){
                for (int i = 0; i < neighborCount; ++i) {
                    int v = neighbors[i];
                    if (!inMST[v]) {
                        int weight = graph.getEdgeWeight(u.vertex, v);
                        if (weight < key[v]) {
                            key[v] = weight;
                            parent[v] = u.vertex;

                            if (minHeap.contains(v)) {
                                minHeap.decreaseKey(v, key[v]);
                            } else {
                                minHeap.insert(v, key[v]);
                            }
                        }
                    }
                }
                delete[] neighbors;
                neighbors = nullptr;
            }
            if (u.vertex != startVertex) { // Avoid adding edge for the starting vertex in the MST
                mst.addDirectedEdge(parent[u.vertex], u.vertex, key[u.vertex]);
                std::cout << "Adding edge to MST: " << parent[u.vertex] << " -> " << u.vertex << " with weight " << key[u.vertex] << std::endl;
            }
        }

        delete[] inMST;
        delete[] key;
        delete[] parent;

        return mst;
    }
    Graph Algorithms::kruskal(const Graph& graph) {
        std::cout << "\nKruskal's algorithm starting\n" << std::endl;
        
        int V = graph.getVertices();
        Graph mst(V);
        
        // Get all edges from the graph
        int maxPossibleEdges = V * (V - 1) / 2; // Max edges in a complete graph
        Edge* edges = new Edge[maxPossibleEdges];
        int edgeCount = 0;
        
        // Extract all edges from the graph
        for (int i = 0; i < V; i++) {
            int size;
            int* neighbors = graph.getNeighbors(i, size);
            
            if (neighbors && size > 0) {
                for (int j = 0; j < size; j++) {
                    int neighbor = neighbors[j];
                    
                    // To avoid duplicate edges in undirected graph, add only if i < neighbor
                    if (i < neighbor) {
                        int weight = graph.getEdgeWeight(i, neighbor);
                        edges[edgeCount++] = createEdge(i, neighbor, weight);
                        std::cout << "Found edge: " << i << " - " << neighbor << " with weight " << weight << std::endl;
                    }
                }
                delete[] neighbors;
            }
        }
        
        // Sort edges by weight
        sortEdges(edges, edgeCount);
        std::cout << "Sorted " << edgeCount << " edges by weight\n" << std::endl;
        
        // Create Union-Find data structure
        UnionFind unionFind(V);
        
        // Process edges in sorted order
        int edgesAdded = 0;
        for (int i = 0; i < edgeCount && edgesAdded < V - 1; i++) {
            int src = edges[i].src;
            int dest = edges[i].dest;
            int weight = edges[i].weight;
            
            // Check if adding this edge creates a cycle
            if (!unionFind.connected(src, dest)) {
                // Add edge to MST
                mst.addDirectedEdge(src, dest, weight);
                unionFind.unite(src, dest);
                edgesAdded++;
                std::cout << "Adding edge to MST: " << src << " - " << dest << " with weight " << weight << std::endl;
            } else {
                std::cout << "Skipping edge " << src << " - " << dest << " (would create cycle)" << std::endl;
            }
        }
        
        delete[] edges;
        
        std::cout << "Kruskal's algorithm completed. Added " << edgesAdded << " edges to MST." << std::endl;
        return mst;
    }
}