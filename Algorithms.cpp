// Mail : yakirli45@gmail.com
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

    /**
     * @brief Performs a Breadth-First Search (BFS) on the given graph starting from a specified vertex.
     * 
     * This function traverses the graph using the BFS algorithm, starting from the given vertex.
     * It constructs a BFS tree, calculates distances from the start vertex to all other vertices,
     * and prints traversal details, including visited vertices, discovered edges, and distances.
     * 
     * @param graph The input graph to perform BFS on. It must provide methods to retrieve the number 
     *              of vertices, neighbors of a vertex, and edge weights.
     * @param startVertex The vertex from which BFS traversal begins. Must be within the valid range 
     *                    of vertices in the graph.
     * @return Graph A BFS tree represented as a graph, where edges correspond to parent-child 
     *               relationships discovered during the BFS traversal.
     * 
     * @throws std::out_of_range If the startVertex is invalid (less than 0 or greater than or equal 
     *                           to the number of vertices in the graph).
     */
    Graph Algorithms::bfs(const Graph& graph, int startVertex)
    {
        std::cout << "\nBFS starting from vertex " << startVertex << std::endl;

        if (startVertex < 0 || startVertex >= graph.getVertices()) {
            throw std::out_of_range("Invalid start vertex");
        }

        Graph bfsTree(graph.getVertices());
        
        // BFS initialization
        int numVertices = graph.getVertices();
        bool* visited = new bool[numVertices]();
        int* queue = new int[numVertices];
        int* parent = new int[numVertices];
        int* distance = new int[numVertices];

        for (int i = 0; i < numVertices; ++i) {
            parent[i] = -1;
            distance[i] = -1; 
        }


        int front = 0, rear = 0;
        queue[rear++] = startVertex; // Add start vertex to queue
        visited[startVertex] = true; // Mark it as visited
        distance[startVertex] = 0; 

        // Breadth-First Search (BFS) main loop:
        // Continues processing vertices in the queue until all reachable vertices are visited.
        while (front < rear) {
            int currentVertex = queue[front++];
            std::cout << "Visited: " << currentVertex << std::endl;

            int size = 0;
            int* neighbors = graph.getNeighbors(currentVertex, size);

            if (neighbors) {
                for (int i = 0; i < size; ++i) {
                    int neighbor = neighbors[i];

                    if (!visited[neighbor]) {
                        parent[neighbor] = currentVertex;
                        distance[neighbor] = distance[currentVertex] + 1;
                        visited[neighbor] = true;
                        queue[rear++] = neighbor;
                        std::cout << "Discovered: " << neighbor << " from " << currentVertex 
                                << ", distance: " << distance[neighbor] << std::endl;
                    }
                }
                delete[] neighbors;
            }
        }

        // After BFS is complete, build the BFS tree by connecting each vertex to its parent.
        // This creates a new graph representing the traversal path taken during BFS.
        for (int v = 0; v < numVertices; ++v) {
            if (parent[v] != -1) {
                int weight = graph.getEdgeWeight(parent[v], v);
                bfsTree.addDirectedEdge(parent[v], v, weight);
                std::cout << "Adding edge to BFS tree: " << parent[v] << " -> " << v << std::endl;
            }
        }

        // Print the distances from the start vertex to all other vertices.
        std::cout << "\nDistance from start vertex " << startVertex << ":\n";
        for (int v = 0; v < numVertices; ++v) {
            std::cout << "Vertex " << v << ": ";
            if (distance[v] != -1)
                std::cout << distance[v];
            else
                std::cout << "unreachable";
            std::cout << std::endl;
        }

        delete[] visited;
        delete[] queue;
        delete[] parent;
        delete[] distance;

        std::cout << "BFS completed." << std::endl;
        return bfsTree;
}
    /**
     * @brief Performs a Depth-First Search (DFS) traversal on the given graph starting from a specified vertex.
     * 
     * This function explores the graph using the DFS algorithm and constructs a DFS tree
     * representing the traversal. It uses an iterative approach with a stack to avoid recursion.
     * 
     * @param graph The input graph on which DFS is performed. It must provide methods to get the 
     *              number of vertices and the neighbors of a vertex.
     * @param startVertex The vertex from which the DFS traversal begins. Must be within the range 
     *                    [0, graph.getVertices() - 1].
     * @return Graph A new graph representing the DFS tree constructed during the traversal.
     * 
     * @throws std::out_of_range If the startVertex is invalid (less than 0 or greater than or equal 
     *                           to the number of vertices in the graph).
     * 
     * @note The function dynamically allocates memory for the visited array and stack, which are 
     *       properly deallocated before the function returns.
     * 
     * @details The DFS traversal visits each vertex exactly once and explores as far as possible 
     *          along each branch before backtracking. The traversal order is printed to the console.
     *          The DFS tree contains directed edges representing the traversal path.
     */
    Graph Algorithms::dfs(const Graph& graph, int startVertex)
    {
        std::cout << "\nDFS starting from vertex " << startVertex << std::endl;
        if (startVertex < 0 || startVertex >= graph.getVertices()) {
            throw std::out_of_range("Invalid start vertex");
        }
        // Initialize the DFS tree and visited array
        Graph dfsTree(graph.getVertices());
        bool* visited = new bool[graph.getVertices()];
        int* stack = new int[graph.getVertices()];
        int top = -1;

        for (int i = 0; i < graph.getVertices(); i++) {
            visited[i] = false;
        }

        stack[++top] = startVertex;
        visited[startVertex] = true;
        // DFS main loop:
        // Continues as long as there are vertices in the stack to explore
        while (top >= 0) {
            int currentVertex = stack[top--];
            std::cout << "Visited: " << currentVertex << std::endl;

            int size;
            int* neighbors = graph.getNeighbors(currentVertex, size);
            for (int i = 0; i < size; i++) {
                int neighbor = neighbors[i];
                if (!visited[neighbor]) {
                    dfsTree.addDirectedEdge(currentVertex, neighbor);
                    stack[++top] = neighbor; // Push the neighbor to the stack for further exploration
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
    /**
     * @brief Implements Dijkstra's algorithm to find the shortest paths from a starting vertex
     *        to all other vertices in a graph.
     * 
     * @param graph The input graph represented as an instance of the Graph class.
     * @param startVertex The starting vertex for the algorithm. Must be within the range [0, graph.getVertices()).
     * 
     * @return A Graph object representing the shortest path tree (SPT) from the starting vertex.
     * 
     * @throws std::out_of_range If the startVertex is invalid (less than 0 or greater than or equal to the number of vertices).
     * 
     * @details
     * - The function initializes distances to all vertices as infinity and marks all vertices as unvisited.
     * - It uses a MinHeap to efficiently extract the vertex with the smallest distance.
     * - For each vertex, it relaxes all its edges and updates the shortest path tree accordingly.
     * - The shortest path tree is represented as a directed graph with edges corresponding to the shortest paths.
     * - The function outputs intermediate steps, such as processing vertices, updating distances, and modifying the shortest path tree.
     * - At the end, it prints the shortest distances and parent vertices for all vertices.
     * 
     * @note The caller is responsible for ensuring that the input graph and startVertex are valid.
     */
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
    /**
     * @brief Implements Prim's algorithm to find the Minimum Spanning Tree (MST) of a graph.
     * 
     * This function takes a graph as input and returns a new graph representing the MST.
     * It uses a MinHeap to efficiently extract the minimum weight edge and updates the MST
     * iteratively by adding the smallest edge that connects a new vertex to the MST.
     * 
     * @param graph The input graph for which the MST is to be computed. It must be a connected graph.
     * @return Graph The resulting graph representing the MST.
     * 
     * @throws std::invalid_argument If the input graph is empty.
     * 
     * @details
     * - The algorithm initializes all vertices with infinite key values and sets the starting vertex's key to 0.
     * - A MinHeap is used to prioritize vertices based on their key values.
     * - For each vertex extracted from the heap, its neighbors are examined, and the key values are updated
     *   if a smaller edge weight is found.
     * - The parent array is used to track the edges included in the MST.
     * - The function dynamically allocates memory for arrays (inMST, key, parent) and ensures proper cleanup.
     * 
     * @note The graph is assumed to have non-negative edge weights.
     */
    Graph Algorithms::prim(const Graph& graph) {
        std::cout << "\nPrim's algorithm starting\n" << std::endl;

        int V = graph.getVertices();
        if (V == 0) {
            throw std::invalid_argument("Graph is empty");
        }
         // Initialize structures for Prim's algorithm:
        // key[i] holds the minimum weight to connect vertex i to the MST
        // parent[i] stores the parent of vertex i in the MST
        // inMST[i] indicates whether vertex i is already included in the MS
        int startVertex = 0; 
        bool* inMST = new bool[V];
        int* key = new int[V];
        int* parent = new int[V];
        Graph mst(V);
        MinHeap minHeap(V);

        // Initialize all keys as infinite and all vertices as not included in MST
        for (int i = 0; i < V; i++) {
            key[i] = INT_MAX;
            inMST[i] = false;
            parent[i] = -1;
        }

        key[startVertex] = 0;
        minHeap.insert(startVertex, 0);

        // Main loop of Prim's algorithm:
        // Repeatedly extract the minimum weight vertex and update its neighbors
        while (!minHeap.isEmpty()) {
            HeapNode u = minHeap.extractMin();
            inMST[u.vertex] = true;

            int neighborCount;
            int* neighbors = graph.getNeighbors(u.vertex, neighborCount);

            // Relax edges for all unvisited neighbors
            if (neighbors && neighborCount > 0){
                for (int i = 0; i < neighborCount; ++i) {
                    int v = neighbors[i];
                    if (!inMST[v]) {
                        int weight = graph.getEdgeWeight(u.vertex, v);
                        // If a better edge is found, update key and parent
                        if (weight < key[v]) {
                            key[v] = weight;
                            parent[v] = u.vertex;

                            // Update priority in heap
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
    /**
     * @brief Implements Kruskal's algorithm to find the Minimum Spanning Tree (MST) of a graph.
     * 
     * This function takes a graph as input and returns a new graph representing the MST.
     * It uses the Union-Find data structure to detect cycles and ensures that the MST
     * contains the minimum possible weight for all edges.
     * 
     * @param graph The input graph for which the MST is to be computed.
     *              It is assumed to be an undirected graph.
     * 
     * @return A Graph object representing the MST of the input graph.
     * 
     * @details
     * - The algorithm extracts all edges from the input graph and sorts them by weight.
     * - It processes the edges in sorted order, adding them to the MST if they do not
     *   form a cycle.
     * - The Union-Find data structure is used to efficiently check for cycles and
     *   merge connected components.
     * 
     * @note
     * - The input graph must provide methods to retrieve vertices, neighbors, and edge weights.
     * - The returned MST is represented as a directed graph, where each edge is added
     *   as a directed edge.
     */
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