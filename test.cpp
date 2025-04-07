// Mail : yakirli45@gmail.com
#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"
#include "Algorithms.hpp"
#include "Graph.hpp"
#include "MinHeap.hpp"
#include "UnionFind.hpp"

using namespace graph;

TEST_CASE("Graph edge cases") {
    // Test graph with no edges
    Graph g1(3);
    CHECK(g1.edgeExists(0, 1) == false);
    CHECK(g1.edgeExists(1, 2) == false);

    // Test adding duplicate edges
    Graph g2(3);
    g2.addEdge(0, 1);
    g2.addEdge(0, 1); // Duplicate edge
    CHECK(g2.edgeExists(0, 1) == true);

    // Test self-loops
    Graph g3(3);
    g3.addEdge(0, 0);
    CHECK(g3.edgeExists(0, 0) == true);
}
TEST_CASE("Graph with weights - edge cases") {
    // Test graph with no weights
    Graph g1(3);
    CHECK_THROWS_AS(g1.getEdgeWeight(0, 1), std::runtime_error);

    // Test negative weights
    Graph g2(3);
    //g2.addEdge(0, 1, -2);
    CHECK_THROWS_AS(g2.addEdge(0, 1, -2), std::invalid_argument);

    // Test updating weights
    g2.addEdge(0, 1, 5.0); // Update weight
    CHECK(g2.getEdgeWeight(0, 1) == 5.0);
}
TEST_CASE("Delete edges") {
    Graph g(5);
    g.addEdge(0, 1);
    g.addEdge(1, 2);
    g.addEdge(2, 3);
    g.addEdge(3, 4);

    // Check if edges exist before deletion
    CHECK(g.edgeExists(0, 1) == true);
    CHECK(g.edgeExists(1, 2) == true);

    // Delete edges
    g.removeEdge(0, 1);
    g.removeEdge(1, 2);

    // Check if edges exist after deletion
    CHECK(g.edgeExists(0, 1) == false);
    CHECK(g.edgeExists(1, 2) == false);
}
TEST_CASE("Dijkstra's Algorithm") {
    Graph g(5);
    g.addEdge(0, 1, 2);
    g.addEdge(0, 2, 4);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 4);
    g.addEdge(3, 4, 5);

    auto dijkstraTree = Algorithms::dijkstra(g, 0);

    // Check if relevant edges exist in the Dijkstra tree
    CHECK(dijkstraTree.edgeExists(0, 1) == true);
    CHECK(dijkstraTree.edgeExists(1, 2) == true);
    CHECK(dijkstraTree.edgeExists(1, 3) == true);
    CHECK(dijkstraTree.edgeExists(3, 4) == true);

    // Ensure irrelevant edges are not in the Dijkstra tree
    CHECK(dijkstraTree.edgeExists(0, 2) == false);
}
TEST_CASE("Breadth-First Search (BFS)") {
    Graph g(6);
    g.addEdge(0, 1);
    g.addEdge(0, 2, 6);
    g.addEdge(1, 2, 2);
    g.addEdge(2, 3, 2);
    g.addEdge(0, 3, 7);

    auto bfsResult = Algorithms::bfs(g, 0);

    // Check BFS traversal order
    CHECK(bfsResult.edgeExists(0, 1) == true);
    CHECK(bfsResult.edgeExists(0, 2) == true); 
    CHECK(bfsResult.edgeExists(1, 2) == false); // 1-2 is not in the BFS tree
    CHECK(bfsResult.edgeExists(2, 3) == false); // 2-3 is not in the BFS tree
    CHECK(bfsResult.edgeExists(0, 3) == true); 
}
TEST_CASE("Depth-First Search (DFS)") {
    Graph g(6);
    g.addEdge(0, 1);
    g.addEdge(0, 2);
    g.addEdge(1, 2);
    g.addEdge(1, 3);
    g.addEdge(2, 3);
    g.addEdge(2, 4);
    g.addEdge(4, 5);

    auto dfsResult = Algorithms::dfs(g, 0);

    CHECK(dfsResult.edgeExists(0, 1) == true);
    CHECK(dfsResult.edgeExists(0, 2) == true); 
    CHECK(dfsResult.edgeExists(1, 2) == false); // 1-2 is not in the DFS tree
    CHECK(dfsResult.edgeExists(1, 3) == true); 
    CHECK(dfsResult.edgeExists(2, 3) == false); // 2-3 is not in the DFS tree
    CHECK(dfsResult.edgeExists(2, 4) == true);
    CHECK(dfsResult.edgeExists(4, 5) == true);
}
TEST_CASE("Prim's Algorithm") {
    Graph g(5);
    g.addEdge(0, 1, 2);
    g.addEdge(0, 2, 4);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 7);
    g.addEdge(2, 3, 3);
    g.addEdge(3, 4, 5);

    auto primTree = Algorithms::prim(g);

    // Check if relevant edges exist in the Prim tree
    CHECK(primTree.edgeExists(0, 1) == true);
    CHECK(primTree.edgeExists(1, 2) == true);
    CHECK(primTree.edgeExists(2, 3) == true);
    CHECK(primTree.edgeExists(3, 4) == true);

    // Ensure irrelevant edges are not in the Prim tree
    CHECK(primTree.edgeExists(0, 2) == false);
    CHECK(primTree.edgeExists(1, 3) == false);
}
TEST_CASE("Kruskal's Algorithm") {
    Graph g(5);
    g.addEdge(0, 1, 2);
    g.addEdge(0, 2, 4);
    g.addEdge(1, 2, 1);
    g.addEdge(1, 3, 7);
    g.addEdge(2, 3, 3);
    g.addEdge(3, 4, 5);

    auto kruskalTree = Algorithms::kruskal(g);

    // Check if relevant edges exist in the Kruskal tree
    CHECK(kruskalTree.edgeExists(0, 1) == true);
    CHECK(kruskalTree.edgeExists(1, 2) == true);
    CHECK(kruskalTree.edgeExists(2, 3) == true);
    CHECK(kruskalTree.edgeExists(3, 4) == true);

    // Ensure irrelevant edges are not in the Kruskal tree
    CHECK(kruskalTree.edgeExists(0, 2) == false);
    CHECK(kruskalTree.edgeExists(1, 3) == false);
}
TEST_CASE("Union-Find - unit and find") {
    UnionFind uf(5);
    uf.unite(0, 1);
    uf.unite(1, 2);
    uf.unite(3, 4);

    CHECK(uf.connected(0, 2) == true);
    CHECK(uf.connected(0, 3) == false);
    CHECK(uf.connected(3, 4) == true);
}
TEST_CASE("Union-find - find with path compression") {
    UnionFind uf(5);
    uf.unite(0, 1);
    uf.unite(1, 2);
    uf.unite(3, 4);

    CHECK(uf.find(2) == uf.find(0)); // Check if path compression works
    CHECK(uf.find(4) == uf.find(3)); // Check if path compression works
}
TEST_CASE("MinHeap") {
    MinHeap heap(10);
    heap.insert(0, 5);
    heap.insert(1, 3);
    heap.insert(2, 8);

    CHECK(heap.extractMin().vertex == 1);
    CHECK(heap.extractMin().vertex == 0);
    CHECK(heap.extractMin().vertex == 2);
}
TEST_CASE("MinHeap edge cases") {
    MinHeap heap(3);
    heap.insert(0, 5);
    heap.insert(1, 3);

    CHECK(heap.extractMin().vertex == 1);
    CHECK(heap.extractMin().vertex == 0);

    // Check if the heap is empty after extracting all elements
    CHECK(heap.isEmpty() == true);

    // Attempt to extract from an empty heap
    CHECK_THROWS_AS(heap.extractMin(), std::underflow_error);
}
TEST_CASE("MinHeap with invalid vertex") {
    MinHeap heap(3);
    heap.insert(0, 5);
    heap.insert(1, 3);

    // Attempt to decrease key of a non-existent vertex
    CHECK_THROWS_AS(heap.decreaseKey(2, 1), std::out_of_range);
}
TEST_CASE("MinHeap with invalid distance") {
    MinHeap heap(3);
    heap.insert(0, 5);
    heap.insert(1, 3);

    // Attempt to decrease key with an invalid distance
    CHECK_THROWS_AS(heap.decreaseKey(0, -1), std::invalid_argument);
}
TEST_CASE("Graph with no vertices") {
    Graph g(0);
    CHECK(g.getVertices() == 0);
    CHECK_THROWS_AS(g.addEdge(0, 1), std::out_of_range);
    CHECK_THROWS_AS(g.removeEdge(0, 1), std::out_of_range);
}
TEST_CASE("Graph - printing") {
    Graph g(3);
    g.addEdge(0, 1);
    g.addEdge(1, 2);

    std::ostringstream oss;
    std::streambuf* oldCoutBuffer = std::cout.rdbuf(oss.rdbuf());

    g.printGraph();

    std::cout.rdbuf(oldCoutBuffer); // Restore original cout buffer

    std::string output = oss.str();
    CHECK(output.find("vertex 0: -> (1, 1)") != std::string::npos);
    CHECK(output.find("vertex 1: -> (0, 1) -> (2, 1)") != std::string::npos);
}
TEST_CASE("Graph - get neighbors") {
    Graph g(3);
    g.addEdge(0, 1);
    g.addEdge(1, 2);

    int size = 0;
    int* neighbors = g.getNeighbors(1, size);

    CHECK(size == 2);
    CHECK(neighbors[0] == 0);
    CHECK(neighbors[1] == 2);

    delete[] neighbors;
}
