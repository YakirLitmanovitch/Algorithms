#include <iostream>
#include "Graph.hpp"
#include "Algorithms.hpp"

int main() {
    using namespace graph;

    // יצירת גרף עם 5 קודקודים
    Graph g(10);

    // הוספת צלעות
    g.addEdge(0, 1, 5);
    g.addEdge(1, 2);
    g.addEdge(1, 3, 3);
    g.addEdge(2, 5, 2);
    g.addEdge(2, 4);
    g.addEdge(2, 6, 6);
    g.addEdge(3, 5, 2);
    g.addEdge(3, 6, 1);
    g.addEdge(5, 7, 4);
    g.addEdge(6, 9, 1);
    g.addEdge(9, 8, 6);
    


    std::cout << "The graph after the edges adding:" << std::endl;
    g.printGraph();

    // // מחיקת צלע קיימת
    // std::cout << "\n Delete the edge - (1-2)..." << std::endl;
    // g.removeEdge(1, 2);
    // g.printGraph();

    // // מחיקת צלע שאינה קיימת
    // std::cout << "\n Delete unexist edge" << std::endl;
    // g.removeEdge(3, 4);
    // g.printGraph();
    
    graph::Graph bfsTree = graph::Algorithms::bfs(g, 0);
    bfsTree.printGraph();

    graph::Graph dfsTree = graph::Algorithms::dfs(g, 0);
    dfsTree.printGraph();

    graph::Graph dijkstraTree = graph::Algorithms::dijkstra(g, 0);
    std::cout << "\nDijkstra's algorithm tree:" << std::endl;
    dijkstraTree.printGraph();

    graph::Graph primTree = graph::Algorithms::prim(g);
    std::cout << "\nPrim's algorithm tree:" << std::endl;
    primTree.printGraph();

    graph::Graph kruskalTree = graph::Algorithms::kruskal(g);
    std::cout << "\nKruskal's algorithm tree:" << std::endl;
    kruskalTree.printGraph();
    return 0;
}
