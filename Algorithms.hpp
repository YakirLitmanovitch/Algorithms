// Mail : yakirli45@gmail.com
#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP
#include "Graph.hpp"

namespace graph
{
    /**
     * @brief Class containing various graph algorithms.
     * 
     * This class provides static methods for performing common graph algorithms
     * such as BFS, DFS, Dijkstra's algorithm, Prim's algorithm, and Kruskal's algorithm.
     */
    class Algorithms
    {
        public:
            static Graph bfs(const Graph& graph, int startVertex);
            static Graph dfs(const Graph& graph, int startVertex);
            static Graph dijkstra(const Graph& graph, int startVertex);
            static Graph prim(const Graph& graph);
            static Graph kruskal(const Graph& graph);
    };
}
#endif