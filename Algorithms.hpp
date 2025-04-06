#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP
#include "Graph.hpp"

namespace graph
{
    
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