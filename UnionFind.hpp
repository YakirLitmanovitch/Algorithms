// Mail : yakirli45@gmail.com
#ifndef UNION_FIND_HPP
#define UNION_FIND_HPP

namespace graph {
    /**
     * @brief A class representing a Union-Find data structure.
     * 
     * This class provides methods for union and find operations, which are useful
     * for Kruskal's algorithm and other applications involving disjoint sets.
     */
    class UnionFind {
    private:
        int* parent;
        int* rank;
        int size;

    public:
        UnionFind(int n);
        ~UnionFind();

        int find(int x);
        void unite(int x, int y);
        bool connected(int x, int y);
    };

}

#endif
