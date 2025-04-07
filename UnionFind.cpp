// Mail : yakirli45@gmail.com
#include "UnionFind.hpp"

namespace graph {
    /// Constructor for UnionFind
    UnionFind::UnionFind(int n) : size(n) {
        parent = new int[n];
        rank = new int[n];
        for (int i = 0; i < n; ++i) {
            parent[i] = i;
            rank[i] = 0;
        }
    }

    /// Destructor for UnionFind
    UnionFind::~UnionFind() {
        delete[] parent;
        delete[] rank;
    }

    /// Find the root of the set containing x
    int UnionFind::find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]); // path compression
        }
        return parent[x];
    }

    /// Union two sets
    void UnionFind::unite(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);
        if (rootX == rootY) return;

        // union by rank
        if (rank[rootX] < rank[rootY]) {
            parent[rootX] = rootY;
        } else if (rank[rootX] > rank[rootY]) {
            parent[rootY] = rootX;
        } else {
            parent[rootY] = rootX;
            rank[rootX]++;
        }
    }
    /// Check if two elements are in the same set
    bool UnionFind::connected(int x, int y) {
        return find(x) == find(y);
    }

}
