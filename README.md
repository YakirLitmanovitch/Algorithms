# Graph Algorithms Project

## Overview
This project implements various graph algorithms using C++ to demonstrate concepts such as graph traversal, shortest path finding, and minimum spanning tree construction. The project is modular, with separate classes for graph representation, algorithms, and utility structures like heaps and union-find. It is designed for educational purposes and provides a clear and efficient implementation of fundamental graph algorithms.

## Features
The project includes the following graph algorithms:
1. **Breadth-First Search (BFS)**: Traverses the graph level by level starting from a given vertex.
2. **Depth-First Search (DFS)**: Explores as far as possible along each branch before backtracking.
3. **Dijkstra's Algorithm**: Finds the shortest path from a source vertex to all other vertices in a weighted graph.
4. **Prim's Algorithm**: Constructs a Minimum Spanning Tree (MST) using a greedy approach.
5. **Kruskal's Algorithm**: Constructs an MST by sorting edges and using a union-find data structure.

## Project Structure
The project is organized into the following files:

### 1. **`main.cpp`**
The entry point of the program. It demonstrates the usage of the implemented graph algorithms by:
- Creating a graph.
- Adding edges to the graph.
- Running BFS, DFS, Dijkstra's, Prim's, and Kruskal's algorithms.
- Printing the results of each algorithm.

### 2. **`Graph.hpp` and `Graph.cpp`**
These files define and implement the `Graph` class, which represents a graph using an adjacency list.

#### Key Features:
- **Graph Representation**: Uses an adjacency list for efficient storage and traversal.
- **Edge Management**: Functions to add and remove edges (directed and undirected).
- **Utility Functions**:
  - `getNeighbors`: Returns the neighbors of a given vertex.
  - `getEdgeWeight`: Retrieves the weight of an edge.
  - `printGraph`: Prints the graph's adjacency list.

### 3. **`Algorithms.hpp` and `Algorithms.cpp`**
These files define and implement the `Algorithms` class, which contains the graph algorithms.

#### Implemented Algorithms:
- **BFS**: Traverses the graph level by level and constructs a BFS tree.
- **DFS**: Recursively explores the graph and constructs a DFS tree.
- **Dijkstra's Algorithm**: Finds the shortest path from a source vertex to all other vertices.
- **Prim's Algorithm**: Constructs an MST using a priority queue (min-heap).
- **Kruskal's Algorithm**: Constructs an MST using edge sorting and a union-find data structure.

### 4. **`MinHeap.hpp` and `MinHeap.cpp`**
These files define and implement the `MinHeap` class, which is used in Dijkstra's and Prim's algorithms.

#### Key Features:
- **Heap Operations**:
  - `insert`: Inserts a vertex with a given priority.
  - `extractMin`: Removes and returns the vertex with the smallest priority.
  - `decreaseKey`: Updates the priority of a vertex.
- **Utility Functions**:
  - `isEmpty`: Checks if the heap is empty.
  - `contains`: Checks if a vertex is in the heap.

### 5. **`UnionFind.hpp` and `UnionFind.cpp`**
These files define and implement the `UnionFind` class, which is used in Kruskal's algorithm.

#### Key Features:
- **Union-Find Operations**:
  - `find`: Finds the root of a set containing a given element.
  - `unite`: Merges two sets.
  - `connected`: Checks if two elements are in the same set.

### 6. **Makefile**
The Makefile automates the build process and provides the following targets:
- `all`: Builds the entire project.
- `run`: Builds and runs the project.
- `valgrind`: Runs the project with `valgrind` to check for memory leaks.
- `clean`: Removes all build artifacts.

## How to Build and Run
1. **Clone the Repository**:
   ```bash
   git clone <repository-url>
   cd FirstProg
   ```

2. **Build the Project**:
   ```bash
   make
   ```

3. **Run the Project**:
   ```bash
   make Main
   ```

4. **Check for Memory Leaks**:
   ```bash
   make valgrind
   ```

5. **Clean Build Files**:
   ```bash
   make clean
   ```

## Example Output
Here is an example of the program's output:
...
### The graph after the edges adding:

![Adding edges](https://github.com/YakirLitmanovitch/Algorithms/blob/main/images/addingEdges.png)

...
### BFS starting from vertex 0
![Bfs process](https://github.com/YakirLitmanovitch/Algorithms/blob/main/images/bfsOp1.png)
![Bfs process](https://github.com/YakirLitmanovitch/Algorithms/blob/main/images/bfsOp2.png)
![Bfs process](https://github.com/YakirLitmanovitch/Algorithms/blob/main/images/bfsTree.png)

...

### DFS starting from vertex 0

![Bfs process](https://github.com/YakirLitmanovitch/Algorithms/blob/main/images/dfs.png)

...

### Dijkstra's algorithm starting from vertex 0:

![Dijkstra Pros](https://github.com/YakirLitmanovitch/Algorithms/blob/main/images/DijkstraPros.png)

![Dijkstra Pros](https://github.com/YakirLitmanovitch/Algorithms/blob/main/images/DijkstraTree.png)

...

### Prim's algorithm starting

![Dijkstra Pros](https://github.com/YakirLitmanovitch/Algorithms/blob/main/images/primMST.png)

...

### Kruskal's algorithm starting

![Dijkstra Pros](https://github.com/YakirLitmanovitch/Algorithms/blob/main/images/kruskalMST.png)

...


## Memory Management
The project uses dynamic memory allocation for graph representation and algorithm data structures. All dynamically allocated memory is properly freed, as verified by `valgrind`:
```
==451254== All heap blocks were freed -- no leaks are possible
==451254== ERROR SUMMARY: 0 errors from 0 contexts
```

