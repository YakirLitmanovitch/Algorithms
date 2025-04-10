# Mail : yakirli45@gmail.com
CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -g
OBJS =main.o Graph.o Algorithms.o MinHeap.o UnionFind.o
TARGET = graph_program
TEST_TARGET = test_program
TEST_SRC = test.cpp

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

main.o: main.cpp Graph.hpp Algorithms.hpp
	$(CXX) $(CXXFLAGS) -c main.cpp

Graph.o: Graph.cpp Graph.hpp
	$(CXX) $(CXXFLAGS) -c Graph.cpp

Algorithms.o: Algorithms.cpp Algorithms.hpp Graph.hpp MinHeap.hpp UnionFind.hpp
	$(CXX) $(CXXFLAGS) -c Algorithms.cpp

MinHeap.o: MinHeap.cpp MinHeap.hpp
	$(CXX) $(CXXFLAGS) -c MinHeap.cpp

UnionFind.o: UnionFind.cpp UnionFind.hpp
	$(CXX) $(CXXFLAGS) -c UnionFind.cpp

Main: all
	./$(TARGET)


valgrind: all
	valgrind --leak-check=full --track-origins=yes ./$(TARGET)

test: $(TEST_SRC) Graph.o Algorithms.o MinHeap.o UnionFind.o
	$(CXX) $(CXXFLAGS) -o $(TEST_TARGET) $(TEST_SRC) Graph.o Algorithms.o MinHeap.o UnionFind.o
	./$(TEST_TARGET)

clean:
	rm -f $(TARGET) $(TEST_TARGET) $(OBJS)
