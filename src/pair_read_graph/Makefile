SHELL=bash

all: pair_read_graph

pair_read_graph: pair_read_graph.o
	g++ -std=c++11 pair_read_graph.o -o pair_read_graph

pair_read_graph.o: pair_read_graph.cpp
	g++ -std=c++11 -c pair_read_graph.cpp -o pair_read_graph.o

graph.png: graph.dot
	dot graph.dot -Tpng -o graph.png
