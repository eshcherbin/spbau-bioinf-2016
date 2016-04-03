SHELL=bash

HTSLIB=-lhts

all: pair_read_graph

pair_read_graph: pair_read_graph.o
	g++ -L /usr/local/lib pair_read_graph.o -o pair_read_graph $(HTSLIB)

pair_read_graph.o: pair_read_graph.cpp
	g++ -L /usr/local/lib -c pair_read_graph.cpp -o pair_read_graph.o
