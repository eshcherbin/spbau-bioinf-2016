SHELL=bash
FLAGS=-O2 -g -DEBUG

all: main

main: main.o pair_read_graph.o
	g++ $(FLAGS) -std=c++11 main.o pair_read_graph.o -o  pair_read_graph

main.o: main.cpp pair_read_graph.h
	g++ $(FLAGS) -std=c++11 -Wall -c main.cpp -o main.o

pair_read_graph: pair_read_graph.o
	g++ $(FLAGS) -std=c++11 pair_read_graph.o -o

pair_read_graph.o: pair_read_graph.cpp pair_read_graph.h
	g++ $(FLAGS) -std=c++11 -c pair_read_graph.cpp -o pair_read_graph.o

graph.png: graph.dot
	dot graph.dot -Tpng -o graph.png

.PHONY: clean
clean:
	rm -rf *.o pair_read_graph graph.png
