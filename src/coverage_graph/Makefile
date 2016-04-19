CPP=g++
FLAGS=-O2

all: main

main: main.o coverage_graph.o
	$(CPP) $(FLAGS) -std=c++11 main.o coverage_graph.o -o coverage_graph

coverage_graph.o: coverage_graph.cpp coverage_graph.h
	$(CPP) $(FLAGS) -std=c++11 -Wall -c coverage_graph.cpp -o coverage_graph.o

main.o: main.cpp coverage_graph.h
	$(CPP) $(FLAGS) -std=c++11 -Wall -c main.cpp -o main.o

.PHONY: clean
clean:
	rm -f *.o coverage_graph
