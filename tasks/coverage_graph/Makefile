HTSLIB=-lhts
CPP=g++

all: coverage_graph

coverage_graph: main.o
	$(CPP) main.o -o coverage_graph $(HTSLIB)

main.o: main.cpp
	$(CPP) -Wall -c main.cpp -o main.o

.PHONY: clean
clean:
	rm -f main.o coverage_graph plot.{dat,png}
