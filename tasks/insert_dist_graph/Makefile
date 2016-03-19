SHELL=bash

HTSLIB=-lhts

all: insert_dist_graph plot.png

insert_dist_graph: insert_dist.o
	g++ -L /usr/local/lib insert_dist.o -o insert_dist_graph $(HTSLIB)

insert_dist.o: insert_dist.cpp
	g++ -L /usr/local/lib -c insert_dist.cpp -o insert_dist.o

plot.dat: insert_dist_graph input
	./insert_dist_graph input 1 >plot.dat

plot.png: plot.dat
	gnuplot <<< "set term png; plot 'plot.dat' using 1:2 with lines" >plot.png

.PHONY: clean
.DELETE_ON_ERROR: plot.png plot.dat
clean:
	rm -f *.o insert_dist_graph *.dat *.png
