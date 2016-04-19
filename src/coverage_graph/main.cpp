#include <bits/stdc++.h>
#include "coverage_graph.h"

using namespace std;

const int kDefaultStepSize = 1000;

void print_help() {
    printf("usage: coverage_graph alignment [step]\n");
    exit(0);
}

int main(int argc, char **argv) {
    // process arguments
    if (argc == 1 || !strcmp(argv[1], "-h"))
        print_help();

    int step_size;
    if (argc == 3)
        step_size = atoi(argv[2]);
    else
        step_size = kDefaultStepSize;

    // build graph
    CoverageGraph graph(argv[1], step_size);
    graph.BuildGraph();

    return 0;
}
