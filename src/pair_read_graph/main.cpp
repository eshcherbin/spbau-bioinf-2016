#include "pair_read_graph.h"

using namespace seqan;
using namespace std;

int main(int argc, char **argv) {
  PairReadGraph prg;

  int pos = 1;

  while (pos < argc) {
    prg.add_reads_to_graph(argv[pos], argv[pos + 1], atoi(argv[pos + 2]));
    pos += 3;
  }

  prg.write_graph();
  return 0;
}
