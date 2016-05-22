#include "pair_read_graph.h"

using namespace seqan;
using namespace std;

/**
 * [name SAM file for left read], [name SAM file for right read], [min edge wight for accept], [dist between pair read]
 */

int main(int argc, char **argv) {
  PairReadGraph prg;

  int pos = 1;

  while (pos < argc) {
    prg.add_reads_to_graph(argv[pos], argv[pos + 1], atoi(argv[pos + 2]), atoi(argv[pos + 3]));
    pos += 4;
  }

  prg.write_graph();
  return 0;
}
