#include "pair_read_graph.h"

using namespace seqan;
using namespace std;

/**
 * [GRAPH/HISTOGRAM][name SAM file for left read], [name SAM file for right read],
 * [min edge wight for accept], [dist between pair read], [contig len]
 * if first key word "HISTOGRAM" only one library and [name SAM file for left read], [name SAM file for right read],
 * [dist between pair read], [contig len]
 */

int main(int argc, char **argv) {
  PairReadGraph prg;

  int pos = 2;

  if (argv[1][0] == 'G') {

    while (pos < argc) {
      cerr << pos << " " << argc << endl;
      prg.setMin_contig_len(atoi(argv[pos + 4]));
      prg.read_and_filter_reads(argv[pos], argv[pos + 1], atoi(argv[pos + 3]));
      prg.add_reads_to_graph(argv[pos], atoi(argv[pos + 2]));
      pos += 5;
    }

    prg.write_graph();
  } else {
    prg.setMin_contig_len(atoi(argv[pos + 3]));
    prg.read_and_filter_reads(argv[pos], argv[pos + 1], atoi(argv[pos + 2]));
    prg.histogram("histogram");
  }
  return 0;
}
