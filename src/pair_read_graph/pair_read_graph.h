#ifndef SPBAU_BIOINF_2016_PAIR_READ_GRAPH_H
#define SPBAU_BIOINF_2016_PAIR_READ_GRAPH_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include <seqan/graph_types.h>

using namespace std;
using namespace seqan;

class PairReadGraph {
 private:
  static const int DEFAULT_MIN_COUNT = 1000;

  typedef unsigned int Wight;

  typedef Graph<Directed<Wight>> DirG;
  typedef VertexDescriptor<DirG>::Type DirVert;

  map<CharString, int> target_id;
  vector<CharString> target_name;
  String<CharString> vmp;
  String<CharString> emp;

  vector<DirVert> vertexById;
  DirG g;

  map<CharString, int> read1_pos;
  map<pair<DirVert, DirVert>, int> cnt;

  BamFileIn fp;

  int pair_target(int x);

  void read_header_init();

  void process_one_first_read(BamAlignmentRecord read);

  void first_reads(char *file_name);

  pair<CharString, int> process_one_second_read(BamAlignmentRecord read);

  void add_edge_to_graph(CharString read_name, int target_id, int min_count);

  void second_reads(char *file_name, int min_count);

  void write_graph();

 public:
  int main(char *file_name1, char *file_name2, int min_count = DEFAULT_MIN_COUNT);

  PairReadGraph() { }
};

#endif //SPBAU_BIOINF_2016_PAIR_READ_GRAPH_H