#ifndef SPBAU_BIOINF_2016_PAIR_READ_GRAPH_H
#define SPBAU_BIOINF_2016_PAIR_READ_GRAPH_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include <seqan/graph_types.h>

using namespace std;
using namespace seqan;

class PairReadGraph {
 private:
  static const int DEFAULT_MIN_CONTIG_LEN = 6000;
  static const int DEFAULT_MIN_COUNT = 1000;
  const double DEFAULT_DEF = 0.179;
  const int DEFAULT_MAX_CNT_EDGE = 3;


  typedef unsigned int Wight;

  typedef Graph<Directed<Wight>> DirG;
  typedef VertexDescriptor<DirG>::Type DirVert;

  vector <pair <int, int> > G[200];
  int count[200];

  map<DirVert, int> vertId;




  map<CharString, int> target_id;
  vector<CharString> target_name;
  String<CharString> vmp;
  String<CharString> emp;

  vector<DirVert> vertexById;
  DirG g;

  map<CharString, int> read1_pos;
  vector< unordered_map<int, int> > cnt;

  vector<int> max_edge;

  vector<double> target_coverage;
  vector<int> contig_len;

  BamFileIn fp;

  int pair_target(int x);

  void read_header_init();

  void process_one_first_read(BamAlignmentRecord read);

  void first_reads(char *file_name);

  pair<CharString, int> process_one_second_read(BamAlignmentRecord read);

  void inc_edge_weight(CharString read_name, int target_id);

  void add_edges(int min_count, CharString color, char* file_name);

  void second_reads(char *file_name, int min_count);

  CharString gen_random_color();

  CharString append_info(CharString property, char *lib_name, int x);

  void write_full_graph();

 public:
  int add_reads_to_graph(char *file_name1, char *file_name2, int min_count = DEFAULT_MIN_COUNT);

  void appendCoverageToMap();

  void write_graph();

  PairReadGraph() { }
};

#endif //SPBAU_BIOINF_2016_PAIR_READ_GRAPH_H
