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

  vector <vector <pair <int, int> > > G; //all edges in full graph
  vector<int> count; //count of edges from vertex i in result graph.

  int average_dist; //if distance between reads more that average_dist we don't take it in our graph.

  map<DirVert, int> vertId;

  map<CharString, int> target_id;
  vector<CharString> target_name;
  vector<double> target_coverage;
  vector<int> target_len;

  String<CharString> vmp;
  String<CharString> emp;

  vector<DirVert> vertexById;
  DirG g;

  map<CharString, int> read1_target;
  map<CharString, int> read1_dist_to_end;
  map<CharString, int> read2_dist_to_end;
  vector< unordered_map<int, int> > cnt; //cnount of edges between two target.

  BamFileIn fp;

  int pair_target(int x);

  void read_header_init();

  void process_one_first_read(BamAlignmentRecord read);

  void first_reads(char *file_name, int dist);

  pair<CharString, int> process_one_second_read(BamAlignmentRecord read);

  void inc_edge_weight(CharString read_name, int target_id);

  void add_edges(int min_count, char* file_name);

  void second_reads(char *file_name); //read and handling second reads.

  CharString gen_random_color();

  CharString append_info(CharString property, char *lib_name, int x);

  void write_full_graph();

  void resize_vectors_on_init(size_t len);

  void append_vertex_label(CharString name, int len);

  void add_vertex(int i, CharString name, int len);

  int cnt_edges_before_break(int v, vector<pair<int, int> > edges);

  void init_average_dist(int dist);

  int read_dist(BamAlignmentRecord read);

 public:
  int add_reads_to_graph(char *file_name, int min_count);

  void appendCoverageToMap();

  void write_graph();

  void read_and_filter_reads(char *file_name1, char *file_name2, int dist);

  void histogram(const char *file_out_name);
  PairReadGraph() { }
};

#endif //SPBAU_BIOINF_2016_PAIR_READ_GRAPH_H
