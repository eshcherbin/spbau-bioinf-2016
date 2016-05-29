#ifndef SPBAU_BIOINF_2016_PAIR_READ_GRAPH_H
#define SPBAU_BIOINF_2016_PAIR_READ_GRAPH_H

#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include <seqan/graph_types.h>

using namespace std;
using namespace seqan;

class PairReadGraph {
public:
    PairReadGraph() = default;
    void setMinContigLength(int min_contig_length);
    int addReadsToGraph(char *file_name, int min_count);
    void appendCoverageToMap();
    void writeGraph();
    void readAndFilterReads(char *file_name1, char *file_name2, int dist);
    void histogram(const char *file_out_name);

private:
    typedef unsigned int Weight;
    typedef Graph<Directed<Weight> > DirG;
    typedef VertexDescriptor<DirG>::Type DirVert;

    static const int DEFAULT_MIN_CONTIG_LENGTH = 6000;
    static const int DEFAULT_MIN_COUNT = 1000;
    static constexpr double DEFAULT_DEF = 0.179;
    static const int DEFAULT_MAX_COUNT_EDGE = 3;

    int pair_target(int x);
    void readHeaderInit();
    void processOneFirstRead(BamAlignmentRecord read);
    void firstReads(char *file_name, int dist);
    pair<CharString, int> processOneSecondRead(BamAlignmentRecord read);
    void incEdgeWeight(CharString read_name, int target_id);
    void addEdges(int min_count, char* file_name);
    void secondReads(char *file_name); //read and handling second reads.
    CharString genRandomColor();
    CharString appendInfo(CharString property, char *lib_name, int x);
    void writeFullGraph();
    void resizeVectorsOnInit(size_t length);
    void appendVertexLabel(CharString name, int length);
    void addVertex(int i, CharString name, int length);
    int countEdgesBeforeBreak(int v, vector<pair<int, int> > edges);
    void init_average_dist(int dist);
    int readDist(BamAlignmentRecord read);

    vector<vector<pair<int, int> > > G; //all edges in full graph
    vector<int> count; //count of edges from vertex i in result graph.

    // if distance between reads more that average_dist
    // we don't take it in our graph.
    int average_dist; 

    int min_contig_length = DEFAULT_MIN_CONTIG_LENGTH;
    map<DirVert, int> vert_id;
    map<CharString, int> target_id;
    vector<CharString> target_name;
    vector<double> target_coverage;
    vector<int> target_length;
    String<CharString> vmp;
    String<CharString> emp;
    vector<DirVert> vertex_by_id;
    DirG graph;
    map<CharString, int> read1_target;
    map<CharString, int> read1_dist_to_end;
    map<CharString, int> read2_dist_to_end;
    //count of edges between two target.
    vector< unordered_map<int, int> > edge_count;
    BamFileIn bam_file;
};

#endif //SPBAU_BIOINF_2016_PAIR_READ_GRAPH_H
