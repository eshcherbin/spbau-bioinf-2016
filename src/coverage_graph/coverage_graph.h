#ifndef COVERAGE_GRAPH_H_
#define COVERAGE_GRAPH_H_

#include <bits/stdc++.h>
#include <seqan/bam_io.h>

using namespace std;
using namespace seqan;

// Creates a coverage graph of a target or several targets
// based on a alignment results stored in a SAM/BAM file
class CoverageGraph {
public:
    CoverageGraph(const char *filename, int step_size);

    void BuildGraph();
private:
    void ProcessRecord();
    void OutputStatistics();

    typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

    int step_size_;
    BamFileIn bam_file_;
    BamHeader bam_header_;
    TBamContext bam_context_;
    BamAlignmentRecord bam_record_;
    int num_targets;
    // how many reads cover a base in a particular target
    vector<vector<int> > reads_per_base; 
};

#endif
