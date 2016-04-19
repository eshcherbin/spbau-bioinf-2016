#include <bits/stdc++.h>
#include "coverage_graph.h"

using namespace std;

CoverageGraph::CoverageGraph(const char *filename, int step_size)
        : step_size_(step_size), bam_file_(filename) {
    readHeader(bam_header_, bam_file_);
    bam_context_ = context(bam_file_);
    num_targets = length(contigNames(bam_context_));
    reads_per_base.resize(num_targets);
    for (int i = 0; i < num_targets; ++i)
        reads_per_base[i].resize(contigLengths(bam_context_)[i]);
}

void CoverageGraph::BuildGraph() {
    // process reads from SAM file while there are any
    while (!atEnd(bam_file_)) {
        readRecord(bam_record_, bam_file_);
        if (hasFlagUnmapped(bam_record_)) // skip if the aligner failed to align
            continue;
        ProcessRecord();
    }
    OutputStatistics();
}

// Outputs plot data for target "<target_name>" to file "plot_<target_name>.dat"
// A png plot can be built using plot.py
void CoverageGraph::ProcessRecord() {
    const static string cigar_operation_query = "MIS=X",
          cigar_operation_reference = "MDN=X";

    int pos = bam_record_.beginPos;
    for (auto &cigar_element : bam_record_.cigar) {
        // if operation consumes reference
        if (cigar_operation_reference.find(cigar_element.operation) != 
                string::npos) {
            //if operation consumes query
            if (cigar_operation_query.find(cigar_element.operation) != 
                    string::npos) {
                for (int j = 0; j < static_cast<int>(cigar_element.count); ++j)
                    ++reads_per_base[bam_record_.rID].at(pos++);
            }
            else
                pos += cigar_element.count;
        }
    }
}

void  CoverageGraph::OutputStatistics() {
    for (int target = 0; target < num_targets; ++target) {
        string filename = "plot_" + 
            string(toCString(contigNames(bam_context_)[target])) + 
            ".dat";
        ofstream data(filename.c_str());
        cout << "reference sequence " << contigNames(bam_context_)[target] <<
            ":\n";
        data << contigNames(bam_context_)[target] << '\n';

        int max_coverage = 0;
        double coverage = 0, average_coverage = 0;

        // compute statistics
        for (int i = 0; i < contigLengths(bam_context_)[target]; ++i) {
            if (reads_per_base[target][i] > max_coverage)
                max_coverage = reads_per_base[target][i];
            average_coverage += reads_per_base[target][i];
            if (reads_per_base[target][i])
                coverage += 1;
        }
        average_coverage /= contigLengths(bam_context_)[target];
        coverage /= contigLengths(bam_context_)[target];

        cout << "maximal coverage: " << max_coverage << " reads\n";
        cout << "average coverage: " << fixed << setprecision(2) <<
            average_coverage << " reads/position\n";
        cout << fixed << setprecision(2) << coverage * 100 << 
            "%% positions covered\n";

        data << max_coverage << ' ' << fixed << setprecision(2) << 
            average_coverage << ' ' << coverage << '\n';
        for (int l = 0; l < contigLengths(bam_context_)[target]; 
                l += step_size_) {
            int r = min(l + step_size_, contigLengths(bam_context_)[target]);
            double average_coverage_step = 0;
            for (int i = l; i < r; ++i)
                average_coverage_step += reads_per_base[target][i];
            average_coverage_step /= (r - l);
            data << l << ' ' << fixed << setprecision(5) << 
                average_coverage_step << '\n';
        }
    }
}
