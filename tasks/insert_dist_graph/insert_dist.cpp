#include "htslib/sam.h"
#include "stdio.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

using namespace std;
int main(int argc, char **argv) { 
    char *fn = argv[1];
    samFile *fp = sam_open(fn, "r");

    int delta = 1;
    if (argc >= 3) {
        delta = atoi(argv[2]);
    }

    bam_hdr_t *sam_hdr = sam_hdr_read(fp);

    bam1_t* read = bam_init1();
   
    vector <int> reads_dist;
    while (sam_read1(fp, sam_hdr, read) == 0) {
        if ((read->core.flag)&BAM_FMUNMAP || (read->core.flag)&BAM_FUNMAP || (read->core).mpos == 0) {
	    continue;
	}

        if ((read->core).pos >= (read->core).mpos) {
	    reads_dist.push_back((read->core).pos + (read->core).l_qseq - (read->core).mpos);
	}
    }

    sort(reads_dist.begin(), reads_dist.end());

    int max_pos = reads_dist[reads_dist.size() - 1];
    int dist_count[max_pos/delta + 1];

    for (int i = 0; i < max_pos/delta + 1; ++i) {
        dist_count[i] = 0;
    }

    for (int i = 0; i < (int)reads_dist.size(); ++i) {
        ++dist_count[reads_dist[i]/delta];
    }

    for (int i = 0; i <= max_pos/delta; ++i) {
        printf("%d %d\n", i*delta, dist_count[i]);
    }

    sam_close(fp);
    return 0;
}

