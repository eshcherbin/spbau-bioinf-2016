#include "htslib/sam.h"
#include "stdio.h"
#include <vector>
#include <algorithm>

using namespace std;

int main(int argc, char **argv) { 
    char *fn = argv[1];
    samFile *fp = sam_open(fn, "r");

    bam_hdr_t *sam_hdr = sam_hdr_read(fp);

    bam1_t* read = bam_init1();
    
    vector<int> reads_dist;

    while (sam_read1(fp, sam_hdr, read) == 0) {
       if ((read->core).mpos != 0 && (read->core).mpos >= (read->core).pos) {
           reads_dist.push_back((read->core).mpos - (read->core).pos);
       }
    }

    sort(reads_dist.begin(), reads_dist.end());

    for (int i = 0; i < (int)reads_dist.size(); ++i) {
        printf("%d ", reads_dist[i]);
    }

    sam_close(fp);
    return 0;
}

