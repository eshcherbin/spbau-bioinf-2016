#include "htslib/sam.h"
#include "stdio.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>

using namespace std;

class InsertDist {
private:
    int delta;
    samFile *fp;
    vector <int> reads_dist;
    vector <int> dist_count;
public:
    InsertDist(char *fn, int _delta = 1) {
        fp = sam_open(fn, "r");

        delta = _delta;
        bam_hdr_t *sam_hdr = sam_hdr_read(fp);

        bam1_t* read = bam_init1();
   
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
        dist_count.resize(max_pos/delta + 1, 0);

        for (int i = 0; i < (int)reads_dist.size(); ++i) {
            ++dist_count[reads_dist[i]/delta];
        }
    }

    void printDistCount() {
        for (int i = 0; i < (int)dist_count.size(); ++i) {
            printf("%d %d\n", i*delta, dist_count[i]);
        }
    }

    ~InsertDist() {
       sam_close(fp);
    }
};

int main(int argc, char **argv) { 
    char *fn = argv[1];

    int delta = 1;
    if (argc >= 3) {
        delta = atoi(argv[2]);
    }

    InsertDist id = InsertDist(fn, delta);
    id.printDistCount();

    return 0;
}

