#include "htslib/sam.h"
#include "stdio.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>

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

    int getMedDist() {
        int start_pos = 0, end_pos = (int)dist_count.size() - 1;
        int cur_val1 = dist_count[start_pos], cur_val2 = dist_count[end_pos];
        while (start_pos != end_pos) {
	    if (cur_val1 > cur_val2) {
	        cur_val1 -= cur_val2;
		--end_pos;
		cur_val2 = dist_count[end_pos];
	    } else {
	        cur_val2 -= cur_val1;
		++start_pos;
		cur_val1 = dist_count[start_pos];
	    }
	}
	return start_pos * delta;	
    }

    int getAverage() {
        long long int sum = 0;
	int cnt = 0;
	for (int i = 0; i < (int)dist_count.size(); ++i) {
	    sum += dist_count[i] * (long long)i * delta;
	    cnt += dist_count[i];
	}
	return sum/cnt;
    }

    int getDeviation() {
        int aver = getAverage();
	long long sum = 0;
	int cnt = 0;
	for (int i = 0; i < (int)dist_count.size(); ++i) {
	    cnt += dist_count[i];
	    sum += (i*delta - aver)*(long long)(i*delta - aver)*dist_count[i];
	}
	return sqrt((double)sum/cnt);
    }

    void cutPartOfSegment() {
        int dev = getDeviation();
	int med = getMedDist();
	int sp = (med - 5*dev)/delta, ep = (med + 5*dev)/delta;
	for (int i = 0; i < sp; ++i) {
	    dist_count[i] = 0;
	}
	dist_count.resize(ep, 0);
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
    for (int i = 0; i < 5; ++i) {
        id.cutPartOfSegment();
    }

    id.printDistCount();

    cerr << "MedDist: " << id.getMedDist() << endl;
    cerr << "Average: " << id.getAverage()  << endl;
    cerr << "Deviation: " << id.getDeviation() << endl;
    return 0;
}

