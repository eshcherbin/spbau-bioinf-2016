#ifndef INSERT_DIST_H_
#define INSERT_DIST_H_

#include "htslib/sam.h"
#include <vector>

using namespace std;

class InsertDist {
public:
    InsertDist(char *fn, int _delta);
    ~InsertDist();
    void printDistCount();
    int getMedDist();
    int getAverage();
    int getDeviation();
    void cutPartOfSegment();
private:
    int getLengthFromCIGAR(bam1_t* read);

    int delta;
    samFile *fp;
    vector <int> reads_dist;
    vector <int> dist_count;
};

#endif
