#include "htslib/sam.h"
#include "stdio.h"
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include <cmath>
#include <map>

using namespace std;

map <string, int> target_id;
vector<string> target_name; 

map<string, int> read1_pos;

vector< vector<int> > G;

int pair_target(int x) {
    if (x & 1) {
        return x + 1;
    } else {
        return x - 1;
    }
}

int main(int argc, char **argv) {
    char *fn = argv[1];

    samFile *fp;
    fp = sam_open(fn, "r");

    bam_hdr_t *sam_hdr = sam_hdr_read(fp);

    int len = sam_hdr->n_targets;

    G.resize(2*len);

    for (int i = 0; i < len; ++i) {
        string name = string(sam_hdr->target_name[i]);

	target_name.push_back(name);
	target_name.push_back(name);

	target_id[name] = target_name.size() - 2;
    }

    bam1_t* read = bam_init1();

    while (sam_read1(fp, sam_hdr, read) == 0) {
        string read_name = string(bam_get_qname(read));
	bool is_rev = bam_is_rev(read);

	int target_id = 2*(read->core).tid;
	
	if (target_id < 0) {
	    continue;
	}

	if (is_rev) {
	    target_id++;
	}
	read1_pos[read_name] = target_id;
    }
    sam_close(fp);

    fn = argv[2];
    fp = sam_open(fn, "r");

    sam_hdr = sam_hdr_read(fp);

    read = bam_init1();

    while (sam_read1(fp, sam_hdr, read) == 0) {
        string read_name = string(bam_get_qname(read));
	bool is_rev = bam_is_rev(read);

	int target_id = 2*(read->core).tid;
	
	if (target_id < 0) {
	    continue;
	}

	if (!is_rev) {
	    target_id++;
	}
	if (read1_pos.count(read_name) == 0) {
	    continue;
	} else {
            G[read1_pos[read_name]].push_back(target_id);
	    G[pair_target(target_id)].push_back(pair_target(read1_pos[read_name]));
        }
    }
}
