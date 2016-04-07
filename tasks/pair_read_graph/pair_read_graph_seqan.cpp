#include <bits/stdc++.h>
#include <seqan/bam_io.h>
#include <seqan/file.h>
#include <seqan/graph_types.h>
#include <seqan/graph_algorithms.h>

using namespace std;
using namespace seqan;

const int MIN_COUNT = 1000;

typedef Graph<Directed < > > DirG;
typedef VertexDescriptor<DirG>::Type DirVert;
typedef Iterator<DirG, EdgeIterator>::Type EdgeIt;

map <CharString, int> target_id;
vector<CharString> target_name; 

vector<DirVert> vertexById;
DirG G;

map<CharString, int> read1_pos;

map<pair<DirVert, DirVert>, int> cnt;
map<pair<DirVert, DirVert>, int> cntRm;

int pair_target(int x) {
    return x^1;
}

int main(int argc, char **argv) {
    BamFileIn fp(argv[1]);

    BamHeader sam_hdr;
    readHeader(sam_hdr, fp);

    typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

    TBamContext const & bamContext = context(fp);

    int len = length(contigNames(bamContext));

    vertexById.resize(2*len);

    for (int i = 0; i < len; ++i) {
        CharString name = contigNames(bamContext)[i];

	target_name.push_back(name);
	target_name.push_back(name);

	target_id[name] = target_name.size() - 2;

	vertexById[2*i] = addVertex(G);
	vertexById[2*i + 1] = addVertex(G);
    }

    BamAlignmentRecord read; 

    while(!atEnd(fp)) {
        readRecord(read, fp);

        CharString read_name = read.qName;
	bool is_rev = hasFlagRC(read);

	int target_id = 2*(read.rID);
	
	if (target_id < 0) {
	    continue;
	}

	if (is_rev) {
	    target_id++;
	}
	read1_pos[read_name] = target_id;
    }
    close(fp);

    open(fp, argv[2]);

    readHeader(sam_hdr, fp);

    while (!atEnd(fp)) {
        readRecord(read, fp);

        CharString read_name = read.qName;
	bool is_rev = hasFlagRC(read);

	int target_id = 2*(read.rID);
	
	if (target_id < 0) {
	    continue;
	}

	if (is_rev) {
	    target_id++;
	}
	if (read1_pos.count(read_name)) {
            if (read1_pos[read_name] == target_id) {
	        continue;
	    }
	    cnt[make_pair(vertexById[read1_pos[read_name]], vertexById[target_id])]++;
	    cnt[make_pair(vertexById[pair_target(target_id)], vertexById[pair_target(read1_pos[read_name])])]++;
	    
	    addEdge(G, vertexById[read1_pos[read_name]], vertexById[target_id]);
	    addEdge(G, vertexById[pair_target(target_id)], vertexById[pair_target(read1_pos[read_name])]);
        }
    }

    EdgeIt it(G);
    for (; !atEnd(it); goNext(it)) {
        DirVert u = sourceVertex(it), v = targetVertex(it);
        if (cnt[make_pair(u, v)] < MIN_COUNT) { 
	    removeEdge(G, getValue(it));
	} else if (cntRm[make_pair(u, v)] + 1 < cnt[make_pair(u, v)]) {
            removeEdge(G, getValue(it));
	    cntRm[make_pair(u, v)]++;	    
	}
    }
   
    std::ofstream dotFile("graph.dot");
    writeRecords(dotFile, G, DotDrawing());
    dotFile.close();
    return 0;
}
