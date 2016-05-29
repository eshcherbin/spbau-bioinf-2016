#include "pair_read_graph.h"

int PairReadGraph::pair_target(int x) {
    return x ^ 1;
}

void PairReadGraph::init_average_dist(int dist) {
    average_dist = dist;
}

void PairReadGraph::resizeVectorsOnInit(size_t length) {
    vertex_by_id.resize(length);
    edge_count.resize(length);
    target_coverage.resize(length);
    target_length.resize(length);
    count.resize(length);
    G.resize(length);
}

void PairReadGraph::appendVertexLabel(CharString name, int length) {
    String<char> label_text("label = \" name: ");
    append(label_text, name);
    append(label_text, "\n length = ");
    append(label_text, to_string(length));
    append(label_text, "\"");
    appendValue(vmp, label_text);
}

void PairReadGraph::addVertex(int i, CharString name, int length) {
    target_name.push_back(name);
    target_name.push_back(name);
    target_id[name] = static_cast<int>(target_name.size()) - 2;
    vertex_by_id[2 * i] = seqan::addVertex(graph);
    vertex_by_id[2 * i + 1] = seqan::addVertex(graph);
    vert_id[vertex_by_id[2 * i]] = 2 * i;
    vert_id[vertex_by_id[2 * i + 1]] = 2 * i + 1;
    target_length[2 * i] = length;
    target_length[2 * i + 1] = length;
}

void PairReadGraph::readHeaderInit() {
    typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

    BamHeader sam_header;
    readHeader(sam_header, bam_file);
    TBamContext const &bamContext = context(bam_file);
    size_t contig_num = length(contigNames(bamContext));
    resizeVectorsOnInit(2 * contig_num);
    CharString name;
    for (int i = 0; i < static_cast<int>(contig_num); ++i) {
        int length = contigLengths(bamContext)[i];
        if (length < min_contig_length) {
            target_name.push_back("");
            target_name.push_back("");
            continue;
        }
        name = contigNames(bamContext)[i];
        addVertex(i, name, length);
        appendVertexLabel(name, length);
        append(name, "-rev");
        appendVertexLabel(name, length);
    }
}


int PairReadGraph::readDist(BamAlignmentRecord read) {
    if (hasFlagRC(read) == false) {
        return (target_length[2 * read.rID] - read.beginPos);
    } else {
        return (read.beginPos + read.tLen);
    }
}

void PairReadGraph::processOneFirstRead(BamAlignmentRecord read) {
    readRecord(read, bam_file);
    CharString read_name = read.qName;
    if (length(read_name) > 1) {
        if (read_name[length(read_name) - 2] == '/' &&
                read_name[length(read_name) - 1] == '1') {
            resize(read_name, length(read_name) - 2);
        }
    }
    assert(read1_target.count(read_name) == 0);
    bool is_rev = hasFlagRC(read);
    int target = 2 * (read.rID);
    if (target < 0 || target_name[target] == "") {
        return;
    }
    if (is_rev) {
        target++;
    }
    read1_target[read_name] = target;
    read1_dist_to_end[read_name] = readDist(read);
    auto read_length = getAlignmentLengthInRef(read);
    auto contig_length = getContigLength(read, bam_file);
    target_coverage[target] += static_cast<double>(read_length) / contig_length;
    target_coverage[pair_target(target)] += static_cast<double>(read_length) /
        contig_length;
}

void PairReadGraph::firstReads(char *file_name, int dist) {
    open(bam_file, file_name);
    init_average_dist(dist);
    if (vertex_by_id.size() == 0) {
        readHeaderInit();
    } else {
        BamHeader sam_header;
        readHeader(sam_header, bam_file);
    }
    BamAlignmentRecord read;
    while (!atEnd(bam_file)) {
        processOneFirstRead(read);
    }
    close(bam_file);
}

pair<CharString, int>
PairReadGraph::processOneSecondRead(BamAlignmentRecord read) {
    readRecord(read, bam_file);
    CharString read_name = read.qName;
    if (length(read_name) > 1) {
        if (read_name[length(read_name) - 2] == '/' &&
                read_name[length(read_name) - 1] == '2') {
            resize(read_name, length(read_name) - 2);
        }
    }
    bool is_rev = hasFlagRC(read);
    int target = 2 * (read.rID);
    if (target < 0 || target_name[target] == "") {
        return make_pair("", -1);
    }
    if (!is_rev) {
        target++;
    }
    read2_dist_to_end[read_name] = readDist(read);
    auto read_length = getAlignmentLengthInRef(read);
    auto contig_length = getContigLength(read, bam_file);
    target_coverage[target] += static_cast<double>(read_length) / contig_length;
    target_coverage[pair_target(target)] += static_cast<double>(read_length) /
        contig_length;
    return make_pair(read_name, target);
}

CharString PairReadGraph::appendInfo(CharString property, char *lib_name, int x) {
    CharString lib = " label = \" library: " + string(lib_name) + "\n weight = "
        + to_string(x) + "\"";
    append(property, lib);
    return property;
}

void PairReadGraph::incEdgeWeight(CharString read_name, int target) {
    if (read1_target.count(read_name)) {
        if (read1_target[read_name] == target || 
                read1_target[read_name] == pair_target(target)) {
            return;
        }
        if (read1_dist_to_end[read_name] + read2_dist_to_end[read_name] > 
                average_dist) {
            return;
        }
        int verFID = read1_target[read_name], verSID = target, 
            verRFID = pair_target(verFID), verRSID = pair_target(verSID);
        edge_count[verFID][verSID]++;
        edge_count[verRSID][verRFID]++;
    }
}

void PairReadGraph::secondReads(char *file_name) {
    open(bam_file, file_name);
    BamHeader sam_header;
    readHeader(sam_header, bam_file);
    BamAlignmentRecord read;
    pair<CharString, int> read_info;
    while (!atEnd(bam_file)) {
        read_info = processOneSecondRead(read);
        if (read_info.second == -1) {
            continue;
        }
        incEdgeWeight(read_info.first, read_info.second);
    }
    close(bam_file);
}

int PairReadGraph::countEdgesBeforeBreak(int v, vector<pair<int, int> > edges) {
    int cnt = 1;
    int max_val = 0;
    if (edges.size() > 0) {
        max_val = edges[edges.size() - 1].first;
    }
    for (int i = static_cast<int>(edges.size()) - 2; i >= 0; --i) {
        int w0 = edges[i + 1].first, w1 = edges[i].first;
        if ((w0 - w1) < max_val * DEFAULT_DEF) {
            ++cnt;
        } else {
            break;
        }
    }
    if (cnt > DEFAULT_MAX_COUNT_EDGE) {
        cnt = 0;
    }
    count[v] = cnt;
    return cnt;
}

void PairReadGraph::addEdges(int min_count, char *file_name) {
    CharString color = genRandomColor();
    for (int v = 0; v < static_cast<int>(target_name.size()); ++v) {
        if (target_name[v] == "") {
            continue;
        }
        vector<pair<int, int> > edges;
        for (auto it = edge_count[v].begin(); it != edge_count[v].end(); ++it) {
            edges.push_back(make_pair(it->second, it->first));
        }
        if (edges.size() == 0) {
            continue;
        }
        sort(edges.begin(), edges.end());
        int cnt = countEdgesBeforeBreak(v, edges);
        for (int i = static_cast<int>(edges.size()) - 1; 
                i >= static_cast<int>(edges.size()) - cnt; --i) {
            if (edges[i].first < min_count) {
                break;
            }
            CharString property = appendInfo(color, file_name, edges[i].first);
            addEdge(graph, vertex_by_id[v], vertex_by_id[edges[i].second]);
            appendValue(emp, property);
        }
        for (int i = 0; i < static_cast<int>(edges.size()); ++i) {
            G[v].push_back(edges[i]);
        }
    }
}

void PairReadGraph::appendCoverageToMap() {
    int i = 0;
    for (int target = 0; 
            target < static_cast<int>(target_name.size()); target++) {
        if (target_coverage[target] == 0) {
            continue;
        }
        ostringstream coverage;
        coverage << setprecision(2) << fixed << target_coverage[target];
        eraseBack(vmp[i]);
        append(vmp[i], "\n coverage = " + coverage.str() + "\"");
        ++i;
    }
}

void PairReadGraph::writeFullGraph() {
    cerr << "start print full graph" << endl;
    ofstream out("full_graph.out");
    for (int i = 0; i < static_cast<int>(G.size()); ++i) {
        if (G[i].size() > 0) {
            out << i << " " << target_name[i] << ":\n";
            sort(G[i].rbegin(), G[i].rend());
            for (int j = 0; j < static_cast<int>(count[i]); ++j) {
                out << "    (" << G[i][j].first << ", " << 
                    G[i][j].second << ") \n";
            }
            out << "---" << endl;
            for (int j = count[i]; j < static_cast<int>(G[i].size()); ++j) {
                out << "    (" << G[i][j].first 
                    << ", " << G[i][j].second << ") \n";
            }
        }
    }
    out.close();
}

void PairReadGraph::writeGraph() {
    appendCoverageToMap();
    writeFullGraph();
    std::ofstream dotFile("graph.dot");
    writeRecords(dotFile, graph, vmp, emp, DotDrawing());
    dotFile.close();
}

void PairReadGraph::readAndFilterReads(char *file_name1, char *file_name2, 
        int dist) {
    cerr << "START" << endl;
    firstReads(file_name1, dist);
    cerr << "After first reads" << endl;
    secondReads(file_name2);
    cerr << "After second reads" << endl;
}

void PairReadGraph::histogram(const char *file_out_name) {
    ofstream out(file_out_name);
    int max_count = 0;
    for (int v = 0; v < static_cast<int>(edge_count.size()); ++v) {
        for (auto it = edge_count[v].begin(); it != edge_count[v].end(); ++it) {
            max_count = max(max_count, it->second);
        }
    }
    vector<int> histogram(static_cast<unsigned long>(max_count + 1));
    for (int v = 0; v < static_cast<int>(edge_count.size()); ++v) {
        for (auto it = edge_count[v].begin(); it != edge_count[v].end(); ++it) {
            histogram[it->second]++;
        }
    }
    for (int i = 0; i <= max_count; ++i) {
        if (histogram[i] != 0) {
            out << i << " " << histogram[i] << "\n";
        }
    }
    out.close();
}

int PairReadGraph::addReadsToGraph(char *file_name, int min_count) {
    addEdges(min_count, file_name);
    read1_target.clear();
    for (int i = 0; i < static_cast<int>(edge_count.size()); ++i) {
        edge_count[i].clear();
    }
    return 0;
}

CharString PairReadGraph::genRandomColor() {
    int color[3] = {rand() % 256, rand() % 256, rand() % 256};
    string res = "#";
    for (int i = 0; i < 3; ++i) {
        if (color[i] / 16 < 10) {
            res += (color[i] / 16) + '0';
        } else {
            res += (color[i] / 16) - 10 + 'a';
        }

        if (color[i] % 16 < 10) {
            res += (color[i] % 16) + '0';
        } else {
            res += (color[i] % 16) - 10 + 'a';
        }
    }
    res = "color = \"" + res + "\"";
    return CharString(res);
}


void PairReadGraph::setMinContigLength(int _min_contig_length) {
    this->min_contig_length = _min_contig_length;
}