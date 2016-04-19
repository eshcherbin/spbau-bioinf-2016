#include "pair_read_graph.h"

using namespace std;
using namespace seqan;

int PairReadGraph::pair_target(int x) {
  return x ^ 1;
}

void PairReadGraph::read_header_init() {
  BamHeader sam_hdr;
  readHeader(sam_hdr, fp);

  typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;

  TBamContext const &bamContext = context(fp);

  size_t len = length(contigNames(bamContext));

  vertexById.resize(2 * len);

  CharString name;

  for (int i = 0; i < len; ++i) {
    name = contigNames(bamContext)[i];

    target_name.push_back(name);
    target_name.push_back(name);

    target_id[name] = (int) target_name.size() - 2;

    vertexById[2 * i] = addVertex(g);
    vertexById[2 * i + 1] = addVertex(g);
  }

  resizeVertexMap(vmp, g);

  for (int i = 0; i < 2*len; ++i) {
    String<char> label_text("label = ");
    append(label_text, target_name[i]);
    assignProperty(vmp, vertexById[i], label_text);
  }
}

void PairReadGraph::process_one_first_read(BamAlignmentRecord read) {
  readRecord(read, fp);

  CharString read_name = read.qName;
  if (length(read_name) > 1) {
    if (read_name[length(read_name) - 2] == '/' &&
        read_name[length(read_name) - 1] == '1') {
      resize(read_name, length(read_name) - 2);
    }
  }

  bool is_rev = hasFlagRC(read);

  int target_id = 2 * (read.rID);

  if (target_id < 0) {
    return;
  }

  if (is_rev) {
    target_id++;
  }
  read1_pos[read_name] = target_id;
}

void PairReadGraph::first_reads(char *file_name) {
  open(fp, file_name);

  read_header_init();

  BamAlignmentRecord read;

  while (!atEnd(fp)) {
    process_one_first_read(read);
  }
  close(fp);
}

pair<CharString, int> PairReadGraph::process_one_second_read(BamAlignmentRecord read) {
  readRecord(read, fp);

  CharString read_name = read.qName;
  if (length(read_name) > 1) {
    if (read_name[length(read_name) - 2] == '/' &&
        read_name[length(read_name) - 1] == '2') {
      resize(read_name, length(read_name) - 2);
    }
  }

  bool is_rev = hasFlagRC(read);

  int target_id = 2 * (read.rID);

  if (target_id < 0) {
    return make_pair("", -1);
  }

  if (!is_rev) {
    target_id++;
  }
  return make_pair(read_name, target_id);
}

void PairReadGraph::add_edge_to_graph(CharString read_name, int target_id, int min_count) {
  if (read1_pos.count(read_name)) {
    if (read1_pos[read_name] == target_id || read1_pos[read_name] == pair_target(target_id)) {
      return;
    }

    DirVert verF = vertexById[read1_pos[read_name]], verS = vertexById[target_id],
        verRF = vertexById[pair_target(read1_pos[read_name])], verRS = vertexById[pair_target(target_id)];

    cnt[make_pair(verF, verS)]++;
    cnt[make_pair(verRS, verRF)]++;

    if (cnt[make_pair(verF, verS)] == min_count) {
      addEdge(g, verF, verS);
      addEdge(g, verRS, verRF);
    }
  }
}

void PairReadGraph::second_reads(char *file_name, int min_count) {
  open(fp, file_name);

  BamHeader sam_hdr;
  readHeader(sam_hdr, fp);

  BamAlignmentRecord read;

  pair<CharString, int> read_info;
  while (!atEnd(fp)) {
    read_info = process_one_second_read(read);
    if (read_info.second == -1) {
      continue;
    }

    add_edge_to_graph(read_info.first, read_info.second, min_count);
  }
}

void PairReadGraph::write_graph() {
  std::ofstream dotFile("graph.dot");
  writeRecords(dotFile, g, vmp, DotDrawing());
  dotFile.close();
}

int PairReadGraph::main(char *file_name1, char *file_name2, int min_count) {
  cerr << "START" << endl;
  first_reads(file_name1);
  cerr << "After first reads" << endl;
  second_reads(file_name2, min_count);
  cerr << "After second reads" << endl;
  write_graph();
  cerr << "THE END" << endl;
  return 0;
}

int main(int argc, char **argv) {
  PairReadGraph prg;

  if (argc == 4) {
    prg.main(argv[1], argv[2], atoi(argv[3]));
  }
  else {
    prg.main(argv[1], argv[2]);
  }
  return 0;
}
