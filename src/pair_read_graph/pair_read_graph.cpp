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


    String<char> label_text("label = \" name: ");
    append(label_text, target_name[i]);
    append(label_text, " len = ");
    int len = contigLengths(bamContext)[i];
    append(label_text, to_string(len));
    append(label_text, "\"");

    appendValue(vmp, label_text);

    String<char> label_text2("label = \" name: ");
    append(label_text2, target_name[i]);
    append(label_text2, "-rev len = ");
    append(label_text2, to_string(len));
    append(label_text2, "\"");

    appendValue(vmp, label_text2);
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

  if (vertexById.size() == 0) {
    read_header_init();
  } else {
    BamHeader sam_hdr;
    readHeader(sam_hdr, fp);
  }

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

CharString PairReadGraph::append_info(CharString property, char *lib_name, int x) {
  CharString lib = " label = \" library: " + string(lib_name) + " wieght = " + to_string(x) + "\"";
  append(property, lib);
  return property;
}

void PairReadGraph::inc_edge_weight(CharString read_name, int target_id) {
  if (read1_pos.count(read_name)) {
    if (read1_pos[read_name] == target_id || read1_pos[read_name] == pair_target(target_id)) {
      return;
    }

    DirVert verF = vertexById[read1_pos[read_name]], verS = vertexById[target_id],
        verRF = vertexById[pair_target(read1_pos[read_name])], verRS = vertexById[pair_target(target_id)];

    cnt[make_pair(verF, verS)]++;
    cnt[make_pair(verRS, verRF)]++;
  }
}

void PairReadGraph::second_reads(char *file_name, int min_count) {
  open(fp, file_name);

  BamHeader sam_hdr;
  readHeader(sam_hdr, fp);

  BamAlignmentRecord read;

  CharString color = gen_random_color();

  pair<CharString, int> read_info;
  while (!atEnd(fp)) {
    read_info = process_one_second_read(read);
    if (read_info.second == -1) {
      continue;
    }

    inc_edge_weight(read_info.first, read_info.second);
  }

  add_edges(min_count, color, file_name);

  close(fp);
}

void PairReadGraph::add_edges(int min_count, CharString color, char* file_name) {
 for (map<pair<DirVert, DirVert>, int>::iterator it = cnt.begin(); it != cnt.end(); ++it) {
   DirVert verF = (*it).first.first, verS = (*it).first.second;
   DirVert verRF = vertexById[pair_target((verF))], verRS = vertexById[pair_target((verS))];

   CharString property = append_info(color, file_name, cnt[make_pair(verF, verS)]);
   if (cnt[make_pair(verF, verS)] >= min_count) {
     addEdge(g, verF, verS);
     appendValue(emp, property);
     addEdge(g, verRS, verRF);
     appendValue(emp, property);
   }
 }
}

void PairReadGraph::write_graph() {
  std::ofstream dotFile("graph.dot");
  writeRecords(dotFile, g, vmp, emp, DotDrawing());
  dotFile.close();
}

int PairReadGraph::add_reads_to_graph(char *file_name1, char *file_name2, int min_count) {
  cerr << "START" << endl;
  first_reads(file_name1);
  cerr << "After first reads" << endl;
  second_reads(file_name2, min_count);
  cerr << "After second reads" << endl;
  read1_pos.clear();
  cnt.clear();
  return 0;
}

CharString PairReadGraph::gen_random_color() {
  int color[3] = { rand() % 256, rand() % 256, rand() % 256 };
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
