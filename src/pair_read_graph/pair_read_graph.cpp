#include "pair_read_graph.h"

using namespace std;
using namespace seqan;

int PairReadGraph::pair_target(int x) {
  return x ^ 1;
}

void PairReadGraph::init_average_dist(int dist) {
  average_dist = dist;
}

void PairReadGraph::resize_vectors_on_init(size_t len) {
  vertexById.resize(len);
  cnt.resize(len);
  target_coverage.resize(len);
  target_len.resize(len);

  count.resize(len);
  G.resize(len);
}

void PairReadGraph::append_vertex_label(CharString name, int len) {
  String<char> label_text("label = \" name: ");
  append(label_text, name);
  append(label_text, "\n len = ");
  append(label_text, to_string(len));
  append(label_text, "\"");

  appendValue(vmp, label_text);
}

void PairReadGraph::add_vertex(int i, CharString name, int len) {
  target_name.push_back(name);
  target_name.push_back(name);

  target_id[name] = (int) target_name.size() - 2;

  vertexById[2 * i] = addVertex(g);
  vertexById[2 * i + 1] = addVertex(g);

  vertId[vertexById[2 * i]] = 2 * i;
  vertId[vertexById[2 * i + 1]] = 2 * i + 1;

  target_len[2 * i] = len;
  target_len[2 * i + 1] = len;
}

void PairReadGraph::read_header_init() {
  BamHeader sam_hdr;
  readHeader(sam_hdr, fp);

  typedef FormattedFileContext<BamFileIn, void>::Type TBamContext;
  TBamContext const &bamContext = context(fp);

  size_t len = length(contigNames(bamContext));
  resize_vectors_on_init(2 * len);

  CharString name;
  for (int i = 0; i < len; ++i) {
    int length = contigLengths(bamContext)[i];

    if (length < DEFAULT_MIN_CONTIG_LEN) {
      target_name.push_back("");
      target_name.push_back("");
      continue;
    }

    name = contigNames(bamContext)[i];

    add_vertex(i, name, length);

    append_vertex_label(name, length);
    append(name,"-rev");
    append_vertex_label(name, length);
  }
}


int PairReadGraph::read_dist(BamAlignmentRecord read) {
  if (hasFlagRC(read) == false) {
    return (target_len[2*read.rID] - read.beginPos);
  } else {
    return (read.beginPos + read.tLen);
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


  assert(read1_target.count(read_name) == 0);
  assert(read.rNextId == -1);

  bool is_rev = hasFlagRC(read);

  int target_id = 2 * (read.rID);

  if (target_id < 0 || target_name[target_id] == "") {
    return;
  }

  if (is_rev) {
    target_id++;
  }
  read1_target[read_name] = target_id;
  read1_dist_to_end[read_name] = read_dist(read);

  auto read_length = getAlignmentLengthInRef(read);
  auto contig_length = getContigLength(read, fp);
  target_coverage[target_id] += static_cast<double>(read_length) / contig_length;
  target_coverage[pair_target(target_id)] += static_cast<double>(read_length) / contig_length;
}

void PairReadGraph::first_reads(char *file_name, int dist) {
  open(fp, file_name);

  init_average_dist(dist);
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

  if (target_id < 0 || target_name[target_id] == "") {
    return make_pair("", -1);
  }

  if (!is_rev) {
    target_id++;
  }

  read2_dist_to_end[read_name] = read_dist(read);

  auto read_length = getAlignmentLengthInRef(read);
  auto contig_length = getContigLength(read, fp);
  target_coverage[target_id] += static_cast<double>(read_length) / contig_length;
  target_coverage[pair_target(target_id)] += static_cast<double>(read_length) / contig_length;

  return make_pair(read_name, target_id);
}

CharString PairReadGraph::append_info(CharString property, char *lib_name, int x) {
  CharString lib = " label = \" library: " + string(lib_name) + "\n weight = " + to_string(x) + "\"";
  append(property, lib);
  return property;
}

void PairReadGraph::inc_edge_weight(CharString read_name, int target_id) {
  if (read1_target.count(read_name)) {
    if (read1_target[read_name] == target_id || read1_target[read_name] == pair_target(target_id)) {
      return;
    }

    if (read1_dist_to_end[read_name] + read2_dist_to_end[read_name] > average_dist) {
      return;
    }

    int verFID = read1_target[read_name], verSID = target_id, verRFID = pair_target(verFID), verRSID = pair_target(verSID);

    cnt[verFID][verSID]++;
    cnt[verRSID][verRFID]++;
  }
}

void PairReadGraph::second_reads(char *file_name, int min_count) {
  printf("%s\n", file_name);
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

int PairReadGraph::cnt_edges_before_break(int v, vector<pair<int, int> > edges) {
  int cnt = 1;

  int maxVal = 0;
  if (edges.size() > 0) {
    maxVal = edges[edges.size() - 1].first;
  }

  for (int i = (int)edges.size() - 2; i >= 0; --i) {
    int w0 = edges[i + 1].first, w1 = edges[i].first;
    if ((w0 - w1) < maxVal * DEFAULT_DEF) {
      ++cnt;
    } else {
      break;
    }
  }

  if (cnt > DEFAULT_MAX_CNT_EDGE) {
    cnt = 0;
  }

  count[v] = cnt;
  return cnt;
}

void PairReadGraph::add_edges(int min_count, CharString color, char* file_name) {
  for (int v = 0; v < target_name.size(); ++v) {
    if (target_name[v] == "") {
      continue;
    }

    vector< pair<int, int> > edges;

    for (auto it = cnt[v].begin(); it != cnt[v].end(); ++it) {
      edges.push_back(make_pair(it->second, it->first));
    }

    sort(edges.begin(), edges.end());
    int cnt = cnt_edges_before_break(v, edges);

    for (int i = (int)edges.size() - 1; i >= (int)edges.size() - cnt; --i) {

      if (edges[i].first < min_count) {
        break;
      }

      CharString property = append_info(color, file_name, edges[i].first);
      addEdge(g, vertexById[v], vertexById[edges[i].second]);
      appendValue(emp, property);
    }

    for (int i = 0; i < (int)edges.size(); ++i) {
      G[v].push_back(edges[i]);
    }
  }
}

void PairReadGraph::appendCoverageToMap()
{
  int i = 0;
  for (int target = 0; target < target_name.size(); target++)
  {
    if (target_coverage[target] == 0) {
      continue;
    }

    ostringstream coverage;
    coverage << setprecision(2) << fixed << target_coverage[target];
    eraseBack(vmp[i]);
    append(vmp[i], "\n coverage = " +
            coverage.str() + "\"");
    ++i;
  }
}

void PairReadGraph::write_full_graph() {
  cerr << "satrt print full graph" << endl;
  ofstream out("full_graph.out");

  for (int i = 0; i < 100; ++i) {
    if (G[i].size() > 0) {
      out << i << " " << target_name[i] << ":\n";
      sort(G[i].rbegin(), G[i].rend());
      for (int j = 0; j < (int)count[i]; ++j) {
        out << "    (" << G[i][j].first << ", " << G[i][j].second << ") \n";
      }
      out << "---" << endl;
      for (int j = count[i]; j < (int)G[i].size(); ++j) {
        out << "    (" << G[i][j].first << ", " << G[i][j].second << ") \n";
      }
    }
  }

  out.close();
}

void PairReadGraph::write_graph() {
  appendCoverageToMap();

  write_full_graph();

  std::ofstream dotFile("graph.dot");
  writeRecords(dotFile, g, vmp, emp, DotDrawing());
  dotFile.close();
}

int PairReadGraph::add_reads_to_graph(char *file_name1, char *file_name2, int min_count, int dist) {
  cerr << "START" << endl;
  first_reads(file_name1, dist);
  cerr << "After first reads" << endl;
  second_reads(file_name2, min_count);
  cerr << "After second reads" << endl;
  read1_target.clear();
  for (int i = 0; i < (int)cnt.size(); ++i) {
    cnt[i].clear();
  }
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
