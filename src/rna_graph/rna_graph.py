#!/usr/bin/env python
import random
import argparse
from collections import Counter
import networkx
from alignment import AlignmentRecord


DEFAULT_DOT_FILENAME = 'graph.dot'
DEFAULT_MIN_IDY = 0
DEFAULT_MIN_ALIGNMENT_LEN = 0
DEFAULT_MIN_REFERENCE_LEN = 0

COLOR_RANDOM_SEED = 117  # to make the same colors for different runs


def get_random_color():
    return '{:.2} {:.2} {:.2}'.format(random.random(), random.random(),
                                      random.random())


def toggle_rev(string):
    return string[:-4] if string.endswith('-rev') else string + '-rev'


class RnaGraphBuilder:
    def __init__(self):
        self.ref_lengths = None
        self.graph = None
        self.weights = None

    def build(self, coords_file, min_idy, min_alignment_len,
              min_reference_len):
        self.ref_lengths = dict()
        self.graph = networkx.DiGraph()
        self.weights = Counter()
        with open(coords_file) as coords_file:
            one_query_records = []
            for line in coords_file.readlines():
                record = AlignmentRecord.from_line(line)
                # filtering
                if record.idy < min_idy or\
                        min(record.len1, record.len2) < min_alignment_len or\
                        record.lenr < min_reference_len:
                    continue
                if one_query_records and\
                        record.query != one_query_records[-1].query:
                    self._process_one_query(one_query_records)
                    one_query_records = []
                one_query_records.append(record)
            if one_query_records:
                self._process_one_query(one_query_records)

        for ref, length in self.ref_lengths.items():
            self.graph.add_node(ref,
                                label='{}\nlength = {}'.format(ref, length))
        random.seed(COLOR_RANDOM_SEED)
        color = get_random_color()
        for (first_ref, second_ref), weight in self.weights.items():
            self.graph.add_edge(first_ref, second_ref, color=color)

        return self

    def save(self, filename):
        vertices = list(sorted(networkx.weakly_connected_components(self.graph), key=len,
                     reverse=True))[0]
        output_graph = networkx.DiGraph()
        output_graph.add_nodes_from([(node, data) for node, data in self.graph.nodes(data=True) if node in vertices])
        output_graph.add_edges_from([(src, dest, data) for src, dest, data in self.graph.edges(data=True) if src in vertices and dest in vertices])
        networkx.drawing.nx_pydot.write_dot(output_graph, filename)
        return self

    def analyze(self):
        print(len(max(networkx.strongly_connected_components(self.graph))))
        print(' '.join(str(len(wcc)) for wcc in
                       sorted(networkx.weakly_connected_components(self.graph),
                              key=len)))
        print('\n'.join(sorted(list(sorted(networkx.weakly_connected_components(self.graph), key=len,
                     reverse=True))[0])))
        return self

    def _add_reference(self, record):
        self.ref_lengths[record.ref] = record.lenr
        self.ref_lengths[record.ref + '-rev'] = record.lenr

    def _process_one_query(self, records):
        assert records
        for record in records:
            if record.ref not in self.ref_lengths:
                self._add_reference(record)
        for i in range(len(records) - 1):
            j = i + 1
            first_ref = records[i].ref +\
                ('-rev' if records[i].s2 > records[i].e2 else '')
            second_ref = records[j].ref +\
                ('-rev' if records[j].s2 > records[j].e2 else '')
            if first_ref == second_ref or first_ref == toggle_rev(second_ref):
                continue
            self.weights[(first_ref, second_ref)] += 1
            self.weights[(toggle_rev(second_ref), toggle_rev(first_ref))] += 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('coords_file')
    parser.add_argument('-o', '--output', dest='dot_filename',
                        default=DEFAULT_DOT_FILENAME,
                        metavar='output_dot_file')
    parser.add_argument('-i', '--idy', dest='min_idy', type=float,
                        default=DEFAULT_MIN_IDY, metavar='min_idy')
    parser.add_argument('-l', '--minlength', dest='min_alignment_len',
                        type=int, default=DEFAULT_MIN_ALIGNMENT_LEN,
                        metavar='min_alignment_len')
    parser.add_argument('-r', '--refminlength', dest='min_reference_len',
                        type=int, default=DEFAULT_MIN_REFERENCE_LEN,
                        metavar='min_reference_len')
    args = parser.parse_args()
    RnaGraphBuilder().build(args.coords_file, args.min_idy,
                            args.min_alignment_len, args.min_reference_len)\
        .save(args.dot_filename).analyze()


if __name__ == '__main__':
    main()
