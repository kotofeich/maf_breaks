#!/usr/bin/env python

from argparse import ArgumentParser
import core
import utils

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('maf', help='path to maf file')
    parser.add_argument('specie1', help='reference specie name')
    parser.add_argument('specie2', help='target specie name')
    parser.add_argument('--prefix', nargs='?', help='prefix to output files')
    args = parser.parse_args()
    print utils.get_time()
    print 'parsing maf'
    maf_entries_sp1, maf_entries_sp2 = core.parse_maf(args.maf, args.specie1, args.specie2)
    print utils.get_time()
    print 'building graph from blocks', args.specie1
    graph_sp1 = core.build_graph_from_blocks(maf_entries_sp1)
    print utils.get_time()
    print 'building graph from blocks', args.specie2
    graph_sp2 = core.build_graph_from_blocks(maf_entries_sp2)
    print utils.get_time()
    print 'processing breakpoints'
    prefix = ''
    if args.prefix:
        prefix = args.prefix
    core.output_breakpoints_by_maf_id(graph_sp1, graph_sp2, args.specie1, args.specie2, prefix)
    print utils.get_time()
    print 'storing graph information for', args.specie1
    core.store_graph(graph_sp1, args.specie1, prefix)
    print utils.get_time()
    print 'storing graph information for', args.specie2
    core.store_graph(graph_sp2, args.specie2, prefix)
