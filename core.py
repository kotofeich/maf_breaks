from collections import defaultdict
import os

from model import MafGraph
from model import MafEntry

def init_graph(maf_entries_sp):
    graph = {}
    for seq in maf_entries_sp:
        for i in range(len(maf_entries_sp[seq])):
            maf_id,entry = maf_entries_sp[seq][i]
            if i == 0:
                prev = -1
            else:
                prev = maf_entries_sp[seq][i-1][0]
            if i == len(maf_entries_sp[seq]) - 1:
                next = -1
            else:
                next = maf_entries_sp[seq][i+1][0]
            graph[maf_id] = MafGraph(entry, prev, next)
    return graph

def build_graph_from_blocks(maf_entries_sp):
    for seq in maf_entries_sp:
        maf_entries_sp[seq] = sorted(maf_entries_sp[seq], key=lambda x: x[1].global_start)

    #maf_entries_sp1[scaffold10] == [(11,MafEntry()),(12,MafEntry()),(13,MafEntry())...]
    return init_graph(maf_entries_sp)

def store_graph(graph, specie, prefix=''):
    with open(os.path.join(prefix,specie+'.graph'), 'w') as f:
        f.write('#MAF_BLOCK_ID SEQ GLOBAL_START GLOBAL_END PREV_MAF_BLOCK_ID NEXT_MAF_BLOCK_ID\n')
        for id in graph:
            f.write(str(id)+' ')
            f.write(graph[id].log_out())


def output_breakpoints_by_maf_id(graph_sp1, graph_sp2, specie1, specie2, prefix=''):
    assert graph_sp1.keys() == graph_sp2.keys()
    #common_keys = graph_sp1.keys() + graph_sp2.keys()
    with open(os.path.join(prefix,'breakpoints.txt'), 'w') as f:
        for key in graph_sp1.keys():
            next_key_sp1 = graph_sp1[key].next_id
            prev_key_sp1 = graph_sp1[key].prev_id
            next_key_sp2 = graph_sp2[key].next_id
            prev_key_sp2 = graph_sp2[key].prev_id
            if next_key_sp1 != next_key_sp2 and \
                            next_key_sp1 != prev_key_sp2 and \
                                next_key_sp2 != prev_key_sp1:
                f.write('A '+str(key)+'\n')
                f.write('B '+graph_sp1[key].log_out())
                if prev_key_sp1 > -1:
                    f.write('C '+graph_sp1[prev_key_sp1].log_out())
                else:
                    f.write('D '+'end of chromosome\n')
                if prev_key_sp2 > -1:
                    f.write('E as ' + specie2 + ' ' + graph_sp1[prev_key_sp2].log_out())

                if next_key_sp1 > -1:
                    f.write('F '+graph_sp1[next_key_sp1].log_out())
                else:
                    f.write('G '+'end of chromosome\n')
                if next_key_sp2 > -1:
                    f.write('H as ' + specie2 + ' ' + graph_sp1[next_key_sp2].log_out())

                f.write('I '+graph_sp2[key].log_out())
                if prev_key_sp2 > -1:
                    f.write('J '+graph_sp2[prev_key_sp2].log_out())
                else:
                    f.write('K '+'end of chromosome\n')
                if prev_key_sp1 > -1:
                    f.write('L as ' + specie1 + ' ' + graph_sp2[prev_key_sp1].log_out())

                if next_key_sp2 > -1:
                    f.write('M '+graph_sp2[next_key_sp2].log_out())
                else:
                    f.write('N '+'end of chromosome\n')
                if next_key_sp1 > -1:
                    f.write('O as ' + specie1 + ' ' + graph_sp2[next_key_sp1].log_out())




def output_breakpoints_by_sp1():
    pass


#skips duplications!
def parse_maf(maf_file, specie1, specie2):
    maf_id = 0
    maf_entries_sp1 = defaultdict(list)
    maf_entries_sp2 = defaultdict(list)
    with open(maf_file) as maf:
        maf_entries = []
        for line in maf:
            line = line.strip()
            if not line or '#' in line:
                continue
            if line[0] == 'a':
                maf_id += 1
                if maf_entries:
                    entry_sp1 = filter(lambda x: x.get_specie()==specie1, maf_entries)
                    entry_sp2 = filter(lambda x: x.get_specie()==specie2, maf_entries)
                    #print len(entry_sp1), len(entry_sp2)
                    if len(entry_sp1) != 1 or len(entry_sp2) != 1 :
                        maf_entries = []
                        continue
                    maf_entries_sp1[entry_sp1[0].get_seq_id()].append((maf_id,entry_sp1[0]))
                    maf_entries_sp2[entry_sp2[0].get_seq_id()].append((maf_id,entry_sp2[0]))
                maf_entries = []
            else:
                line = line.split()
                line = line[1:]
                genome = line[0].split('.')[0]
                chrom = '.'.join(line[0].split('.')[1:])
                if len(line) == 4:
                    entry = MafEntry(genome, chrom, int(line[1]), int(line[2]), line[3], -1, "")
                else:
                    entry = MafEntry(genome, chrom, int(line[1]), int(line[2]),\
                                             line[3], int(line[4]), line[5])
                maf_entries.append(entry)

    if maf_entries:
        entry_sp1 = filter(lambda x: x.get_specie()==specie1, maf_entries)
        entry_sp2 = filter(lambda x: x.get_specie()==specie2, maf_entries)
        if len(entry_sp1) != 1 or len(entry_sp2) != 1 :
            return maf_entries_sp1,maf_entries_sp2
        maf_entries_sp1[entry_sp1[0].get_seq_id()].append((maf_id+1,entry_sp1[0]))
        maf_entries_sp2[entry_sp2[0].get_seq_id()].append((maf_id+1,entry_sp2[0]))
    return maf_entries_sp1,maf_entries_sp2