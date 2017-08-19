
class MafEntry:
    def __init__(self, genome, chrom, start, length, strand, all_length, seq):
        self.genome = genome
        self.chrom = chrom
        self.start = start
        self.length = length
        self.strand = strand
        self.all_length = all_length
        self.seq = seq
        self.block_id = 0
        if self.strand == '+':
            self.global_start = self.start
            self.global_end = self.start + self.length
        elif self.strand == '-':
            self.global_end = self.all_length - self.start
            self.global_start = self.all_length - self.start - self.length

    def set_block_id(self,block_id):
        self.block_id = block_id

    def get_block_id(self):
        return self.block_id

    def get_specie(self):
        return self.genome

    def get_seq_id(self):
        return self.chrom

    def get_chrom(self):
        return self.chrom

    def get_start(self):
        return self.start

    def get_end(self):
        return self.start + self.length

    def get_global_start(self):
         return self.global_start

    def get_global_end(self):
         return self.global_end

    def print_out(self):
        #print ' '.join(map(str,['s', self.genome + '.' + self.chrom, self.start,\
        #                        self.length, self.strand, self.all_length, self.seq]))
        print ' '.join(map(str,['s', self.genome + '.' + self.chrom, self.global_start,\
                                self.global_end, self.strand, self.all_length, self.seq]))
    def print_out_short(self):
        print ' '.join(map(str,[self.genome + '.' + self.chrom, self.global_start,\
                                self.global_end]))

    def log_out_short(self):
        return ' '.join(map(str,[self.genome + '.' + self.chrom, self.global_start,\
                                self.global_end]))
    #this prints canonical maf entry
    def print_out_local_coords(self):
        print ' '.join(map(str,['s', self.genome + '.' + self.chrom, self.start,\
                                self.length, self.strand, self.all_length, self.seq]))

class MAF_Block:
    def __init__(self,maf_id,entries,isdup=False):
        self.maf_id = maf_id
        self.entries = entries
        self.isdup = isdup
        self.entries = self.set_block_ids(maf_id)

    def set_block_ids(self,maf_id):
        for e in self.entries:
            e.set_block_id(maf_id)
        return self.entries

    def print_out(self):
        print self.maf_id
        for e in self.entries:
            e.print_out()


class MafGraph:
    def __init__(self, maf_entry, next_id, prev_id):
        self.maf_entry = maf_entry
        self.next_id = next_id
        self.prev_id = prev_id

    def log_out(self):
        prev = 'END'
        if self.prev_id != -1:
            prev = self.prev_id
        next = 'END'
        if self.next_id != -1:
            next = self.next_id
        return self.maf_entry.log_out_short() + ' PREV: '\
               + str(prev) + ' NEXT: ' + str(next) + '\n'

'''
def parse_maf(maf_file,filter_dups=False,number_genomes=0):
    maf_id = 0
    blocks = []
    with open(maf_file) as maf:
        maf_entries = []
        for line in maf:
            line = line.strip()
            if not line or '#' in line:
                continue
            if line[0] == 'a':
                if maf_entries:
                    maf_id += 1
                    isdup = False
                    genomes = (lambda x:x.genome,maf_entries)
                    if not filter_dups or\
                        filter_dups and len(genomes)==len(set(genomes))==number_genomes:
                        isdup = True
                    blocks.append(MAF_Block(maf_id,maf_entries,isdup))
                maf_entries = []
            else:
                line = line.split()
                line = line[1:]
                genome = line[0].split('.')[0]
                chrom = '.'.join(line[0].split('.')[1:])
                maf_entries.append(MAF_Entry(genome, chrom, int(line[1]), int(line[2]),\
                                             line[3], int(line[4]), line[5]))
    if maf_entries:
        genomes = (lambda x:x.genome,maf_entries)
        isdup = False
        if not filter_dups or\
            filter_dups and len(genomes) == len(set(genomes)) == number_genomes:
                isdup = True
        blocks.append(MAF_Block(maf_id+1,maf_entries,isdup))
    return blocks

def check_maf_for_no_overlaps(maf_entries):
    #i don't keep any storage for processed entries
    #in order not to load the whole maf into memory
    for e1 in maf_entries:
        for e2 in maf_entries:
            if e1 == e2 or e1.genome != e2.genome or e1.chrom != e2.chrom:
                continue
            if e1.start < e2.start and e1.start + e1.length > e2.start:
                print 'e2 overlaps e1:'
                e1.print_out_local_coords()
                e2.print_out_local_coords()
            if e1.global_start < e2.global_start and e1.global_end > e2.global_start:
                print 'e2 overlaps e1 global:'
                e1.print_out()
                e2.print_out()
            if e2.start < e1.start and e2.start + e2.length > e1.start:
                print 'e1 overlaps e2:'
                e1.print_out_local_coords()
                e2.print_out_local_coords()
            if e2.global_start < e1.global_start and e2.global_end > e1.global_start:
                print 'e1 overlaps e2 global:'
                e1.print_out()
                e2.print_out()

'''
