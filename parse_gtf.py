import os.path as op
import parse_fasta as pf

class ParseGtf:
    """utility functions to parse gtf"""

    def __init__(self, gtf_file):
        assert isinstance(gtf_file, str)
        assert op.isfile(gtf_file)
        self.gtf_path = gtf_file
        self.trans_exons = {}
        self.trans_seqs = {}

    def set_genome_fasta(self, genome_file):
        self.fasta_parser = pf.ParseFasta(genome_file)

    def get_folds(self):
        self.get_trans_seqs()


    def get_trans_seqs(self):
        assert self.fasta_parser.has_fasta()
        self.get_trans_exons()
        for trans in self.trans_exons.keys():
            seq = ''
            for exon in self.trans_exons[trans]:
                chr_pos = exon[0] + ':' + exon[1] + '-' + exon[2]
                seq += self.fasta_parser.get_fasta(chr_pos, exon[3])
            self.trans_seqs[trans] = seq

    def get_trans_exon(self):
        # gtf must be sorted by transcript exons
        for linea in open(self.gtf_path).readlines():
            splat = linea.rstrip('\n').split('\t')
            if len(splat) < 8 or splat[2] != 'exon':
                continue
            attr = self.parse_attributes(splat[8])
            if 'transcript_id' not in attr:
                continue
            trans_id = attr['transcript_id']
            chr_pos = [splat[0], splat[3], splat[4], splat[6]]
            self.add_exon(trans_id, chr_pos)
            assert len(self.trans_exons[trans_id]) == int(attr['exon_number'])

    def add_exon(self, trans_id, chr_pos):
        if trans_id in self.trans_exons:
            self.trans_exons[trans_id].append(chr_pos)
            assert chr_pos[3] == self.trans_exons[trans_id][0][3]
        else:
            self.trans_exons[trans_id] = [chr_pos]

    def parse_attributes(self, attrs):
        attrs_splat = attrs.split(';')
        attr_dic = {}
        for attr in attrs_splat:
            if not attr:
                continue
            attr = attr.lstrip(' ')
            attr_key = attr.split(' ')[0]
            attr_item = attr.split('"')[1]
            attr_dic[attr_key] = attr_item
        return attr_dic
