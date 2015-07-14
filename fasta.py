import os.path as op
import subprocess as sp

class ParseFasta:
    """utility functions to parse fasta"""

    def __init__(self, fasta_file):
        assert isinstance(fasta_file, str)
        assert op.isfile(fasta_file)
        assert op.isfile(fasta_file + '.fai')
        self.fasta_file = fasta_file

    def get_trans_seqs(self, trans_exons):
        assert self.has_fasta()
        trans_seqs = {}
        for trans in trans_exons.keys():
            seq = ''
            for exon in trans_exons[trans]:
                chr_pos = exon[0] + ':' + exon[1] + '-' + exon[2]
                seq += self.get_fasta(chr_pos, exon[3])
            trans_seqs[trans] = seq
        return trans_seqs

    def get_fasta(self, chr_pos, strand):
        assert self.has_fasta()
        samtools_cmm = ['samtools', 'faidx', self.fasta_file, chr_pos]
        samtools_out = sp.check_output(samtools_cmm)
        seq = self.fasta2seq(samtools_out)
        if strand == '-':
            seq = self.comp_rev(seq)
        return seq.upper()

    def fasta2seq(self, fasta):
        assert fasta.startswith('>')
        samtools_splat = fasta.split('\n')
        return ''.join(samtools_splat[1:])

    def comp_rev(self, seq):
        seq = seq.lower()
        seq = seq[::-1]
        seq = seq.replace('a', 'T')
        seq = seq.replace('t', 'A')
        seq = seq.replace('c', 'G')
        seq = seq.replace('g', 'C')
        return seq

    def has_fasta(self):
        return len(self.fasta_file) > 0