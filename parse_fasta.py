import os.path as op
import subprocess as sp


class ParseFasta:
    """utility functions to parse fasta"""

    def __init__(self, fasta_file = ''):
        if len(fasta_file) > 0:
            assert isinstance(fasta_file, str)
            assert op.isfile(fasta_file)
            assert op.isfile(fasta_file + '.fai')
            self.fasta_file = fasta_file

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