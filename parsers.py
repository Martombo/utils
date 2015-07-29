import os.path as op
import numpy as np
import subprocess as sp
import pysam as ps
import re

class Bam:
    """utility functions to parse bam"""

    def __init__(self, bam_file, reads_orientation = 'forward'):
        """
        bam file must be indexed
        :param reads_orientation: either 'forward' or 'reverse'
        """
        assert op.isfile(bam_file)
        assert op.isfile(bam_file + '.bai')
        self.bam_file = ps.AlignmentFile(bam_file,'rb')
        self.reads_orientation = reads_orientation

    def get_coverage(self, chrom, start, stop, min_qual=40):
        """
        get the number of reads in region
        :param min_qual: default only uniquely mapped reads
        """
        fetch = self.bam_file.fetch(chrom, start, stop)
        return len(fetch)                               # does it work?

    def get_splice_sites(self, min_qual=40):
        """
        get splice sites counts as dictionary
        :param min_qual: default only uniquely mapped reads
        :return: dict: {"chr1_12038_12759_+": 56, ...}
        """
        self.splices_dic = {}
        for read in self.bam_file.fetch():
            if read.mapq < min_qual:
                continue
            read_splicer = self._read_splicer(read)
            read_splice_sites = read_splicer.get_sites()
            if not read_splice_sites:
                continue
            chrom = self.sam_file.getrname(read.reference_id)
            for read_splice_site in read_splice_sites:
                strand = self._determine_strand(read)
                self._add2dict(read_splice_site, chrom, strand)
        return self.splices_dic


    def _add2dict(self, splice_site, chrom, strand):
        locus_string = '_'.join([str(x) for x in [chrom, splice_site[0], splice_site[1], strand]])
        if locus_string in self.splices_dic:
            self.splices_dic[locus_string] += 1
        else:
            self.splices_dic[locus_string] = 1


    def _determine_strand(self, read):
        strand_bool = True
        if read.is_reverse:
            strand_bool = not strand_bool
        if self.reads_orientation == 'reverse':
            strand_bool = not strand_bool
        if read.is_read2:
            strand_bool = not strand_bool
        return '+' if strand_bool else '-'

    class _read_splicer:
        """
        handles read info about splicing (junction) gaps
        """
        def __init__(self, read):
            self.before_splice = True
            self.pos = read.reference_start
            self.splice_sites = []
            self.start_site = 0

        def get_sites(self):
            """
            computes the donor and acceptor sites (junction) of a read
            :param read: a pysam read
            :return: list of donor and acceptor position: [[152683, 153107],[153194, 153867]]
                     None if read overlaps no junction
            """
            for cigar_part in self.read.cigar:
                if cigar_part[0] == 0:          # match
                    if not self.before_splice:
                        self.splice_sites.append((self.start_site,self.pos))
                        self.before_splice = True
                    self.pos += cigar_part[1]
                elif cigar_part[0] == 3:        # non-match
                    if self.before_splice:
                        self.before_splice = False
                        start_site = self.pos
                    self.pos += cigar_part[1]
                elif cigar_part[0] == 2:        # deletion
                    self.pos += cigar_part[1]
            return self.splice_sites


class Gtf:
    """utility functions to parse gtf"""

    def __init__(self, gtf_file):
        assert isinstance(gtf_file, str)
        assert op.isfile(gtf_file)
        self.gtf_path = gtf_file

    def get_trans_exon(self):
        """gets transcript exons coordinates
        gtf must be sorted by transcript exons number
        :returns {'ENST000024631: [['chr1', '12623', '13486', '+'], [...]]}
        """
        trans_exons = {}
        for linea in open(self.gtf_path).readlines():
            dic = self._parse_line(linea, feature='exon', annot='havana')
            if not dic:
                continue
            chr_pos = [dic['chr'], dic['start'], dic['stop'], dic['strand']]
            assert 'transcript_id' in dic['attr']
            trans_id = dic['attr']['transcript_id']
            trans_exons = self._add_exon(trans_exons, trans_id, chr_pos)
            assert len(trans_exons[trans_id]) == int(dic['attr']['exon_number'])
        return trans_exons

    def _parse_line(self, linea, feature='.', annot='.'):
        """parses a gtf line, with attributes (see self._parse_attributes())
        :returns {chr, annot, feat, start, ..., {attrs}}"""
        feat_re = re.compile(feature)
        annot_re = re.compile(annot)
        splat = linea.rstrip('\n').split('\t')
        if len(splat) < 8:
            return None
        if not feat_re.search(splat[2]) or not annot_re.search(splat[1]):
            return None
        attr = self._parse_attributes(splat[8])
        return dict(chr=splat[0], annot=splat[1], feat=splat[2], start=int(splat[3]), stop=int(splat[4]),
                    score=splat[5], strand=splat[6], frame=splat[7], attr=attr)

    def _add_exon(self, trans_exons, trans_id, chr_pos):
        if trans_id in trans_exons:
            trans_exons[trans_id].append(chr_pos)
            assert chr_pos[3] == trans_exons[trans_id][0][3]
        else:
            trans_exons[trans_id] = [chr_pos]
        return trans_exons

    def _parse_attributes(self, attrs):
        """parses the attributes in a dict
        :returns {'gene_id': 'ENSG000000123', 'gene_version' = '1', ...}
        """
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

    def trans_exon2np(self, trans_exon):
        """converts trans_exon (see self.get_trans_exon()) to numpy array
        which is sorted by position (chr, start)
        :param trans_exon: {'ENST000024631: [['chr1', '12623', '13486', '+'], [...]]}
        :returns np.array: chr,start,stop,strand,trans
        """
        if not trans_exon:
            trans_exon = self.get_trans_exon()
        arr = []
        for trans, exons in trans_exon.items():
            for exon in exons:
                exon.append(trans)
                arr.append(tuple(exon))
        dtype = [('chr', 'U15'), ('start', int), ('stop', int), ('strand', 'U1'), ('trans', 'U15')]
        np_arr = np.array(arr, dtype)
        np_arr = np.sort(np_arr, order=['chr', 'start'])
        return np.sort(np_arr, order=['chr', 'start'])

class Fasta:
    """utility functions to parse fasta"""

    def __init__(self, fasta_file):
        assert isinstance(fasta_file, str)
        assert op.isfile(fasta_file)
        assert op.isfile(fasta_file + '.fai')
        self.fasta_file = fasta_file

    def get_trans_seqs(self, trans_exons):
        """determines the sequence of a trans_exons dic
        as provided by Gtf.get_trans_exon()
        :param trans_exons {'ENST000024631: [['chr1', '12623', '13486', '+'], [...]]}
        :returns {'ENST000024631: 'CGATCGTTACGCGTATTAG...'}
        """
        assert self._has_fasta()
        trans_seqs = {}
        for trans in trans_exons.keys():
            seq = ''
            for exon in trans_exons[trans]:
                chr_pos = exon[0] + ':' + str(exon[1]) + '-' + str(exon[2])
                seq += self.get_fasta(chr_pos, exon[3])
            trans_seqs[trans] = seq
        return trans_seqs

    def get_fasta(self, chr_pos, strand):
        """finds the sequence of a genomic interval '1:387941-388099', '+'
        must have fasta
        """
        assert self._has_fasta()
        samtools_cmm = ['samtools', 'faidx', self.fasta_file, chr_pos]
        samtools_out = sp.check_output(samtools_cmm).decode()
        seq = self.fasta2seq(samtools_out)
        if strand == '-':
            seq = self.comp_rev(seq)
        return seq.upper()

    def fasta2seq(self, fasta):
        """extracts the fasta sequence from a multi-line fasta format"""
        assert fasta.startswith('>')
        samtools_splat = fasta.split('\n')
        return ''.join(samtools_splat[1:])

    def comp_rev(self, seq):
        """complementary reverse of a sequence (only ACTGN)"""
        seq = seq.lower()
        seq = seq[::-1]
        seq = seq.replace('a', 'T')
        seq = seq.replace('t', 'A')
        seq = seq.replace('c', 'G')
        seq = seq.replace('g', 'C')
        seq = seq.replace('n', 'N')
        return seq

    def _has_fasta(self):
        return len(self.fasta_file) > 0


class Bed:
    """utility functions to parse bed"""

    def __init__(self, bed_file):
        assert isinstance(bed_file, str)
        assert op.isfile(bed_file)
        self.bed_path = bed_file

    def get_first(self):
        """get 1-based start (gtf-like) and score of bed intervals
        :returns np.array: chr,pos,strand,score
        """
        chr_pos_score = []
        dtype = [('chr', 'U15'), ('pos', int), ('strand', 'U1'), ('score', int)]
        for linea in open(self.bed_path).readlines():
            dic = self._parse_line6(linea)
            chr_pos_score.append((dic['chr'], dic['first'], dic['strand'], dic['score']))
        return np.array(chr_pos_score, dtype)

    def _parse_line6(self, linea):
        """parses a bed line, with at least 6 columns
        :returns {chr, start, stop, name, ..., first_base}
        """
        splat = linea.rstrip('\n').split('\t')
        if len(splat) < 6:
            return None
        first = int(splat[1]) + 1 if splat[5] == '+' else int(splat[2])
        return dict(chr=splat[0], start=int(splat[1]), stop=int(splat[2]), name=splat[3], score=int(splat[4]),
                    strand=splat[5], first=first)

    def _parse_line3(self, linea):
        """parses a bed line, with at least 3 columns
        :returns {chr, start, stop, first_base}
        """
        splat = linea.rstrip('\n').split('\t')
        if len(splat) < 3:
            return None
        first = int(splat[1]) + 1 if splat[5] == '+' else int(splat[2])
        return dict(chr=splat[0], start=int(splat[1]), stop=int(splat[2]), first=first)
