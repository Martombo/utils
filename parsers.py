import functions as fn
import os.path
import numpy as np
import subprocess as sp
import pysam as ps
import re

class Bam:
    """utility functions to parse bam"""

    def __init__(self, bam_path='', sam_data='', reads_orientation='forward'):
        """
        sam data can be provided for testing
        splices_dic contains all splice sites (key) mapped to a list with the count of the site and the intron coverage
        :param reads_orientation: either 'forward' or 'reverse'
        """
        self.reads_orientation = reads_orientation
        self.splices_dic = {}
        if bam_path:
            assert bam_path[-4:] == '.bam'
            assert os.path.isfile(bam_path)
        elif sam_data:
            bam_path = self.Sam(sam_data = sam_data).toIndexedBam()
            os.remove('tmp.sam')
        else:
            return
        if not os.path.isfile(bam_path + '.bai'):
            pIndex = sp.Popen(['samtools', 'index', bam_path])
            pIndex.communicate()
        self.bam_file = ps.AlignmentFile(bam_path,'rb')

    def delete(self):
        os.remove(self.bam_file.filename.decode())
        os.remove(self.bam_file.filename.decode() + '.bai')

    def get_coverage(self, chrom, start, stop, min_qual=40):
        """
        get the number of reads in region
        reads is counted even if only 1 base overlaps region
        :param chrom: str chromosome name
        :param start: int start
        :param stop: int stop
        :param min_qual: default TopHat: only uniquely mapped reads
        """
        fetch = self.bam_file.fetch(chrom, start, stop)
        n_reads = 0
        for read in fetch:
            if read.mapq >= min_qual:
                n_reads += 1
        return n_reads

    def get_splicing_coverage(self, min_qual = 40, max_mm=0, max_sites_within = 0, min_overlap = 10):
        """
        get the number of reads from within splice site
        :param chrom: str chromosome name
        :param start: int start
        :param stop: int stop
        :param min_qual: default TopHat: only uniquely mapped reads
        :param max_mm: max tolerated mismatches on read ( including indel
        :param: min_overlap number of nucleotides inside splice junctions
        """
        for splice_site in self.splices_dic.keys():
            (chrom, start, stop, strand) = splice_site.split('_')
            start = int(start)
            stop = int(stop)
            (count, na) = self._count_coverage(chrom, start, stop, strand, min_qual, max_mm, min_overlap)
            if na > max_sites_within:
                self.splices_dic[splice_site].append('Na')
            else:
                self.splices_dic[splice_site].append(count)

    def _count_coverage(self, chrom, start, stop, strand,  min_qual, max_mm, min_overlap):
        """
        returns the number of reads from within given region
        reads is only counted if fully inside and with less than max_mm mismatches
        """
        fetch = self.bam_file.fetch(chrom, start, stop)
        n_reads = 0
        n_na = 0
        for read in fetch:
            if read.mapq  >= min_qual and strand == self._determine_strand(read) and not self._is_around_site(read, start, stop):
                count = self._read_in_intron(read, start, stop, max_mm, min_overlap)
                if  count <= 1:
                    n_reads += count
                else:
                    n_na += 1
        return (n_reads, n_na)

    def _read_in_intron(self, read, start, stop, max_mm, min_overlap):
        """
        checks if a read is fully inside region and a good enough match
        :param read: read of interest, pysam object
        :param start: satrt of region
        :param stop: end of region
        :param max_mm: tolerate number of mismatches
        :return: is_in_intron 0 not ,1 is in site 3: Na, other splice site within site thus biased coverage
        """
        is_in_intron = 1
        leng = 0
        mm = 0
        has_intron = False

        for c_part in read.cigar:
            if c_part[0] != 1: # insertion = 1 no extension on refference
                leng += c_part[1]
            if c_part[0] != 0:
                mm += c_part[1]
            if c_part[0] == 3:
                has_intron = True
        if start > read.reference_start: # read befor start
            overlap = leng - (start - read.reference_start)
        elif stop < read.reference_start + leng : # read longer than junction
            overlap = stop -read.reference_start
        else:
            overlap = leng
        if overlap < min_overlap:
            is_in_intron = 0
        if has_intron:
            is_in_intron = 2
        return is_in_intron

    def _is_around_site(self, read, start, stop):
        read_splicer = self._read_splicer(read.cigar, read.reference_start)
        read_splice_sites = read_splicer.get_sites()
        for site in read_splice_sites:
            if site[0] <= start and site[1] >= stop:
                return True
        return False

    def get_splice_sites(self, min_qual=40):
        """
        get splice sites counts as dictionary
        :param min_qual: default TopHat only uniquely mapped reads
        :return dict: {"chr1_12038_12759_+": 56, ...}
        """
        for read in self.bam_file.fetch():
            if read.mapq < min_qual:
                continue
            read_splicer = self._read_splicer(read.cigar, read.reference_start)
            read_splice_sites = read_splicer.get_sites()
            if not read_splice_sites:
                continue
            chrom = self.bam_file.getrname(read.reference_id)
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

    class Sam:
        """
        handles Sam data or files.
        raw data should be provided only for tests.
        """

        def __init__(self, sam_data='', sam_path=''):
            assert bool(sam_data) != bool(sam_path)
            if sam_data:
                sam_path = 'tmp.sam'
                with open(sam_path, 'w') as fin:
                    fin.write(sam_data)
            self.sam_path = sam_path

        def toIndexedBam(self):
            bam_file = self._toBam()
            sorted_bam_path = self._sort(bam_file.filename.decode())
            os.remove(bam_file.filename.decode())
            return sorted_bam_path

        def _toBam(self):
            tmp_bam_path = self.sam_path + '.bam'
            pView = sp.Popen(['samtools', 'view', '-Sb', '-o', tmp_bam_path, self.sam_path])
            pView.communicate()
            return ps.AlignmentFile(tmp_bam_path, 'rb')

        def _sort(self, bam_path):
            sorted_bam_path = bam_path[:len(bam_path)-8]
            ps.sort(bam_path, sorted_bam_path)
            return sorted_bam_path + '.bam'


    class _read_splicer:
        """
        retrieves read info about splicing (junction) gaps
        :param cigar: read cigar string
        :param start: read mapping start position
        """
        def __init__(self, cigar, start):
            self.before_splice = True
            self.pos = start
            self.cigar = cigar
            self.splice_sites = []

        def get_sites(self):
            """
            computes the donor and acceptor sites (junction) of a read
            :return list of donor and acceptor position: [(152683, 153107), (153194, 153867)]
                     None if read overlaps no junction
            """
            for cigar_part in self.cigar:
                if cigar_part[0] == 0:
                    self._match(cigar_part[1])
                elif cigar_part[0] == 3:
                    self._non_match(cigar_part[1])
                elif cigar_part[0] == 2:
                    self._move_pos(cigar_part[1])
            return self.splice_sites

        def _match(self, leng):
            if not self.before_splice:
                self.splice_sites.append((self.start_site,self.pos))
                self.before_splice = True
            self._move_pos(leng)

        def _non_match(self, leng):
            if self.before_splice:
                self.before_splice = False
                self.start_site = self.pos
            self._move_pos(leng)

        def _move_pos(self, leng):
            self.pos += leng


class Gtf:
    """utility functions to parse gtf"""

    def __init__(self, gtf_file, test=False):
        if not test:
            assert isinstance(gtf_file, str)
            assert os.path.isfile(gtf_file)
        self.gtf_path = gtf_file

    def get_trans_exon(self, feature='exon', annot='havana'):
        """gets transcript exons coordinates
        gtf must be sorted by transcript exons number
        first exon is first also for reverse strand
        :returns {'ENST000024631: [['chr1', 12623, 13486, '+'], [...]]}
        """
        trans_exons = {}
        for linea in open(self.gtf_path).readlines():
            dic = self._parse_line(linea, feature=feature, annot=annot)
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

    def get_start_stop_codons(self, annot='havana'):
        """gets start and stop codons position
        start codon is always first"""
        trans_start_stop = {}
        for linea in open(self.gtf_path).readlines():
            dic = self._parse_line(linea, feature='start_codon|stop_codon', annot=annot)
            if not dic:
                continue
            trans_start_stop = fn.add2dict(trans_start_stop, dic['attr']['transcript_id'], [dic['start'], dic['stop']])
        return trans_start_stop

class Fasta:
    """utility functions to parse fasta"""

    def __init__(self, fasta_file, test=False):
        if not test:
            assert isinstance(fasta_file, str)
            assert os.path.isfile(fasta_file)
            assert os.path.isfile(fasta_file + '.fai')
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
        """finds the sequence of a genomic interval '1:387941-388099', '+' """
        assert self._has_fasta()
        samtools_cmm = ['samtools', 'faidx', self.fasta_file, chr_pos]
        samtools_out = sp.check_output(samtools_cmm).decode()
        seq = self.fasta2seq(samtools_out)
        if strand == '-':
            seq = fn.comp_rev(seq)
        return seq.upper()

    def fasta2seq(self, fasta):
        """extracts the fasta sequence from a multi-line fasta format"""
        assert fasta.startswith('>')
        samtools_splat = fasta.split('\n')
        return ''.join(samtools_splat[1:])

    def _has_fasta(self):
        return len(self.fasta_file) > 0


class Bed:
    """utility functions to parse bed"""

    def __init__(self, bed_file, test=False):
        if not test:
            assert isinstance(bed_file, str)
            assert os.path.isfile(bed_file)
        self.bed_path = bed_file

    def get_first(self):
        """get 1-based start (gtf-like) and score of bed intervals
        :returns np.array: chr,pos,strand,score
        """
        chr_pos_score = []
        dtype = [('chr', 'U15'), ('pos', int), ('strand', 'U1'), ('score', int)]
        for linea in open(self.bed_path).readlines():
            dic = self._parse_tophat(linea)
            if dic:
                chr_pos_score.append((dic['chr'], dic['first'], dic['strand'], dic['score']))
        return np.array(chr_pos_score, dtype)

    def get_intervals(self):
        """get 1-based start (gtf-like) and score of bed intervals
        :returns np.array: chr,pos,strand,score
        """
        chr_pos_score = []
        dtype = [('chr', 'U15'), ('start', int), ('stop', int), ('strand', 'U1')]
        for linea in open(self.bed_path).readlines():
            dic = self._parse_line6(linea)
            if dic:
                chr_pos_score.append((dic['chr'], dic['start'], dic['stop'], dic['strand']))
        return np.array(chr_pos_score, dtype)

    def _parse_tophat(self, linea):
        """parses a bed line, with 5 columns, strand is before score
        :returns {chr, start, stop, score, strand, first_base}
        """
        splat = linea.rstrip('\n').split('\t')
        if len(splat) < 5:
            return None
        first = int(splat[1]) + 1 if splat[3] == '+' else int(splat[2])
        return dict(chr=splat[0], start=int(splat[1]), stop=int(splat[2]), score=int(splat[4]),
                    strand=splat[3], first=first)

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

class Maf:
    """utility functions to parse maf
    ... it appears to be 0 based ...
    """

    def __init__(self, maf_file, main_genome, genomes, test=False):
        if not test:
            assert isinstance(maf_file, str)
            assert os.path.isfile(maf_file)
        self.maf_file = maf_file
        self.align_dict = {}
        self.main_genome = main_genome
        self.genomes = genomes
        for genome in [main_genome, 'pos'] + genomes:
            self.align_dict[genome] = []

    def get_region(self, start, stop, strand='+'):
        """returns a dict of seqs of a region of the alignment
        the interval is a closed one
        :returns genome_align: {'dm6': 'ACTGTAA---GCTAG'}
        """
        genome_align = self._scroll_region(start, stop)
        genome_align = self._list2seq(genome_align)
        if strand == '-':
            genome_align = self._rev_comp(genome_align)
        return genome_align

    def _scroll_region(self, start, stop):
        k = 0
        pos = self.align_dict['pos'][k]
        genome_align = {x:[] for x in [self.main_genome] + self.genomes}
        start = start - 1
        while pos < start:
            (k, pos) = self._advance_pos(k, pos)
        while pos < stop:
            for genome_i in [self.main_genome] + self.genomes:
                genome_align[genome_i].append(self.align_dict[genome_i][k-1])
            (k, pos) = self._advance_pos(k, pos)
        return genome_align

    def _list2seq(self, genome_align):
        for llave, lista in genome_align.items():
            genome_align[llave] = ''.join(lista).upper()
        return genome_align

    def _rev_comp(self, genome_align):
        for llave, seq in genome_align.items():
            genome_align[llave] = fn.comp_rev(seq)
        return genome_align

    def _advance_pos(self, index, pos):
        pos_or_NA = self.align_dict['pos'][index]
        if pos_or_NA != 'NA':
            pos = pos_or_NA
        index += 1
        return (index, pos)

    def get_alignments(self):
        """stores in align_dict all MAF data of specified genomes
        only one chrom of main_genome must be present
        align_dict is a dict, with genome keys: (eg: 'dm6', 'droSec1')
        which has a list: [['A', 1732], ['-', 'NA'], ['C', 1733], ...]
        careful for memory usage!
        """
        with open(self.maf_file) as fin:
            for linea in fin.readlines():
                align_line = self._parse_line(linea)
                if not align_line or align_line['genome'] not in self.genomes + [self.main_genome]:
                    continue
                self._add2dic(align_line)
        self._fill_missing()

    def _add2dic(self, line):
        if line['genome'] == self.main_genome:
            self._fill_missing()
        pos = line['start']
        for k in line['seq']:
            self.align_dict[line['genome']].append(k)
            if line['genome'] == self.main_genome:
                pos = self._add_pos2dic(pos, k)

    def _add_pos2dic(self, pos, k):
        if k == '-':
            self.align_dict['pos'].append('NA')
        else:
            self.align_dict['pos'].append(pos)
            pos += 1
        return pos

    def _fill_missing(self):
        ref_len = len(self.align_dict[self.main_genome])
        for genome in self.genomes:
            while len(self.align_dict[genome]) < ref_len:
                self.align_dict[genome].append('-')

    def _parse_line(self, line):
        splat = [x for x in line.rstrip('\n').split(' ') if x]
        if not splat or splat[0] != 's':
            return None
        genome_chr = splat[1].split('.')
        assert len(splat) == 7
        assert len(genome_chr) == 2
        return dict(genome=genome_chr[0], chr=genome_chr[1], start=int(splat[2]), strand=splat[4], seq=splat[6])
