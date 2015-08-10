import unittest as ut
import parsers as ps
import unittest.mock as um
import handlers as hn
import functions as fn
import pysam
import os

class TestFasta(ut.TestCase):
    parser = ps.Fasta('/Users/martin/Dropbox/utils/gNome/GRCh38_2.fa')
    seq0 = 'ACTGACTG'
    seq1 = 'CGGAAGTGACCGTGGTGGCTTCAATAAATTTGGTG'
    pos1 = '16:31186802-31186836'
    trans_exons = {'trans123': [['1', 101869, 101949, '+'], ['1', 120394, 120489, '+']]}
    trans_exons_strand = {'trans123': [['1', 120394, 120489, '-'], ['1', 101869, 101949, '-']]}
    seq_exp = 'TGAGGCAGGCAGATCACCTGAGGTCAGGAGTTCCAGACCAGCCTGGCCAACATGGTGAAATCTTGTCTCTCCTACAAATACTTTGTAAG\
GAATTAATAAATAAAAATGTTCTTGAAAGACAGAAATTAATATGCAGTTCATACTGTCAGAATTGCAGGCAATTTATCAAAGTCCCCT'

    def test_comp_rev(self):
        self.assertEquals('CCAAGGTT', self.parser.comp_rev('AACCTTGG'))

    def test_comp_rev_low(self):
        self.assertEquals('CAGTCAGT', self.parser.comp_rev('actgactg'))

    def test_fasta2seq(self):
        fasta = '\n'.join(['>header', self.seq0])
        self.assertEquals(self.seq0, self.parser.fasta2seq(fasta))

    def test_fasta2seq_multiline(self):
        fasta = '\n'.join(['>header', self.seq0, self.seq0])
        self.assertEquals(2 * self.seq0, self.parser.fasta2seq(fasta))

    def test_get_fasta(self):
        self.assertEquals(self.seq1, self.parser.get_fasta(self.pos1, '+'))

    def test_get_fasta_strand(self):
        self.assertEquals(self.parser.comp_rev(self.seq1), self.parser.get_fasta(self.pos1, '-'))

    def test_trans_seq(self):
        seq_out = self.parser.get_trans_seqs(self.trans_exons)
        self.assertEquals(self.seq_exp, seq_out['trans123'])

    def test_trans_seq_strand(self):
        seq_out = self.parser.get_trans_seqs(self.trans_exons_strand)
        self.assertEquals(self.parser.comp_rev(self.seq_exp), seq_out['trans123'])


class TestBed(ut.TestCase):
    parser = ps.Bed('/Users/martin/Dropbox/projects/year2/joel_clip/fold/prova.bed')
    linea = '\t'.join(['1', '10928', '14598', 'lulli', '10', '+'])
    linea_strand = '\t'.join(['1', '10928', '14598', 'lulli', '10', '-'])

    def test_parse_line6(self):
        dic = self.parser._parse_line6(self.linea)
        self.assertEquals('1', dic['chr'])
        self.assertEquals(10928, dic['start'])
        self.assertEquals(14598, dic['stop'])
        self.assertEquals('lulli', dic['name'])
        self.assertEquals(10, dic['score'])
        self.assertEquals('+', dic['strand'])
        self.assertEquals(10929, dic['first'])

    def test_parse_line6_strand(self):
        dic = self.parser._parse_line6(self.linea_strand)
        self.assertEquals(14598, dic['first'])

    def test_get_first(self):
        with um.patch('parsers.open', um.mock_open(read_data=self.linea_strand), create=True) as m:
            np_out = self.parser.get_first()
            self.assertEquals('1', np_out[0]['chr'])
            self.assertEquals(10, np_out[0]['score'])
            self.assertEquals(14598, np_out[0]['pos'])


class TestGtf(ut.TestCase):
    parser = ps.Gtf('/Users/martin/Dropbox/utils/genes/hg38/Homo_sapiens.GRCh38.80toy.gtf')
    attr = 'gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; \
gene_biotype "transcribed_unprocessed_pseudogene"; transcript_id "ENST00000249857";'
    linea = '\t'.join(['1', 'havana', 'gene', '11869', '14409', '.', '+', '.', attr])
    linea_exon = '\t'.join(['1', 'havana', 'exon', '11869', '14409', '.', '+', '.', attr])
    linea_ens_hav = '\t'.join(['1', 'ensembl_havana', 'gene', '11869', '14409', '.', '+', '.', attr])
    mocking = linea_exon + ' exon_number "1";\n' + '\t'.join(['1', 'havana', 'exon', '20394', '21489', '.', '+', '.',
                                                              attr + ' exon_number "2";\n'])

    def test_attr(self):
        dic = self.parser._parse_attributes(self.attr)
        self.assertEquals('ENSG00000223972', dic['gene_id'])
        self.assertEquals('5', dic['gene_version'])
        self.assertEquals('DDX11L1', dic['gene_name'])
        self.assertEquals('havana', dic['gene_source'])
        self.assertEquals('transcribed_unprocessed_pseudogene', dic['gene_biotype'])

    def test_line(self):
        dic = self.parser._parse_line(self.linea)
        self.assertEquals('1', dic['chr'])
        self.assertEquals('havana', dic['annot'])
        self.assertEquals('gene', dic['feat'])
        self.assertEquals(11869, dic['start'])
        self.assertEquals(14409, dic['stop'])
        self.assertEquals('+', dic['strand'])

    def test_line_feat(self):
        dic = self.parser._parse_line(self.linea, feature='gene')
        self.assertIsNotNone(dic)

    def test_line_annot(self):
        dic = self.parser._parse_line(self.linea, annot='havana')
        self.assertIsNotNone(dic)

    def test_line_not_feat(self):
        dic = self.parser._parse_line(self.linea, feature='exon')
        self.assertIsNone(dic)

    def test_line_not_annot(self):
        dic = self.parser._parse_line(self.linea, annot='ensembl')
        self.assertIsNone(dic)

    def test_line_feat_exon(self):
        dic = self.parser._parse_line(self.linea_exon, feature='exon')
        self.assertIsNotNone(dic)

    def test_line_annot2(self):
        dic = self.parser._parse_line(self.linea_ens_hav, annot='havana')
        self.assertIsNotNone(dic)

    def test_trans_exon(self):
        with um.patch('parsers.open', um.mock_open(read_data=self.mocking), create=True) as m:
            trans = 'ENST00000249857'
            trans_exons = self.parser.get_trans_exon()
            self.assertTrue(trans in trans_exons)
            self.assertEquals(1, len(trans_exons))
            self.assertEquals(2, len(trans_exons[trans]))
            self.assertEquals(14409, trans_exons[trans][0][2])

    def test_np_exon(self):
        with um.patch('parsers.open', um.mock_open(read_data=self.mocking), create=True) as m:
            trans = 'ENST00000249857'
            trans_exons = self.parser.get_trans_exon()
            np_exons = self.parser.trans_exon2np(trans_exons)


class TestRnafold(ut.TestCase):

    plfold_out = '1  2  0.1' + '\n' + '1  3  0.2'
    scores_exp = [0.3, 0.1, 0.2]
    handler = hn.RnaFold()
    seq = ''.join([x*10 for x in ['N', 'A', 'N', 'T', 'N']])

    def test_parse_plfold(self):
        with um.patch('handlers.open', um.mock_open(read_data=self.plfold_out), create=True) as m:
            scores = self.handler._parse_plfold(fin = 'mock', len_seq = len(self.scores_exp))
            for i in range(len(self.scores_exp)):
                self.assertAlmostEqual(self.scores_exp[i], scores[i])

    def test_trans_plfold(self):
        seq_scores = self.handler.plfold(self.seq)
        suma10 = sum(seq_scores[:10])
        self.assertEquals(0, suma10)
        suma = sum(seq_scores)
        self.assertLess(suma, 20)
        self.assertGreater(suma, 15)


class TestIntersecter(ut.TestCase):

    sites = [{'chr': 'X', 'pos': 1000, 'strand': '+'}]
    regions = [{'chr': 'X', 'start': 900, 'stop': 1100, 'strand': '+', 'trans': 'trans1'}]

    def test_match(self):
        dummy_counter = fn.Intersecter.Counter()
        intersecter = fn.Intersecter(self.sites, self.regions, dummy_counter)
        intersecter.intersect()
        results = dummy_counter.get_results()
        self.assertEquals(1, len(results))

    def test_nomatch(self):
        dummy_counter = fn.Intersecter.Counter()
        sites = [{'chr': 'X', 'pos': 10000, 'strand': '+'}]
        intersecter = fn.Intersecter(sites, self.regions, dummy_counter)
        intersecter.intersect()
        results = dummy_counter.get_results()
        self.assertEquals(0, len(results))

    def test_strand(self):
        dummy_counter = fn.Intersecter.Counter()
        sites = [{'chr': 'X', 'pos': 1000, 'strand': '-'}]
        intersecter = fn.Intersecter(sites, self.regions, dummy_counter)
        intersecter.intersect()
        results = dummy_counter.get_results()
        self.assertEquals(0, len(results))

    def test_chr(self):
        dummy_counter = fn.Intersecter.Counter()
        sites = self.sites + [{'chr': 'Y', 'pos': 1000, 'strand': '+'}]
        intersecter = fn.Intersecter(sites, self.regions, dummy_counter)
        intersecter.intersect()
        results = dummy_counter.get_results()
        self.assertEquals(1, len(results))

    def test_2matches(self):
        dummy_counter = fn.Intersecter.Counter()
        sites = self.sites + [{'chr': 'X', 'pos': 1010, 'strand': '+'}]
        intersecter = fn.Intersecter(sites, self.regions, dummy_counter)
        intersecter.intersect()
        results = dummy_counter.get_results()
        self.assertEquals(2, len(results))

    def test_site_move(self):
        dummy_counter = fn.Intersecter.Counter()
        sites = [{'chr': 'X', 'pos': 10, 'strand': '+'}] + self.sites
        intersecter = fn.Intersecter(sites, self.regions, dummy_counter)
        intersecter.intersect()
        results = dummy_counter.get_results()
        self.assertEquals(1, len(results))

    def test_region_move(self):
        dummy_counter = fn.Intersecter.Counter()
        regions = [{'chr': 'X', 'start': 90, 'stop': 100, 'strand': '+', 'trans': 'trans1'}] + self.regions
        intersecter = fn.Intersecter(self.sites, regions, dummy_counter)
        intersecter.intersect()
        results = dummy_counter.get_results()
        self.assertEquals(1, len(results))

    def test_chr_switch(self):
        dummy_counter = fn.Intersecter.Counter()
        sites = [{'chr': 'A', 'pos': 10, 'strand': '+'}] + self.sites
        regions = [{'chr': 'A', 'start': 90, 'stop': 100, 'strand': '+', 'trans': 'trans1'}] + self.regions
        intersecter = fn.Intersecter(sites, regions, dummy_counter)
        intersecter.intersect()
        results = dummy_counter.get_results()
        self.assertEquals(1, len(results))

class TestFoldsCounter(ut.TestCase):

    exons = [['1', 1, 10, '+'], ['1', 101, 110, '+']]
    exons_strand = [['1', 211, 230, '-'], ['1', 1, 10, '-']]
    trans_exons = {'trans1': exons}
    trans_exons_strand = {'trans1': exons_strand}
    trans_folds = {'trans1': [0,0,0,1,1,1,1,1,0,0,0,0,1,1,1,1,1,0,0,0]}
    trans_folds_strand = {'trans1': range(30)}
    counter = fn.Intersecter.FoldsCounter(trans_folds, trans_exons, 1)
    counter_strand = fn.Intersecter.FoldsCounter(trans_folds_strand, trans_exons_strand, 1)
    region = {'trans': 'trans1'}

    def test_pos2index0(self):
        i = self.counter._pos2index(self.exons, 1)
        self.assertEquals(0, i)

    def test_pos2index9(self):
        i = self.counter._pos2index(self.exons, 10)
        self.assertEquals(9, i)

    def test_pos2index15(self):
        i = self.counter._pos2index(self.exons, 106)
        self.assertEquals(15, i)

    def test_pos2index0_strand(self):
        i = self.counter._pos2index(self.exons_strand, 230)
        self.assertEquals(0, i)

    def test_pos2index29_strand(self):
        i = self.counter._pos2index(self.exons_strand, 1)
        self.assertEquals(29, i)

    def test_countit_first(self):
        site = {'pos': 2, 'score': 10}
        self.counter.restart()
        self.counter.countit(site, self.region)
        results = self.counter.get_results()
        self.assertEquals([0,0,0], results[0][0])
        self.assertEquals(10, results[0][1])

    def test_countit_last(self):
        site = {'pos': 109, 'score': 10}
        self.counter.restart()
        self.counter.countit(site, self.region)
        results = self.counter.get_results()
        self.assertEquals([0,0,0], results[0][0])

    def test_countit_middle(self):
        site = {'pos': 107, 'score': 10}
        self.counter.restart()
        self.counter.countit(site, self.region)
        results = self.counter.get_results()
        self.assertEquals([1,1,0], results[0][0])

    def test_countit_strand(self):
        site = {'pos': 229, 'score': 10}
        self.counter_strand.restart()
        self.counter_strand.countit(site, self.region)
        results = self.counter_strand.get_results()
        self.assertEquals(range(3), results[0][0])

    def test_countit_strand(self):
        site = {'pos': 212, 'score': 10}
        self.counter_strand.restart()
        self.counter_strand.countit(site, self.region)
        results = self.counter_strand.get_results()
        self.assertEquals(range(17,20), results[0][0])


class TestBam(ut.TestCase):

    bam_file = 'example.bam'
    assert os.path.isfile(bam_file)
    parser = ps.Bam(bam_file)
    read = pysam.AlignedSegment()
    bam_fields = ['asd', '0', 'chr1', '100', '50', '100M', '*', '0', '0', 'A'*100, 'A'*100 + '\n']

    def test_det_strand(self):
        self.read.is_reverse = False
        self.read.is_read2 = False
        strand = self.parser._determine_strand(self.read)
        self.assertEquals('+', strand)

    def test_det_strand_rev(self):
        self.read.is_reverse = True
        self.read.is_read2 = False
        strand = self.parser._determine_strand(self.read)
        self.assertEquals('-', strand)

    def test_det_strand_mate(self):
        self.read.is_reverse = False
        self.read.is_read2 = True
        strand = self.parser._determine_strand(self.read)
        self.assertEquals('-', strand)

    def test_det_strand_rev_mate(self):
        self.read.is_reverse = True
        self.read.is_read2 = True
        strand = self.parser._determine_strand(self.read)
        self.assertEquals('+', strand)

    def test_det_strand_opp(self):
        parser = ps.Bam(self.bam_file, reads_orientation='reverse')
        self.read.is_reverse = False
        self.read.is_read2 = False
        strand = parser._determine_strand(self.read)
        self.assertEquals('-', strand)

    def test_add2dict(self):
        self.parser._add2dict([100,200], 'chr1', '+')
        dic = self.parser.splices_dic
        self.assertEquals(1, len(dic))
        self.assertEquals(1, dic['chr1_100_200_+'])
        self.parser.splices_dic = {}

    def test_add2dict(self):
        self.parser._add2dict([300,500], 'chr2', '-')
        self.parser._add2dict([300,500], 'chr2', '-')
        dic = self.parser.splices_dic
        self.assertEquals(2, dic['chr2_300_500_-'])
        self.parser.splices_dic = {}

    def test_splicer_get_sites(self):
        cigar = [[0, 10], [3, 100], [0, 10]]
        reference_start = 90
        splicer = ps.Bam._read_splicer(cigar, reference_start)
        sites = splicer.get_sites()
        self.assertEquals(1, len(sites))
        self.assertEquals(2, len(sites[0]))
        self.assertEquals(100, sites[0][0])
        self.assertEquals(200, sites[0][1])

    def test_get_splice_sites(self):
        splice_sites = self.parser.get_splice_sites()
        self.assertEquals(3, len(splice_sites))
        self.assertTrue('1_12227_12612_-' in splice_sites)
        self.assertTrue(1, splice_sites['1_12227_12612_-'])

    def test_get_coverage(self):
        cov = self.parser.get_coverage('1', 10000, 11000, min_qual = 0)
        self.assertEquals(cov, 1)

    def test_get_coverage_min(self):
        cov = self.parser.get_coverage('1', 10000, 10020, min_qual = 0)
        self.assertEquals(cov, 1)


if __name__ == '__main__':
    ut.main()
