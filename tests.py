import unittest as ut
import parsers as ps
import unittest.mock as um
import handlers as hn
import functions as fn
import pysam
import os

class TestFasta(ut.TestCase):
    parser = ps.Fasta('/Users/martin/Dropbox/utils/gNome/GRCh38_2.fa', test=True)
    seq0 = 'ACTGACTG'
    seq1 = 'CGGAAGTGACCGTGGTGGCTTCAATAAATTTGGTG'
    pos1 = '16:31186802-31186836'
    trans_exons = {'trans123': [['1', 101869, 101949, '+'], ['1', 120394, 120489, '+']]}
    trans_exons_strand = {'trans123': [['1', 120394, 120489, '-'], ['1', 101869, 101949, '-']]}
    seq_exp = 'TGAGGCAGGCAGATCACCTGAGGTCAGGAGTTCCAGACCAGCCTGGCCAACATGGTGAAATCTTGTCTCTCCTACAAATACTTTGTAAG\
GAATTAATAAATAAAAATGTTCTTGAAAGACAGAAATTAATATGCAGTTCATACTGTCAGAATTGCAGGCAATTTATCAAAGTCCCCT'

    def test_comp_rev(self):
        self.assertEquals('CCAAGGTT', fn.comp_rev('AACCTTGG'))

    def test_comp_rev_low(self):
        self.assertEquals('CAGTCAGT', fn.comp_rev('actgactg'))

    def test_fasta2seq(self):
        fasta = '\n'.join(['>header', self.seq0])
        self.assertEquals(self.seq0, self.parser.fasta2seq(fasta))

    def test_fasta2seq_multiline(self):
        fasta = '\n'.join(['>header', self.seq0, self.seq0])
        self.assertEquals(2 * self.seq0, self.parser.fasta2seq(fasta))

    def test_get_fasta(self):
        self.assertEquals(self.seq1, self.parser.get_fasta(self.pos1, '+'))

    def test_get_fasta_strand(self):
        self.assertEquals(fn.comp_rev(self.seq1), self.parser.get_fasta(self.pos1, '-'))

    def test_trans_seq(self):
        if os.path.isdir('/Users/martin'):
            seq_out = self.parser.get_trans_seqs(self.trans_exons)
            self.assertEquals(self.seq_exp, seq_out['trans123'])

    def test_trans_seq_strand(self):
        if os.path.isdir('/Users/martin'):
            seq_out = self.parser.get_trans_seqs(self.trans_exons_strand)
            self.assertEquals(fn.comp_rev(self.seq_exp), seq_out['trans123'])


class TestBed(ut.TestCase):
    parser = ps.Bed('prova.bed', test=True)
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
    parser = ps.Gtf('/Users/martin/Dropbox/utils/genes/Homo_sapiens.GRCh38.80toy.gtf', test=True)
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
            trans_exons = self.parser.get_trans_exon()
            trans_exoner = fn.TransExons(trans_exons)
            np_exons = trans_exoner.trans_exon2np()
            self.assertIn(20394, np_exons['start'])


class TestRnafold(ut.TestCase):

    plfold_out = '1  2  0.1' + '\n' + '1  3  0.2'
    scores_exp = [0.3, 0.1, 0.2]
    handler = hn.RnaFold()
    rnaPLfolder = handler.PlFold()
    rnaLfolder = handler.Lfold()
    seq = ''.join([x*10 for x in ['N', 'A', 'N', 'T', 'N']])

    def test_compute(self):
        output = self.rnaPLfolder.compute(self.seq)
        os.remove('plfold_basepairs')
        self.assertGreater(sum(output), 1)
        self.assertEquals(len(output), 50)

    def test_trans_fold(self):
        trans = 'trans1'
        output = self.handler.trans_plfolds({trans:self.seq})
        os.remove('plfold_basepairs')
        self.assertIn(trans, output)
        self.assertEquals(1, len(output))
        self.assertGreater(sum(output[trans]), 1)
        self.assertEquals(len(output[trans]), 50)

    def test_Lfold(self):
        seq = 'AAAACAAAAATCGATTTTTTGTTTTT'
        output = self.rnaLfolder.compute(seq)
        self.assertIsNotNone(output['fold'])
        self.assertIsNotNone(output['energy'])
        self.assertIsNotNone(output['pos'])
        self.assertEquals(len(output['fold']), 1)
        self.assertEquals(len(output['fold'][0]), len(seq))
        defold = output['fold'][0].replace('(','').replace(')','').replace('.','')
        self.assertIsNotNone(defold)

    def test_Lfold_2seqs(self):
        seq = 'CCCCAAAAAAGGGGGGTTTTTTN'
        output = self.rnaLfolder.compute(seq)
        self.assertEquals(2, len(output['fold']))


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

    def test_countit_strand2(self):
        site = {'pos': 212, 'score': 10}
        self.counter_strand.restart()
        self.counter_strand.countit(site, self.region)
        results = self.counter_strand.get_results()
        self.assertEquals(range(17,20), results[0][0])


class TestBam(ut.TestCase):

    read = pysam.AlignedSegment()
    bam_fields = ['asd', '0', 'chr1', '100', '50', '100M', '*', '0', '0', 'A'*100, 'A'*100 + '\n']
    splic_fields = ['asd', '0', 'chr1', '100', '50', '50M50N50M', '*', '0', '0', 'A'*100, 'A'*100 + '\n']
    header = '\t'.join(['@SQ','SN:chr1','LN:1000000'])
    parser = ps.Bam()
    sam_data = '\n'.join([header, '\t'.join(bam_fields)])
    splic_data = '\n'.join([header, '\t'.join(splic_fields)])
    parser_sam = ps.Bam(sam_data=sam_data)
    parser_splic = ps.Bam(sam_data=splic_data)


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
        parser = ps.Bam(reads_orientation='reverse')
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

    def test_add2dict2(self):
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

    def test_get_coverage(self):
        cov = self.parser_sam.get_coverage('chr1', 10, 1000, min_qual=0)
        self.assertEquals(cov, 1)
        self.parser_sam.delete()

    def test_get_splice_sites(self):
        splice_sites = self.parser_splic.get_splice_sites()
        self.assertEquals(1, len(splice_sites))
        self.assertTrue('chr1_149_199_+' in splice_sites)
        self.assertTrue(1, splice_sites['chr1_149_199_+'])


class TestMaf(ut.TestCase):
    dro_genomes = ['dro' + x for x in ['Bia2', 'Ele2', 'Ere2', 'Eug2', 'Moj3', 'Per1', 'Rho2',
                                       'Sec1', 'Sim1', 'Suz1', 'Tak2', 'Vir3', 'Yak3']]
    parser = ps.Maf('prova.maf', main_genome='dm6', genomes=dro_genomes, test=True)
    header = '##maf version=1 scoring=blastz\n'
    score = 'a score=1\n'
    align1 = 's dm6.chr2L        10 4 + 2 ACCG'
    maf_block = 's dm6.chr2L        442 268 + 23513712 AACCGCAAACCCAA---atcgacaatgcacgaca\n'
    maf_block +='s droSec1.super_14  39 262 +  2068291 AACCGCAAACCCGACCGAAtgccaatactcgaca\n'

    def test_parse_line_header(self):
        parsed = self.parser._parse_line(self.header)
        self.assertIsNone(parsed)

    def test_parse_line_score(self):
        parsed = self.parser._parse_line(self.score)
        self.assertIsNone(parsed)

    def test_parse_line_align1(self):
        parsed = self.parser._parse_line(self.align1)
        self.assertIsNotNone(parsed)
        self.assertEquals('dm6', parsed['genome'])
        self.assertEquals('chr2L', parsed['chr'])
        self.assertEquals(10, parsed['start'])
        self.assertEquals('+', parsed['strand'])
        self.assertEquals('ACCG', parsed['seq'])

    def test_genome_align(self):
        with um.patch('parsers.open', um.mock_open(read_data=self.maf_block), create=True) as m:
            parser = ps.Maf('', main_genome='dm6', genomes=self.dro_genomes, test=True)
            parser.get_alignments()
            len_dm6 = len(parser.align_dict['dm6'])
            for dro_genome in self.dro_genomes:
                self.assertEquals(len_dm6, len(parser.align_dict[dro_genome]))

    def test_get_region(self):
        with um.patch('parsers.open', um.mock_open(read_data=self.maf_block), create=True) as m:
            parser = ps.Maf('', main_genome='dm6', genomes=self.dro_genomes, test=True)
            parser.get_alignments()
            (start, stop) = (443, 462)
            region = parser.get_region(start, stop, '+')
            len_region = stop - start + 1
            dm6 = region['dm6']
            self.assertEquals('ACCGCAAACCCAA---atcgaca', dm6)
            self.assertEquals('ACCGCAAACCCGACCGAAtgcca', region['droSec1'])
            self.assertGreaterEqual(len(dm6), len_region)
            for k,i in region.items():
                self.assertEquals(len(dm6), len(i))
            dm6 = dm6.replace('-', '')
            self.assertEquals(len_region, len(dm6))

    def test_get_rev_region(self):
        with um.patch('parsers.open', um.mock_open(read_data=self.maf_block), create=True) as m:
            parser = ps.Maf('', main_genome='dm6', genomes=self.dro_genomes, test=True)
            parser.get_alignments()
            (start, stop) = (443, 462)
            region = parser.get_region(start, stop, '-')
            self.assertEquals('TGTCGAT---TTGGGTTTGCGGT', region['dm6'])
            self.assertEquals('TGGCATTCGGTCGGGTTTGCGGT', region['droSec1'])



class TransExons(ut.TestCase):

    def test_rel_pos_trans(self):
        trans_exon = {'trans1':[['chr1',100,200,'+']]}
        pos = [150]
        te_manager = fn.TransExons(trans_exon)
        rel_pos = te_manager.rel_pos_trans('trans1',pos)
        self.assertEquals(rel_pos, [50])

    def test_rev_rel_pos(self):
        trans_exon = {'trans1':[['chr1',100,200,'-']]}
        pos = [180]
        te_manager = fn.TransExons(trans_exon)
        rel_pos = te_manager.rel_pos_trans('trans1',pos)
        self.assertEquals(rel_pos, [20])

    def test_2exons_rel_pos(self):
        trans_exon = {'t1':[['1',100,199,'+'],['1',400,500,'+']]}
        pos = [450]
        te_manager = fn.TransExons(trans_exon)
        rel_pos = te_manager.rel_pos_trans('t1',pos)
        self.assertEquals(rel_pos, [150])

    def test_2exons_rev_pos(self):
        trans_exon = {'t1':[['1',401,500,'-'],['1',100,200,'-']]}
        pos = [150]
        te_manager = fn.TransExons(trans_exon)
        rel_pos = te_manager.rel_pos_trans('t1',pos)
        self.assertEquals(rel_pos, [150])


if __name__ == '__main__':
    ut.main()
