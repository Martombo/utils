import unittest as ut
import parsers as ps
import unittest.mock as um
import handlers as hn

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

    def test_match(self):
        pass


if __name__ == '__main__':
    ut.main()
