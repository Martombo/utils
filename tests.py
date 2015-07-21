import unittest as ut
import parsers as ps

class TestFasta(ut.TestCase):
    parser = ps.Fasta('/Users/martin/Dropbox/utils/gNome/GRCh38_2.fa')
    seq0 = 'ACTGACTG'
    seq1 = 'CGGAAGTGACCGTGGTGGCTTCAATAAATTTGGTG'
    pos1 = '16:31186802-31186836'

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

        # def test_get_trans_seq


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


class TestGtf(ut.TestCase):
    parser = ps.Gtf('/Users/martin/Dropbox/utils/genes/hg38/Homo_sapiens.GRCh38.80toy.gtf')
    attr = 'gene_id "ENSG00000223972"; gene_version "5"; gene_name "DDX11L1"; gene_source "havana"; ' \
           'gene_biotype "transcribed_unprocessed_pseudogene";'
    linea = '\t'.join(['1', 'havana', 'gene', '11869', '14409', '.', '+', '.', attr])
    linea_exon = '\t'.join(['1', 'havana', 'exon', '11869', '14409', '.', '+', '.', attr])
    linea_ens_hav = '\t'.join(['1', 'ensembl_havana', 'gene', '11869', '14409', '.', '+', '.', attr])

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
            pass

if __name__ == '__main__':
    ut.main()
