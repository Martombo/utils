import os.path as op

class ParseGtf:
    """utility functions to parse gtf"""

    def __init__(self, gtf_file):
        assert isinstance(gtf_file, str)
        assert op.isfile(gtf_file)
        self.gtf_path = gtf_file

    def get_trans_exon(self):
        """ gets transcript exons coordinates
        gtf must be sorted by transcript exons
        :return {'ENST000024631: [['chr1', '12623', '13486', '+']]}"""
        trans_exons = {}
        for linea in open(self.gtf_path).readlines():
            splat = linea.rstrip('\n').split('\t')
            if len(splat) < 8 or splat[2] != 'exon':
                continue
            attr = self.parse_attributes(splat[8])
            if 'transcript_id' not in attr:
                continue
            trans_id = attr['transcript_id']
            chr_pos = [splat[0], splat[3], splat[4], splat[6]]
            trans_exons = self.add_exon(trans_exons, trans_id, chr_pos)
            assert len(trans_exons[trans_id]) == int(attr['exon_number'])
        return trans_exons

    def add_exon(self, trans_exons, trans_id, chr_pos):
        if trans_id in trans_exons:
            trans_exons[trans_id].append(chr_pos)
            assert chr_pos[3] == trans_exons[trans_id][0][3]
        else:
            trans_exons[trans_id] = [chr_pos]
        return trans_exons

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
