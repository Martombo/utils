import os.path as op


class ParseBed:
    """utility functions to parse bed"""

    def __init__(self, bed_file):
        assert isinstance(bed_file, str)
        assert op.isfile(bed_file)
        self.bed_path = bed_file

    def get_start1_score(self):
        """get 1-based start (gtf-like) and score of bed intervals
        :return dict: ['chr1_340985_+', 24]"""
        chr_pos_score = []
        for linea in open(self.bed_path).readlines():
            splat = linea.rstrip('\n').split('\t')
            if splat[5] == '+':
                chr_pos = '_'.join([splat[0], str(int(splat[1]) + 1), splat[5]])
            else:
                chr_pos = '_'.join([splat[0], splat[2], splat[5]])
            chr_pos_score.append([chr_pos, int(splat[4])])
        return chr_pos_score
