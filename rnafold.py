import subprocess as sp


class HandleRNAfold:
    """handler for RNAfold programs"""

    def __init__(self):
        pass

    def get_plfolds(self, trans_seqs, wind_size=70):
        trans_folds = {}
        for trans, seq in trans_seqs.items():
            trans_folds[trans] = self.plfold(seq, wind_size)
        return trans_folds

    def plfold(self, seq, wind_size):
        len_seq = len(seq)
        wind_size = wind_size if wind_size < len_seq else len_seq
        p1 = sp.Popen(['echo', seq], stdout=sp.PIPE)
        p2 = sp.Popen(['RNAplfold', '-W', str(wind_size), '-o'], stdin=p1.stdout)
        p1.stdout.close()
        p2.communicate()
        return self.parse_plfold('plfold_basepairs', len_seq)

    def parse_plfold(self, fin, len_seq):
        fold_array = [0] * len_seq
        for linea in open(fin).readlines():
            splat = [x for x in linea.split(' ') if x]
            assert len(splat) == 3
            fold_array[int(splat[0]) - 1] += round(float(splat[2]), 3)
            fold_array[int(splat[1]) - 1] += round(float(splat[2]), 3)
        return fold_array