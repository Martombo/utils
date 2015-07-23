import subprocess as sp

class RnaFold:
    """handler for RNAfold programs"""

    def __init__(self):
        pass

    def trans_plfolds(self, trans_seqs, wind_size=70):
        """determines plfold scores of trans_seq
        :param trans_seqs: {'ENST000024631': 'CGATCGTTACGCGTATTAG...'}
        :returns {'ENST000024631: [0.123, 1499541: 0.352, ...]}
        """
        trans_folds = {}
        for trans, seq in trans_seqs.items():
            folds = self.plfold(seq, wind_size)
            trans_folds[trans] = 1
        return trans_folds

    def plfold(self, seq, wind_size=70):
        """computes plfold of sequence
        :returns sum of plfold scores array
        """
        len_seq = len(seq)
        wind_size = min(wind_size, len_seq)
        p1 = sp.Popen(['echo', seq], stdout=sp.PIPE)
        p2 = sp.Popen(['RNAplfold', '-W', str(wind_size), '-o'], stdin=p1.stdout)
        p1.stdout.close()
        p2.communicate()
        return self._parse_plfold('plfold_basepairs', len_seq)

    def _parse_plfold(self, fin, len_seq):
        fold_array = [0] * len_seq
        for linea in open(fin).readlines():
            splat = [x for x in linea.split(' ') if x]
            assert len(splat) == 3
            fold_array[int(splat[0]) - 1] += float(splat[2])
            fold_array[int(splat[1]) - 1] += float(splat[2])
        return [round(x, 3) for x in fold_array]
