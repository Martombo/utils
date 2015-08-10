import subprocess as sp

class RnaFold:
    """handler for RNAfold programs"""

    def __init__(self, msg = None):
        if msg:
            self.msg = msg
        else:
            self.msg = lambda x: print(x)

    def trans_plfolds(self, trans_seqs, wind_size=70):
        """determines plfold scores of trans_seq
        :param trans_seqs: {'ENST000024631': 'CGATCGTTACGCGTATTAG...'}
        :returns {'ENST000024631: [0.123, 0.352, ...]}
        """
        trans_folds = {}
        k = 0
        self.msg('total transcripts: ' + str(len(trans_seqs)))
        for trans, seq in trans_seqs.items():
            trans_folds[trans] = self.plfold(seq, wind_size)
            k += 1
            if k % 1000 == 0:
                self.msg('folded ' + str(k) + ' transcripts')
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
            fold_array = self._add_score(fold_array, splat[0], splat[2])
            fold_array = self._add_score(fold_array, splat[1], splat[2])
        return [round(x, 3) for x in fold_array]

    def _add_score(self, fold_array, pos_str, score):
        pos = int(pos_str) - 1
        if 0 <= pos < len(fold_array):
            fold_array[pos] += float(score)
        return fold_array
