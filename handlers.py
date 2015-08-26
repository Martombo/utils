import subprocess as sp
import numpy as np

class RnaFold:
    """handler for RNAfold programs"""

    def __init__(self, msg=None):
        if msg:
            self.msg = msg
        else:
            self.msg = lambda x: None

    def trans_plfolds(self, trans_seqs, wind_size=70):
        """determines plfold scores of trans_seq
        :param trans_seqs: {'ENST000024631': 'CGATCGTTACGCGTATTAG...'}
        :returns {'ENST000024631: [0.123, 0.352, ...]}
        """
        trans_folds = {}
        k = 0
        self.msg('total transcripts: ' + str(len(trans_seqs)))
        for trans, seq in trans_seqs.items():
            plfold = self.PlFold()
            trans_folds[trans] = plfold.compute(seq)
            k += 1
            if k % 1000 == 0:
                self.msg('folded ' + str(k) + ' transcripts')
        return trans_folds

    class FoldProgram:

        def __init__(self, program, options=[]):
            self._set_program(program)
            self.options = options

        def _set_program(self, program):
            self.program = program

        def _run(self, seq):
            p1 = sp.Popen(['echo', seq], stdout=sp.PIPE)
            p2 = sp.Popen([self.program] + self.options, stdin=p1.stdout, stdout=sp.PIPE)
            p1.stdout.close()
            out_err = p2.communicate()
            return out_err[0].decode()

        def _parse_output(self, output):
            return output

        def compute(self, seq):
            output = self._run(seq)
            return self._parse_output(output)

    class PlFold(FoldProgram):

        def __init__(self, options=[]):
            if '-o' not in options:
                options.append('-o')
            super().__init__('RNAplfold', options)

        def compute(self, seq):
            self.len_seq = len(seq)
            return super().compute(seq)

        def _parse_output(self, output):
            fold_array = [0] * self.len_seq
            for linea in open('plfold_basepairs').readlines():
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

    class Lfold(FoldProgram):

        def __init__(self, options=[]):
            super().__init__('RNALfold', options)

        def _parse_output(self, output):
            fold_en_pos = []
            for linea in output.split('\n'):
                linea = linea.replace('( ', '(')
                splat = [x for x in linea.split(' ') if x]
                if len(splat) == 3:
                    energy = float(splat[1].replace('(','').replace(')',''))
                    fold_en_pos.append((splat[0], energy, int(splat[2])))
            dtype = [('fold', 'U200'), ('energy', float), ('pos', int)]
            return np.array(fold_en_pos, dtype)
