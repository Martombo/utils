import numpy as np
import subprocess as sp

class Intersecter:

    def __init__(self, sites, regions, counter, strand=True, over_check=10):
        """applies counter.countit on sites contained in regions
        a site over two overlapping regions will be assigned to both
        it's also not 100% accurate (over_check)
        :param sites: np.array or dict: chr,pos,[strand]
        :param regions: np.array or dict: chr,start,stop,[strand]
        :param counter: counter class ex: FoldsCounter
        :return: dict with data defined by counter
        """
        self.over_check = over_check
        self.strand = strand
        self.sites = sites
        self.regions = regions
        self.sit_i = self.reg_i = 0
        self.site, self.region = sites[0], regions[0]
        self.n_sites = len(sites)
        self.n_regions = len(regions)
        self.chr = self.site['chr']
        assert self.chr == self.region['chr']
        self.counter = counter
        self.over = False

    def intersect(self):
        if self.strand:
            while not self.over:
                self._match()
        else:
            while not self.over:
                self._match_nostrand()

    def _match(self):
        self._check_chr()
        if self.site['pos'] < self.region['start']:
            self._move_site()
        elif self.site['pos'] < self.region['stop']:
            if self._included(self.site, self.region):
                self.counter.countit(self.site, self.region)
            self._check_next_regions()
            self._move_site()
        else:
            self._move_region()

    def _same_chr(self):
        return self.site['chr'] == self.region['chr']

    def _check_chr(self):
        if self.site['chr'] != self.chr:
            while not (self._same_chr() or self._region_last()):
                self._move_region()
        elif self.region['chr'] != self.chr:
            while not (self._same_chr() or self._site_last()):
                self._move_site()
        self.chr = self.site['chr']

    def _region_last(self):
        return self.reg_i == self.n_regions - 1

    def _site_last(self):
        return self.sit_i == self.n_sites - 1

    def _move_region(self):
        if not self._region_last():
            self.reg_i += 1
            self.region = self.regions[self.reg_i]
        else:
            self.over = True

    def _move_site(self):
        if not self._site_last():
            self.sit_i += 1
            self.site = self.sites[self.sit_i]
        else:
            self.over = True

    def _check_next_regions(self):
        for tmp_reg_i in range(self.reg_i + 1, min(self.reg_i + 1 + self.over_check, self.n_regions)):
            tmp_reg = self.regions[tmp_reg_i]
            if self._included(self.site, tmp_reg):
                self.counter.countit(self.site, tmp_reg)

    def _included(self, site, region):
        return site['chr'] == region['chr'] and site['strand'] == region['strand'] \
               and region['start'] <= site['pos'] <= region['stop']

    def _match_nostrand(self):
        self._check_chr()
        if self.site['pos'] < self.region['start']:
            self._move_site()
        elif self.site['pos'] < self.region['stop']:
            self.counter.countit(self.site, self.region)
            self._check_next_regions()
            self._move_site()
        else:
            self._move_region()

    class Counter():
        """
        dummy counter, countit() only returns pos, trans_id
        """

        def __init__(self):
            self.results = []

        def countit(self, site, region):
            self._add_results([site['pos'], region['trans']])

        def _add_results(self, lista):
            self.results.append(lista)

        def restart(self):
            self.results = []

        def get_results(self):
            return self.results

    class FoldsCounter(Counter):

        def __init__(self, folds, exons, span):
            self.folds = folds
            self.exons = exons
            self.span = span
            super().__init__()

        def countit(self, site, region):
            assert site['score']
            folds = self.folds[region['trans']]
            exons = self.exons[region['trans']]
            i = self._pos2index(exons, site['pos'])
            if i >= self.span and i + self.span < len(folds):
                folds_clips = [folds[i-self.span: i+self.span+1], site['score']]
                super()._add_results(folds_clips)

        def _pos2index(self, exons, pos):
            i = 0
            for exon in exons:
                steps = range(exon[1], exon[2] + 1) if exon[3] == '+' else range(exon[2], exon[1] - 1, -1)
                for base in steps:
                    if pos == base:
                        return i
                    i += 1


class TransExons():
    """
    manages a trans_exons dictionary (see parsers.Gtf)
    {'ENST000024631: [['chr1', 12623, 13486, '+'], [...]]}
    """
    def __init__(self, trans_exons):
        self.trans_exons = trans_exons
        self.check_start = None

    def rel_pos_trans(self, trans, poss):
        """finds relative position on trans (0 based) of a chromosomal location
        poss must be reverse sorted if trans is on - strand
        :param trans: trans id 'ENSG0000005007'
        :param poss: list sorted by position on gene [1209345, 1259054]
        :return: relative positions [0, 193, 5001]
        """
        k = 0
        rel_poss = []
        trans_exon = self.trans_exons[trans]
        strand = trans_exon[0][3]
        for exon in trans_exon:
            self._check_exon(exon,strand)
            for base in range(exon[1], exon[2] + 1)[::int(strand + '1')]:
                for pos in poss:
                    if base == pos:
                        rel_poss.append(k)
                k += 1
        self.check_start = None
        return rel_poss

    def _check_exon(self,exon,strand):
        if self.check_start:
            if strand == '+':
                assert exon[1] > self.check_start
            else:
                assert exon[1] < self.check_start
        self.check_start = exon[1]

    def trans_exon2np(self):
        """converts trans_exon to numpy array
        which is sorted by position (chr, start)
        :param trans_exon: {'ENST000024631: [['chr1', '12623', '13486', '+'], [...]]}
        :returns np.array: chr,start,stop,strand,trans
        """
        arr = []
        for trans, exons in self.trans_exons.items():
            for exon in exons:
                exon.append(trans)
                arr.append(tuple(exon))
        dtype = [('chr', 'U15'), ('start', int), ('stop', int), ('strand', 'U1'), ('trans', 'U15')]
        np_arr = np.array(arr, dtype)
        np_arr = np.sort(np_arr, order=['chr', 'start'])
        return np.sort(np_arr, order=['chr', 'start'])

class Fold():
    """manages and computes features of RNA secondary structures"""

    def __init__(self, fold, seq=''):
        assert fold.count('(') == fold.count(')')
        self.fold = fold
        self.len = len(fold)
        self.seq = seq.upper()

    def plot(self):
        assert self.seq
        seq_fold = '\n'.join(['>' + self.seq[0:15], self.seq, self.fold, ''])
        p1 = sp.Popen(['echo', seq_fold], stdout=sp.PIPE)
        p2 = sp.Popen('RNAplot', stdin=p1.stdout, stdout=sp.PIPE)
        p1.stdout.close()
        p2.communicate()

    def is_simple(self):
        closing = False
        for f in self.fold:
            if f == ')':
                closing = True
            elif closing and f == '(':
                return False
        return True

    def _get_pair(self, pos):
        assert self._is_match(pos)
        matches = 0
        for k in range(len(self.fold)):
            if self._is_match(k):
                matches += 1
            if k == pos:
                break
        rev_matches = 0
        for k in range(len(self.fold))[::-1]:
            if self._is_match(k):
                rev_matches += 1
            if rev_matches == matches:
                return k
        assert False

    def _is_symm_match(self, posA, posB):
        pairA = self._get_pair(posA)
        pairB = self._get_pair(posB)
        return posA - posB == pairB - pairA

    def _rescue_hairy_back(self, pos, lenAG):
        if lenAG != 2:
            return 0
        if pos - 5 >= 0 and self.seq[pos - 3] in ['A', 'G']:
            if self._is_match(pos - 4) and self._is_symm_match(pos - 4, pos - 2):
                if self._is_match(pos - 5) and self._is_symm_match(pos - 5, pos - 4):
                    return 3
        return 0

    def _rescue_hairy(self, pos, lenAG):
        if lenAG != 2:
            return 0
        if pos + 2 < self.len and self.seq[pos] in ['A', 'G']:
            if self._is_match(pos + 1) and self._is_symm_match(pos - 1, pos + 1):
                if self._is_match(pos + 2) and self._is_symm_match(pos + 1, pos + 2):
                    return 3
        return self._rescue_hairy_back(pos, lenAG)

    def _update_lenAG(self, pos, lenAG):
        if self.seq[pos] not in ['A', 'G']:
            return self._rescue_hairy_back(pos, lenAG)
        if not self._is_match(pos):
            return self._rescue_hairy(pos, lenAG)
        if lenAG == 0:
            return 1
        pair = self._get_pair(pos)
        prev_pair = self._get_pair(pos - 1)
        if pair - prev_pair in [1, -1]:
            return lenAG + 1
        else:
            return 1

    def get_Aprimes(self):
        assert self.seq
        (n_Aprime, lenAG) = [], 0
        for pos in range(self.len):
            lenAG = self._update_lenAG(pos, lenAG)
            if lenAG == 3:
                n_Aprime.append(pos - 2)
        return n_Aprime

    def _is_match(self, pos, match_type='', fold=''):
        if not fold:
            fold = self.fold
        if not match_type:
            match_type = '()'
        if 0 <= pos < self.len:
            if fold[pos] in match_type:
                return True
        return False

    def _reverse(self, pos):
        fold = self.fold.replace('(','+')
        fold = fold.replace(')','(')
        fold = fold.replace('+',')')
        return (self.len - pos - 1, fold[::-1])

    def _get_pair_pos(self, pos, fold=''):
        """it only works for first match '('
        for second matches ')', make a reverse of
        pos, fold and also returned value
        """
        if not fold:
            fold = self.fold
        substr = self.fold[0:pos]
        n_match = substr.count('(')
        (pair_pos, opp_count) = self.len - 1, 0
        while opp_count < n_match:
            if fold[pair_pos] == ')':
                opp_count += 1
            pair_pos -= 1
        return pair_pos


class GappedSeq():

    def get_gapped_i(self, gapped_seq, pos, leng):
        """pos and leng are indexes of an ungapped subsequence of gapped_seq
        :returns the start and stop of the subsequence in gapped_seq
        """
        (ungapped_k, start, stop) = -1, 0, 0
        for k in range(len(gapped_seq)):
            if gapped_seq[k] != '-':
                ungapped_k += 1
            if ungapped_k == pos and not start:
                start = k
            elif ungapped_k == pos + leng and not stop:
                stop = k
            if start and stop:
                return start, stop

    def get_gaps(self, seq_gap, seq):
        """inserts gaps from seq_gap to seq"""
        seq_out = ''
        for k in range(len(seq_gap)):
            if seq_gap[k] != '-':
                seq_out += seq[0]
                seq = seq[1:]
            else:
                seq_out += '.'
        assert len(seq_out) == len(seq_gap)
        return seq_out

    def get_ortho_constr(self, seqA, foldA, seqB):
        """seqA and seqB are two aligned sequences
        foldA is the ungapped folding of seqA
        :returns constraints for seqB, given fold of seqA"""
        assert len(seqA) == len(seqB)
        foldA_gap = self.get_gaps(seqA, foldA)
        foldA_gap = foldA_gap.replace('(', '<').replace(')', '>')
        foldB_out = ''
        for k in range(len(seqA)):
            if seqB[k] != '-':
                foldB_out += foldA_gap[k]
        seqB_out = seqB.replace('-','')
        assert len(seqB_out) == len(foldB_out)
        return seqB_out, foldB_out


class Random():

    def random_seq(self, leng):
        "returns a random nucleotide sequence of the desired length"
        import random as rn

        seq = ''
        for k in range(leng):
            rand = rn.random()
            if rand <= 0.25:
                next_letter = 'A'
            elif rand <= 0.5:
                next_letter = 'C'
            elif rand <= 0.75:
                next_letter = 'T'
            else:
                next_letter = 'G'

            seq += next_letter

        return seq

def add2dict(dict, key, value):
    if key in dict:
        dict[key].append(value)
    else:
        dict[key] = [value]
    return dict

def comp_rev(seq):
    """complementary reverse of a sequence (only ACTGN)"""
    seq = seq.lower()
    seq = seq[::-1]
    seq = seq.replace('a', 'T')
    seq = seq.replace('t', 'A')
    seq = seq.replace('c', 'G')
    seq = seq.replace('g', 'C')
    seq = seq.replace('n', 'N')
    return seq
