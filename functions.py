import numpy as np

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

    def has2_Aprime_twist(self):
        assert self.seq
        n_Aprime = self._get_Aprimes()
        if len(n_Aprime) > 1:
            return True
        return False

    def _get_Aprimes(self):
        (n_Aprime, lenAG) = [], 0
        for pos in range(0,self.len):
            if self._is_match(pos) and self.seq[pos] in ['A', 'G']:
                lenAG = self._checkTG(lenAG, pos)
            else:
                lenAG = 0
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

    def _checkTG(self, lenAG, pos):
        if lenAG == 0:
            return 1
        fold = self.fold
        if self.fold[pos] == ')':
            (pos, fold) = self._reverse(pos)
        pair_pos = self._get_pair_pos(pos, fold)
        if self._is_match(pair_pos + 1, ')', fold):
            return lenAG + 1
        return 1

    def _reverse(self, pos):
        return (self.len - pos - 1, self.fold[::-1])

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
        while opp_count < n_match + 1:
            if fold[pair_pos] == ')':
                opp_count += 1
            pair_pos -= 1
        return pair_pos

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
