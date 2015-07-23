class Intersecter:

    def __init__(self, sites, regions, counter, strand=True, over_check = 10):
        """applies counter.countit on sites contained in regions
        a site over two overlapping regions will be assigned to both
        it's also not 100% accurate (over_check)
        :param sites: np.array or dict: chr,pos,[strand]
        :param regions: np.array or dict: chr,start,stop,[strand]
        :return: dict with ...
        """
        self.over_check = over_check
        self.strand = strand
        self.sites = sites
        self.regions = regions
        self.sit_i = self.reg_i = 0
        self.n_sites = len(sites)
        self.n_regions = len(regions)
        self.chr = sites[0]['chr']
        self._move_site()
        self._move_region()
        self.counter = counter

    def intersect(self):
        if self.strand:
            while not self._over():
                self._match()
        else:
            while not self._over():
                self._match_nostrand()

    def _over(self):
        return self.sit_i == self.n_sites or self.reg_i == self.n_regions

    def _match(self):
        self._check_chr()
        if self.site['pos'] < self.region['start']:
            self._move_site()
        elif self.site['pos'] < self.region['stop']:
            if self.site['strand'] == self.region['strand']:
                self.counter.countit(self.site['pos'], self.region['trans'])
            self._check_next_regions()
            self._move_site()
        else:
            self._move_region()

    def _same_chr(self):
        return self.site['chr'] == self.region['chr']

    def _check_chr(self):
        if self.site['chr'] != self.chr:
            while not self._same_chr() and not self._over():
                self._move_region()
        elif self.region['chr'] != self.chr:
            while not self._same_chr() and not self._over():
                self._move_site()
        self.chr = self.site['chr']

    def _move_region(self):
        self.region = self.regions[self.reg_i]
        self.reg_i += 1

    def _move_site(self):
        self.site = self.sites[self.sit_i]
        self.sit_i += 1

    def _check_next_regions(self):
        tmp_reg_i = self.reg_i
        for tmp_reg_i in range(self.sit_i, min(self.sit_i + self.over_check, self.n_sites)):
            tmp_reg = self.regions[tmp_reg_i]
            if self._included(self.site, tmp_reg):
                self.counter.countit(self.site['pos'], tmp_reg['trans'])

    def _included(self, site, region):
        return site['chr'] == region['chr'] and site['strand'] == region['strand'] \
               and region['start'] <= site['pos'] <= region['stop']


    def _match_nostrand(self):
        self._check_chr()
        if self.site['pos'] < self.region['start']:
            self._move_site()
        elif self.site['pos'] < self.region['stop']:
            self.counter.countit(self.site['pos'], self.region['trans'])
            self._check_next_regions()
            self._move_site()
        else:
            self._move_region()


class FoldsCounter():

    def __init__(self, folds, exons, span):
        self.folds = folds
        self.exons = exons
        self.span = span
        self.folds_clips = [0] * (span + 1)

    def countit(self, pos, trans_id):
        folds = self.folds[trans_id]
        exons = self.exons[trans_id]
        i = self._pos2index(exons, pos)
        if i > self.span and self.span + 201 < len(folds):
            self.folds_clips = [ folds[i-self.span-1: i+self.span]]

    def _pos2index(self, exons, pos):
        i = 0
        for exon in exons:
            steps = range(exon[1], exon[2]) if exon[3] == '+' else range(exon[2], exon[1], -1)
            for base in steps:
                if pos == base:
                    return i
                i += 1


class Random():

    def random_seq(length):
        "returns a random nucleotide sequence of the desired length"
        import random as rn

        seq = ''
        for k in range(length):
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
