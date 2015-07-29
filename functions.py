class Intersecter:

    def __init__(self, sites, regions, counter, strand=True, over_check = 10):
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


class Random():

    def random_seq(leng):
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
