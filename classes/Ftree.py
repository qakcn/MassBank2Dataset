import math

class Ftree:
    fragintreecnt = {}

    linkcnt = {}
    fraginlinkcnt = {}
    alllinkcnt = 0

    ftrees = []
    ftrees_idx = []

    @staticmethod
    def gen_ftree(rootid):
        ftree = Ftree(rootid)
        Ftree.ftrees.append(ftree)
        Ftree.ftrees_idx.append(rootid)
        return ftree

    @staticmethod
    def get_all_ftrees():
        return Ftree.ftrees

    @staticmethod
    def save():
        return dict(
            fragintreecnt=Ftree.fragintreecnt,
            linkcnt=Ftree.linkcnt,
            fraginlinkcnt=Ftree.fraginlinkcnt,
            alllinkcnt=Ftree.alllinkcnt,
            ftrees=Ftree.ftrees
        )
    
    def del_ftree(rootid):
        if rootid in Ftree.ftrees_idx:
            idx = Ftree.ftrees_idx.index(rootid)
            del Ftree.ftrees[idx]
            del Ftree.ftrees_idx[idx]

    @staticmethod
    def load(data):
        Ftree.fragintreecnt = data['fragintreecnt']
        Ftree.linkcnt = data['linkcnt']
        Ftree.fraginlinkcnt = data['fraginlinkcnt']
        Ftree.alllinkcnt = data['alllinkcnt']
        Ftree.ftrees = data['ftrees']

    def __init__(self, rootid):
        self.rootid = rootid
        self.root_in_idx = 0
        self.non_root_frags = []
        
        self.nodes = {
            'id': [],
            'out_idx': [],
            'role': [],
            'attr': [],
        }
        self.edges = {
            'id': [],
            'from_in_idx': [],
            'to_in_idx': [],
            'from_fragid': [],
            'to_fragid': [],
            'attr': [],
            'weight': [],
        }

    def add_frag(self, out_idx, fragid, is_root=False, **kwargs):
        self.nodes['id'].append(fragid)
        self.nodes['out_idx'].append(out_idx)
        self.nodes['attr'].append(kwargs)
        in_idx = self.get_frag_in_idx(out_idx)
        if is_root:
            self.root_in_idx = in_idx
        else:
            self.non_root_frags.append(in_idx)
            Ftree.increase_counter('fragintree', fragid)


    def get_fragid(self, idx, out=True):
        in_idx = self.get_frag_in_idx(idx) if out else idx
        return self.nodes['id'][in_idx]

    def get_frag_in_idx(self, out_idx):
        return self.nodes['out_idx'].index(out_idx)

    def get_frag_attr(self, in_idx):
        return self.nodes['attr'][in_idx]

    def add_link(self, from_out_idx, to_out_idx, from_fragid, to_fragid, linkid, **kwargs):
        from_in_idx = self.get_frag_in_idx(from_out_idx)
        to_in_idx = self.get_frag_in_idx(to_out_idx)
        if from_in_idx == self.root_in_idx:
            if to_in_idx in self.non_root_frags:
                del self.non_root_frags[self.non_root_frags.index(to_in_idx)]
        else:
            Ftree.increase_counter('fraginlink', from_fragid)
        if linkid != 0:
            Ftree.increase_counter('link', linkid)
        Ftree.increase_counter('fraginlink', to_fragid)
        Ftree.alllinkcnt += 1

        self.edges['id'].append(linkid)
        self.edges['from_in_idx'].append(from_in_idx)
        self.edges['to_in_idx'].append(to_in_idx)
        self.edges['from_fragid'].append(from_fragid)
        self.edges['to_fragid'].append(to_fragid)
        self.edges['attr'].append(kwargs)
        self.edges['weight'].append(None)

    @staticmethod
    def increase_counter(ctype, cid=None):
        if ctype in ('fragintree', 'fraginlink', 'link') and cid is not None:
            exec("Ftree.{ctype}cnt[cid]=1 if cid in Ftree.{ctype}cnt else 1".format(ctype=ctype))
        elif ctype == 'alllink':
            Ftree.alllinkcnt += 1
        else:
            raise ValueError('Counter type must be \'fragintree\', \'fraginlink\', \'link\' or \'alllink\'.')

    @staticmethod
    def get_counter(ctype, cid=None):
        if ctype in ('fragintree', 'fraginlink', 'link') and cid is not None:
            return eval("Ftree.{ctype}cnt[cid]".format(ctype=ctype))
        elif ctype == 'alllink':
            return Ftree.alllinkcnt
        else:
            raise ValueError('Counter type must be \'fragintree\', \'fraginlink\', \'link\' or \'alllink\'.')

    def calc_link_weight(self):
        self.create_root_link()
        for i in range(len(self.edges['weight'])):
            from_in_idx = self.edges['from_in_idx'][i]
            to_in_idx = self.edges['to_in_idx'][i]
            to_fragid = self.edges['to_fragid'][i]

            if from_in_idx == self.root_in_idx:
                weight = self.calc_rootlink_weight(to_in_idx, to_fragid)
            else:
                from_fragid = self.edges['from_fragid'][i]
                linkid = self.edges['id'][i]
                weight = self.calc_fraglink_weight(from_fragid, to_fragid, linkid)

            self.edges['weight'][i]=weight

    def create_root_link(self):
        from_out_idx = self.nodes['out_idx'][self.root_in_idx]
        for i in self.non_root_frags:
            to_out_idx = self.nodes['out_idx'][i]
            from_fragid = self.nodes['id'][i]
            to_fragid = self.nodes['id'][i]
            self.add_link(from_out_idx, to_out_idx, from_fragid, to_fragid, 0)

    def calc_fraglink_weight(self, from_fragid, to_fragid, linkid):
        w = Ftree.get_counter('alllink')
        wi = Ftree.get_counter('fraginlink', from_fragid)
        wj = Ftree.get_counter('fraginlink', to_fragid)
        pi = wi / w
        pj = wj / w
        wij = Ftree.get_counter('link', linkid)
        pij = wij / w
        return math.log10(pij / (pi * pj))

    def calc_rootlink_weight(self, frag_in_idx, fragid):
        fij = self.get_frag_attr(frag_in_idx)['relativeIntensity']
        s = len(Ftree.ftrees)
        fj = Ftree.get_counter('fragintree', fragid)
        return fij * math.log10(s / fj)

    def get_data(self):
        self.calc_link_weight()
        return self.nodes['id'], self.nodes['attr'], [self.edges['from_in_idx']+self.edges['to_in_idx'],self.edges['to_in_idx']+self.edges['from_in_idx']], self.edges['weight']*2