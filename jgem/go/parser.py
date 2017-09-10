""" Copyright (c) 2004, Ken Sugino

 GO (Gene Ontology) related stuffs.


"""
import os
import re
try:
    from cPickle import load, dump  # < v3
except:
    from pickle import load, dump  # > v3
try:
    from itertools import izip
except:
    izip = zip

import pandas as PD

# class Tree(object):
    
#     def __init__(self,rootid=None,parents={}):
#         self.rootid = rootid  # id is index
#         self.parents = parents  # id(idx) => parent id(idx) dict
#         self.children = children # id(idx) => [child id(idx), ...] list

class GOSlim(object):
    pbin = ['DNA metabolism',
            'RNA metabolism',
            'cell adhesion',
            'cell cycle and proliferation',
            'cell organization and biogenesis',
            'cell-cell signaling',
            'death',
            'developmental processes',
            'other biological processes',
            'other metabolic processes',
            'protein metabolism',
            'signal transduction',
            'stress response',
            'transport']
    mbin = ['bone, tooth or skin structural activity',
            'chaperone-related activity',
            'cytoskeletal activity',
            'enzyme regulator activity',
            'extracellular structural activity',
            'kinase activity',
            'nucleic acid binding activity',
            'other molecular function',
            'signal transduction activity',
            'transcription regulatory activity',
            'translation activity',
            'transporter activity']
    cbin = ['ER/Golgi',
            'cytoskeleton',
            'cytosol',
            'extracellular matrix',
            'mitochondrion',
            'non-structural extracellular',
            'nucleus',
            'other cellular component',
            'other cytoplasmic organelle',
            'other membranes',
            'plasma membrane',
            'translational apparatus']

    def __init__(self, fname='map2MGIslim.txt'):
        self.fname = fname
        self.df = df = PD.read_table(fname, header=0) # 40809 entries as of 2014-01-01
        # self.s2e = s2e = {s:i for i,s in enumerate(self.pbin+self.mbin+self.cbin)}
        # self.e2s = e2s = {i:s for i,s in enumerate(self.pbin+self.mbin+self.cbin)}
        go2slim = {}
        for w, a in [('gob','P'),('gom','F'),('goc','C')]:
            sub = df[df['aspect']==a][['GO_id', 'GOSlim_bin']]
            tmp = {}
            for g, s in sub.values:
                # tmp.setdefault(g,[]).append(s2e(s))
                tmp.setdefault(g,[]).append(s)
            go2slim[w] = tmp
        self.go2slim = go2slim

    def id2slim(self, id2gos, usecache=True):
        if usecache:
            cache = self.fname + '-id2slim.pic'
            if os.path.exists(cache):
                return load(open(cache,'rb'))
        id2slim = {}
        for w in ['gob','gom','goc']:
            i2g = id2gos[w]
            g2s = self.go2slim[w]
            id2slim[w] = {i: set([y for x in i2g[i] for y in g2s.get(x,[])]) for i in i2g}
        dump(id2slim, open(cache,'wb'))
        return id2slim
        
class MGIGOAssociation(object):
    FIELDS = ['database',
              'accession_id',
              'symbol',
              'not_designation',
              'go_id',
              'mgi_ref_accession_id',
              'evidence_code',
              'inferred_from',
              'ontology', # P=Biological Process, F=Molecular Function,C=Cellular Component    
              'name',
              'synonyms',
              'type', # gene, transcript, protein 
              'taxon',
              'modification_date',
              'assigned_by',
              'annotation_ext',
              'gene_product']

    def __init__(self, fname='gene_association.txt'):
        self.fname = fname
        fobj = open(fname)
        cnt = 0
        for line in fobj:
            cnt += 1
            if line[0]!='!':
                break
        self.df = PD.read_table(fname, names=self.FIELDS, index_col=['go_id','symbol'], comment='!', skiprows=cnt-1)

    def go2symbol(self, goid):
        return list(set(self.df.xs(goid, level='go_id').index))

    def gotree2symbol(self, goid, go):
        gos = go.get_descendants(goid)
        tmp = set()
        for g in gos:
            if g in self.df.index:
                tmp = tmp.union(set(self.df.xs(g, level='go_id').index))
        return list(tmp)

    def symbol2go(self, symbol):
        return list(set(self.df.xs(symbol, level='symbol').index))

    def id2gos(self):
        def make_id2go(ontology):
            gos = self.df[self.df['ontology']==ontology]
            symbols = list(set([x[1] for x in gos.index]))
            return {s:list(set(gos.xs(s, level='symbol').index)) for s in symbols}
        cache = self.fname+'-id2gos.pic'
        if os.path.exists(cache):
            return load(open(cache,'rb'))
        id2gos = {}
        for w,c in zip(['gob','gom','goc'],['P','F','C']):
            print('making id => %s ...' % w)
            id2gos[w] = make_id2go(c)
        dump(id2gos, open(cache,'wb'))
        return id2gos

class GO(object):
    roots = dict(gom='GO:0003674', 
                 goc='GO:0005575', 
                gob='GO:0008150')
    __no_ref_to_other_class__ = True
    
    def __init__(self, obo='gene_ontology.1_0.obo'):
        """ 
        @param obo: path to GO obo file (ver1.0)
        
        """        
        self.obo = obo
        self.dir = os.path.dirname(obo)
        pname = '%s.pic' % obo
        if os.path.exists(pname):
            print( "loading GO from pickled file...")
            dic = load(open(pname,'rb'))
            self.__dict__.update(dic)
        elif not os.path.exists(obo):
            print( "specified GO file %s does not exists, GO object not initialized" % (obo,))
            return
        else:
            print( "parsing GO file...")
            txt = open(obo,'r').read()
            #recs = re.split('^\n[\w+\]\s*\n',txt) # stanza section sep
            item_re = re.compile('\n(\w+):')
            recs = txt.split('\n\n')
            recs = [x for x in recs if '[Term]\nid:' in x]
            num = len(recs)
            id2idx = {}
            ilist = [None]*num
            rlist = [None]*num
            plist = [None]*num
            nlist = [None]*num
            dlist = ['']*num
            clist = [[] for x in range(num)] # need to instantiate independent ones
            for idx, rec in enumerate(recs):
                items = item_re.split(rec)
                is_a = []
                alt_id = []
                for k in range(1,len(items),2):
                    key = items[k].strip()
                    val = items[k+1].strip()
                    if key == 'id':
                        id = val
                    if key == 'name':
                        name = val
                    if key == 'def':
                        defi = val
                    if key == 'is_a':
                        is_a.append(val.split('!')[0].strip())
                    if key == 'alt_id':
                        alt_id.append(val)
                    if key == 'relationship':
                        key, val = val.split()[:2]
                        if key=='part_of':
                            is_a.append(val)
                id2idx[id] = idx
                for x in alt_id:
                    id2idx[x] = idx
                ilist[idx] = id
                rlist[idx] = rec
                nlist[idx] = name
                plist[idx] = is_a
                dlist[idx] = defi
            print( "  finding reverse relation...")
            for idx, parents in enumerate(plist):
                for p in parents:
                    clist[id2idx[p]].append(idx)
                plist[idx] = [id2idx[x] for x in parents]
            print( "  saving to pickled format...")
            self.id2idx = id2idx
            self.ids = ilist
            self.parents = plist
            self.children = clist
            self.names = nlist
            self.terms = rlist
            self.defs = dlist
            self.precalc_descendants()
            dump(self.__dict__, open(pname,'wb'), 2)
            
    def __getstate__(self):
        return self.obo
    
    def __setstate__(self, val):
        self.__init__(val)

    def find(self, regex, fld='terms'):
        match = re.compile(regex).search
        if fld not in ['names','terms', 'ids']:
            raise KeyError
        target = getattr(self, fld)
        return [id for id, val in zip(self.ids, target) if match(val)]

    def get_parents(self, id):
        ids = self.ids
        idx = self.id2idx[id]
        return [ids[x] for x in self.parents[idx]]

    def get_children(self, id):
        ids = self.ids
        idx = self.id2idx[id]
        return [ids[x] for x in self.children[idx]]
    
    def get_ancestors(self, id):
        tmp = set()
        for x in self.get_parents(id):
            tmp.add(x)
            tmp = tmp.union(set(self.get_ancestors(x)))
        return list(tmp)

    def get_descendants(self, id):
        tmp = {id}
        for x in self.get_children(id):
            tmp = tmp.union(set(self.get_descendants(x)))
        return list(tmp)
    
    def _descendants(self, idx):
        tmp = {idx}
        children = self.children
        for x in children[idx]:
            tmp = tmp.union(set(self._descendants(x)))
        return list(tmp)
        
    def precalc_descendants(self):
        self.descendants = [self._descendants(x) for x in range(len(self.ids))]
                
    def get_name(self, id):
        return self.names[self.id2idx[id]]

    def get_term(self, id):
        return self.terms[self.id2idx[id]]

    def get_def(self, id):
        return self.defs[self.id2idx[id]].replace('\\','').replace('"','')

    def get_roots(self):
        tmp = [i for i,x,y in izip(range(len(self.parents)),self.parents,self.children) if (len(x)==0 and len(y)>0)]
        return [self.ids[x] for x in tmp]
    
    def get_singletons(self):
        tmp = [i for i,x,y in izip(range(len(self.parents)),self.parents,self.children) if (len(x)==0 and len(y)==0)]
        return [self.ids[x] for x in tmp]
    
    
    def get_root_id(self, which):
        return self.roots[which]
    
    def get_minlevels(self):
        x2l = getattr(self, 'minlevels', None)
        if x2l:
            return x2l
        x2c = self.children
        i2x = self.id2idx
        x2l = [-1]*len(x2c) # minimum level
        def _calc(node, level):
            if x2l[node]>0: # if already visited through different parent(s)
                minlevel = min(x2l[node], level) # choose smaller
            else:
                minlevel = level # first time
            x2l[node] = minlevel
            for x in x2c[node]:
                _calc(x, minlevel+1)
        for root in self.roots.values():
            _calc(i2x[root], 0)
        self.minlevels = x2l
        return x2l
        
    def get_maxlevels(self):
        x2l = getattr(self, 'maxlevels', None)
        if x2l:
            return x2l
        x2c = self.children
        i2x = self.id2idx
        x2l = [-1]*len(x2c) # maximum level
        def _calc(node, level):
            if x2l[node]>0: # if already visited through different parent(s)
                maxlevel = max(x2l[node], level) # choose bigger
            else:
                maxlevel = level # first time
            x2l[node] = maxlevel
            for x in x2c[node]:
                _calc(x, maxlevel+1)
        for root in self.roots.values():
            _calc(i2x[root], 0)
        self.maxlevels = x2l
        return x2l

        
    def get_spanning_tree_min(self, idx, p2c=True):
        """ 
        Creat a spanning tree by choosing a unique parent by its min level.
        A child will be assigned to a parent corresponding to the minimum level.
        If there are multiple parents with a same minimum level, then intrinsic
        ordering of go.parents (which is determined by geneontology.obo term ordering) 
        will be used.
        
        Returns a tree structure expressed in a dict of id => parent id.
        idx is used for id rather than GOid ("GO:xxxxxxxx").
        To get GO id: self.ids[idx]
        To get GO name: self.names[idx]
        """
        rslt = {}
        x2l = self.get_minlevels()
        minlevel = x2l[idx]
        # choose parent
        candidate = [x for x in self.parents[idx] if x2l[x]==(minlevel-1)]
        if len(candidate)>0: # at least there should be one if not root or isolated node
            rslt[idx] = candidate[-1] # if multiple choose the left most, this should
        else:
            rslt[idx] = None # None indicates root 
        # always the same, if using the same obo file
        # now get idx=> parent map from children and concatenate them
        for c in self.children[idx]:
            rslt.update(self.get_spanning_tree_min(c, False))
        
        if not p2c:
            return rslt
        rslt2 = {}
        for idx, pidx in rslt.iteritems():
            rslt2.setdefault(pidx,[]).append(idx)
        return rslt2
    
    def get_spanning_tree_max(self, idx, p2c=True):
        rslt = {}
        x2l = self.get_maxlevels()
        maxlevel = x2l[idx]
        # choose parent
        candidate = [x for x in self.parents[idx] if x2l[x]==(maxlevel-1)]
        if len(candidate)>0: # at least there should be one if not root or isolated node
            rslt[idx] = candidate[-1] # if multiple choose the left most, this should
        else:
            rslt[idx] = None # None indicates root 
        # always the same, if using the same obo file
        # now get idx=> parent map from children and concatenate them
        for c in self.children[idx]:
            rslt.update(self.get_spanning_tree_max(c, False))
        
        if not p2c:
            return rslt
        rslt2 = {}
        for idx, pidx in rslt.iteritems():
            rslt2.setdefault(pidx,[]).append(idx)
        return rslt2
        
        
        
