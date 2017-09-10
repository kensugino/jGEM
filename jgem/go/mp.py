# Mammalian Phenotype
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



class MP(object):
    roots = dict(mp='MP:0000001')
    __no_ref_to_other_class__ = True
    
    def __init__(self, obo='MPheno_OBO.ontology'):
        """ 
        @param obo: path to MP obo file (ver1.2)
        
        """        
        self.obo = obo
        self.dir = os.path.dirname(obo)
        pname = '%s.pic' % obo
        if os.path.exists(pname):
            print( "loading OBO from pickled file...")
            dic = load(open(pname,'rb'))
            self.__dict__.update(dic)
        elif not os.path.exists(obo):
            print( "specified OBO file %s does not exists, OBO object not initialized" % (obo,))
            return
        else:
            print( "parsing OBO file...")
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