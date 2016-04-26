"""

.. module:: graph
    :synopsis: Graph related stuffs

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
# system imports
from collections import Counter
import subprocess
import multiprocessing
try:
    import cPickle as pickle
except:
    import pickle
import gzip
import os
import time
from functools import reduce
from itertools import repeat
from operator import iadd
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

# 3rd party imports
import pandas as PD
import numpy as N

# library imports
from jgem import utils as UT
from jgem import bedtools as BT


# Graph #######################################################################

## MEGraph version 2 (no strand: start/end based)
class MEGraph2(object):
    
    def __init__(self, sj, me, depth=500):
        if '_id' not in sj.columns:
            UT.set_ids(sj)
        if '_id' not in me.columns:
            UT.set_ids(me)
        if ('st_id' not in sj.columns) or ('st_id' not in me.columns):
            UT.set_pos_info(sj,me) 
        self.sj = sj
        self.me = me
        self.depth=depth
        
        # prepare joined table  
        # exon[ed_id_st,e_id_st,st_id]=>junction(st_id,_id,ed_id)=>exon[ed_id,e_id_ed,st_id_ed]
        metbl = me[['ed_id','_id','st_id']]
        metbl_a = metbl[metbl['ed_id']!=-1].rename(columns={'st_id':'st_id_ed','_id':'e_id_ed'})
        metbl_d = metbl[metbl['st_id']!=-1].rename(columns={'ed_id':'ed_id_st','_id':'e_id_st'})
        sjtbl = sj[['st_id','_id','ed_id']]
        # join on donor
        j1 = PD.merge(metbl_d, sjtbl, how='outer', on='st_id', sort=False)
        j2 = PD.merge(j1, metbl_a, how='outer', on='ed_id', sort=False)              
        # remove dangling exons, junctions
        idxnd = j2['e_id_ed'].notnull()&j2['e_id_st'].notnull()
        if N.sum(idxnd)<len(idxnd):
            j2nd = j2[idxnd].copy()
        else:
            j2nd = j2
        
        # groupby exon id
        j2nd['e_id_st'] = j2nd['e_id_st'].astype(int)
        j2nd['e_id_ed'] = j2nd['e_id_ed'].astype(int)
        self.r = j2nd.groupby('e_id_ed')['e_id_st']
        self.l = j2nd.groupby('e_id_st')['e_id_ed']
        self.j2 = j2
        self.j2nd = j2nd
        self.gr = j2.groupby('ed_id')['e_id_ed'] # groupby junction right
        self.gl = j2.groupby('st_id')['e_id_st'] # groupby junction left
        self.exons = me
        
    def consistent(self):
        # return consistent subset of me, sj (i.e. non-dangling)
        j2 = self.j2
        _stids = j2[j2['e_id_ed'].isnull()]['st_id'].values
        _edids = j2[j2['e_id_st'].isnull()]['ed_id'].values
        me0 = me[(~me['ed_id'].isin(_edid))|(~me['st_id'].isin(_stid))]
        sj0 = sj[(~sj['ed_id'].isin(_edid))|(~sj['st_id'].isin(_stid))]
        return sj0, me0
        
    def connected(self, eid, exx=None, depth=0):
        #recursive version : reaches max limit with 'PR26_@myofiber_m53_1383'
        if depth>self.depth:
            LOG.debug('depth {0}: eid={1}'.format(depth,eid))
        if exx is None:
            exx = set()
        exx.add(eid)
        for e in self.ex_ex(eid):
            if e not in exx:
                exx = self.connected(e, exx, depth+1)
        return exx
        
    def allcomponents(self): # ~45sec (sid1624)
        me = self.me
        #me['_visited'] = False
        #mei = me.set_index('_id')
        #visited = mei['_visited']
        # above method of recording visited node is >x300 slower
        visited = set()
        genes = []
        tot = len(me)
        for i,eid in enumerate(me['_id'].values):
            # if i % 50000 == 0:
            #     LOG.debug('{0}/{1}...'.format(i,tot))
            #if not visited.ix[eid]:
            if eid not in visited:
                exx = self.connected(eid)
                genes.append(exx)
                #visited.loc[exx] = True
                visited.update(exx)
        return genes
        
    def ex_ex(self, eid):
        # return exons connected to eid
        return self.ex_l_ex(eid)+self.ex_r_ex(eid)
        
    def ex_l_ex(self, eid):
        # connected through donor (down-stream)
        try:
            return list(self.l.get_group(eid).values)
        except:
            return []
        
    def ex_r_ex(self, eid):
        # connected through acceptor (up-stream)
        try:
            return list(self.r.get_group(eid).values)
        except:
            return []
        
    def connected_nr(self, eid):
        # non recursive, slow? still not work for PR26 seems to be in infinite loop
        to_visit = [eid]
        exx = set()
        depth=0
        flag = False
        while(len(to_visit)>0):
            depth+=1
            c = to_visit.pop(0)
            if depth>self.depth:
                flag = True
            exx.add(c)
            for e in self.ex_ex(c):
                if (e not in exx) and (e not in to_visit):
                    to_visit.append(e)
        if flag:
            LOG.debug('eid={1} last visit = {2}: depth {0}'.format(depth,eid,c)) 
        return exx
        
    def allcomponents_nr(self): # ~44sec (sid1624)
        me = self.me
        visited = set()
        genes = []
        tot = len(me)
        for i,eid in enumerate(me['_id'].values):
            # if i % 50000 == 0:
            #     LOG.debug('{0}/{1}...'.format(i,tot))
            if eid not in visited:
                exx = self.connected_nr(eid)
                genes.append(exx)
                visited.update(exx)
        return genes    
        
    def sj_leftex(self,st_id,flds=['st','name']):#sjrec):
        try:
            left = set(self.gl.get_group(st_id).values)
            tmp = self.me.ix[left]
            return list(tmp[tmp['_id'].notnull()][flds].values)
        except:
            return []        
    def sj_rightex(self,ed_id,flds=['ed','name']):#sjrec):
        try:
            right = set(self.gr.get_group(ed_id).values)
            tmp = self.me.ix[right]
            return list(tmp[tmp['_id'].notnull()][flds].values)
        except:
            return []

    def sj_ex(self,sjrec,flds=['strand']):
        return self.sj_leftex(sjrec['st_id'],flds)+self.sj_rightex(sjrec['ed_id'],flds)
            
## MEGraph version 3 (stranded: acceptor/donor based)
class MEGraph3(object):
    """ acceptor/donor version, has problem if there's unstranded data """

    def __init__(self, sj, me, depth=500, maxcnt=990):
        UT.set_info(sj,me)
        self.sj = sj
        self.me = me
        self.depth=depth
        self.maxcnt=maxcnt
        
        # prepare joined table
        metbl = me[['a_id','_id','d_id']]
        metbl_a = metbl[metbl['a_id']!=-1].rename(columns={'d_id':'d_id_a','_id':'e_id_a'})
        metbl_d = metbl[metbl['d_id']!=-1].rename(columns={'a_id':'a_id_d','_id':'e_id_d'})
        sjtbl = sj[['d_id','_id','a_id']]
        # join on donor
        j1 = PD.merge(metbl_d, sjtbl, how='outer', on='d_id', sort=False)
        j2 = PD.merge(j1, metbl_a, how='outer', on='a_id', sort=False)              
        # remove dangling exons, junctions
        j2nd = j2[j2['e_id_a'].notnull()&j2['e_id_d'].notnull()].copy()
        
        # groupby exon id
        j2nd['e_id_d'] = j2nd['e_id_d'].astype(int)
        j2nd['e_id_a'] = j2nd['e_id_a'].astype(int)
        self.a = j2nd.groupby('e_id_a')['e_id_d']
        self.d = j2nd.groupby('e_id_d')['e_id_a']
        self.j2 = j2
        self.j2nd = j2nd
        self.ga = j2.groupby('a_id') # groupby acceptor
        self.gd = j2.groupby('d_id') # groupby donor
        self.exons = me
        
    def consistent(self):
        # return consistent subset of me, sj (i.e. non-dangling)
        j2 = self.j2
        _dids = j2[j2['e_id_a'].isnull()]['d_id'].values
        _aids = j2[j2['e_id_d'].isnull()]['a_id'].values
        me0 = me[(~me['a_id'].isin(_aid))|(~me['d_id'].isin(_did))]
        sj0 = sj[(~sj['a_id'].isin(_aid))|(~sj['d_id'].isin(_did))]
        return sj0, me0
        
    def _connected(self, eid, depth=0):
        self._cnt +=1
        #recursive version : reaches max limit with 'PR26_@myofiber_m53_1383'
        if (depth>self.depth):
            if self._pcnt<3:
                if self._pcnt==0:
                    LOG.debug('-'*10)
                LOG.debug('{2} depth {0}: eid={1}'.format(depth,eid,self._curreid))
                self._pcnt +=1
        #if exx is None:
        #    exx = set()
        self._exx.add(eid)
        if self._cnt>self.maxcnt:
            LOG.debug('cnt {0} > {2} aborting eid:{1}:'.format(depth,eid,self.maxcnt))
            return #exx
        for e in self.ex_ex(eid):
            if e not in self._exx:
                self._connected(e, depth+1)
        return 
    def connected_nr(self, eid):
        # non recursive, slow? still not work for PR26 seems to be in infinite loop
        to_visit = [eid]
        exx = set()
        depth=0
        flag = False
        while(len(to_visit)>0):
            depth+=1
            c = to_visit.pop(0)
            if depth>self.depth:
                flag = True
            exx.add(c)
            for e in self.ex_ex(c):
                if (e not in exx) and (e not in to_visit):
                    to_visit.append(e)
        if flag:
            LOG.debug('eid={1} last visit = {2}: depth {0}'.format(depth,eid,c))         
        return exx

    def allcomponents(self): # ~45sec (sid1624)
        me = self.me
        self.visited = visited = set()
        self.genes = genes = []
        tot = len(me)
        for i,eid in enumerate(me['_id'].values):
            # if i % 5000 == 0:
            #     LOG.debug('{0}/{1}...'.format(i,tot))
            if eid not in visited:
                visited.add(eid)
                self._pcnt=0
                self._cnt=0
                self._exx=set()
                self._curreid = eid
                self._connected(eid)
                genes.append(self._exx)
                visited.update(self._exx)
        return genes

    def allcomponents_nr(self): # ~44sec (sid1624)
        me = self.me
        self.visited =visited = set()
        self.genes = genes = []
        tot = len(me)
        for i,eid in enumerate(me['_id'].values):
            # if i % 50000 == 0:
            #     LOG.debug('{0}/{1}...'.format(i,tot))
            if eid not in visited:
                exx = self.connected_nr(eid)
                genes.append(exx)
                visited.update(exx)
        return genes    
                
    def connected(self, eid, depth=0):
        self._pcnt=0
        self._cnt=0
        self._exx=set()
        self._curreid = eid
        self._connected(eid)
        return self._exx
        
    def ex_ex(self, eid):
        # return exons connected to eid
        return self.ex_d_ex(eid)+self.ex_a_ex(eid)
        
    def ex_d_ex(self, eid):
        # connected through donor (down-stream)
        try:
            return list(self.d.get_group(eid).values)
        except:
            return []
        
    def ex_a_ex(self, eid):
        # connected through acceptor (up-stream)
        try:
            return list(self.a.get_group(eid).values)
        except:
            return []
        

    def sj_leftex(self,aid,did,strand): #sjrec):
        #aid,did,strand = sjrec[['a_id','d_id','strand']]
        if strand=='+':
            left = self.gd.get_group(did)
            return self.me.ix[set(left['e_id_d'].values)]
        else:
            left = self.ga.get_group(aid)
            return self.me.ix[set(left['e_id_a'].values)]
        
    def sj_rightex(self,aid,did,strand): #sjrec):
        #aid,did,strand = sjrec[['a_id','d_id','strand']]
        if strand=='+':
            right = self.ga.get_group(aid)
            return self.me.ix[set(right['e_id_a'].values)]
        else:
            right = self.gd.get_group(did)
            return self.me.ix[set(right['e_id_d'].values)]

## MEGraph version 4 (stranded, overlap of exon also counts as connection)
class MEGraph4(MEGraph3):

    def __init__(self, sj, me, filepre, depth=500, maxcnt=10000):
        MEGraph3.__init__(self, sj, me, depth, maxcnt)
        self.pre = filepre
        a = filepre+'ex1.txt.gz'
        b = filepre+'ex2.txt.gz'
        c = filepre+'ov.txt.gz'
        # calculate exon overlap to self 
        cols0 = ['chr','st','ed','strand','_id']
        a = UT.save_tsv_nidx_nhead(me[cols0], a)
        b = UT.save_tsv_nidx_nhead(me[cols0], b)
        c = BT.bedtoolintersect(a,b,c,wao=True)
        cols1 = cols0+['b_'+x for x in cols0]+['ovl']
        self.ov = ov = UT.read_pandas(c, names=cols1)
        # select same strand overlap to non-self
        self.ov1 = ov1 = ov[(ov['_id']!=ov['b__id'])&(ov['strand']==ov['b_strand'])]
        # make connected dictionary _id => [b__id's]
        tmp = ov1.groupby('_id')['b__id'].apply(lambda x: list(x)).reset_index()
        if 'index' in tmp.columns:
            tmp['_id'] = tmp['index']
        #LOG.debug('graph.MEGraph4.__init__: tmp.columns={0}, len(tmp)={1}'.format(tmp.columns, len(tmp))) 
        self.eoe = dict(UT.izipcols(tmp, ['_id','b__id']))
        # cleanup
        os.unlink(a)
        os.unlink(b)
        os.unlink(c)

    def ex_ex(self, eid):
        # return exons connected to eid
        return self.ex_d_ex(eid)+self.ex_a_ex(eid)+self.eoe.get(eid,[])
        
def _worker2(s,e,c):
    mg = MEGraph2(s,e)
    genes = mg.allcomponents()
    #LOG.debug("finished {0}...".format(c))
    return genes
    
def _worker3(s,e,c):
    mg = MEGraph3(s,e)
    #genes = mg.allcomponents()
    genes = mg.allcomponents_nr()
    #LOG.debug("finished {0}...".format(c))
    #cPickle.dump(genes, open('genes-{0}.pic'.format(c),'w'))
    return genes

def _worker4(s,e,c,fp):
    fpc = fp+c+'.'
    mg = MEGraph4(s,e,fpc)
    #genes = mg.allcomponents()
    genes = mg.allcomponents_nr()
    #LOG.debug("finished {0}...".format(c))
    #cPickle.dump(genes, open('genes-{0}.pic'.format(c),'w'))
    return genes    
        
def worker(args):
    func, arg = args
    return func(*arg)

def mcore_allcomponents(sj, me, np=4, version=3, depth=500, maxcnt=10000, chroms=None, filepre=''): # ~31sec (np=1)
    # np number of processes
    # spawn processes to process chrom-wise
    # common id
    if '_id' not in sj.columns:
        UT.set_ids(sj)
    if '_id' not in me.columns:
        UT.set_ids(me)
    if version==2:
        if ('st_id' not in sj.columns) or ('st_id' not in me.columns):
            UT.set_pos_info(sj,me) 
    else:
        if ('a_pos' not in sj.columns) or ('a_pos' not in me.columns):
            UT.set_ad_info(sj, me) 
    if chroms is None:
        chroms = sorted(me['chr'].unique())
    if version==4:
        data = [(sj[sj['chr']==c], me[me['chr']==c], c, filepre) for c in chroms]
    else:
        data = [(sj[sj['chr']==c], me[me['chr']==c], c) for c in chroms]
    data = [x for x in data if (len(x[0])>0 and len(x[1])>0)]
    if np==1:
        rslts = []
        for d in data:
            if version==3:
                s,e,c = d
                LOG.debug('connected component: processing {0}...'.format(c))
                mg = MEGraph3(s,e,depth=depth, maxcnt=maxcnt)
            elif version==4:
                s,e,c,filepre = d
                LOG.debug('connected component: processing {0}...'.format(c))
                mg = MEGraph4(s,e,filepre,depth=depth, maxcnt=maxcnt)
            else:
                s,e,c = d
                LOG.debug('connected component: processing {0}...'.format(c))
                mg = MEGraph2(s,e,depth=depth)
            tmp = mg.allcomponents_nr()
            rslts.append(tmp)
            LOG.debug("finished {0}...".format(c))
            #cPickle.dump(tmp, open('genes-{0}.pic'.format(c),'w'))
        genes = [x for y in rslts for x in y]
    else:
        p = multiprocessing.Pool(np)
        if version==3:
            a = zip(repeat(_worker3), data)
        elif version==4:
            a = zip(repeat(_worker4), data)
        else:
            a = zip(repeat(_worker2), data)
        genes = reduce(iadd, p.map(worker, a))
        p.close()
        #p.join()
    return genes

def mcore_allcomponents2(sj, me, np=4, depth=500, maxcnt=10000, chroms=None):
    return mcore_allcomponents(sj, me, np, 2, depth, maxcnt, chroms)

def mcore_allcomponents3(sj, me, np=4, depth=500, maxcnt=10000, chroms=None):
    return mcore_allcomponents(sj, me, np, 3, depth, maxcnt, chroms)

def mcore_allcomponents4(sj, me, filepre, np=4, depth=500, maxcnt=10000, chroms=None):
    return mcore_allcomponents(sj, me, np, 4, depth, maxcnt, chroms, filepre)


# GENE FINDER ###############################################################

def find_genes3(sj, ae, cachename=None, np=1, override=False, depth=140):
    """ 
    Adds _gidx column to ae
    Connection: by junctions

    Returns genes [set([_id,..]), ...]
    """
    if '_id' not in ae.columns:
        LOG.debug('setting ex _id...')
        UT.set_ids(ae)
    if '_id' not in sj.columns:
        LOG.debug('setting sj _id...')
        UT.set_ids(sj)
    if 'cat' not in ae.columns:
        UT.set_exon_category(sj,ae)

    ### FIND GENES
    if cachename and os.path.exists(cachename) and not override:
        LOG.debug('loading cached genes (connected components)...')
        genes = pickle.load(open(cachename, 'rb'))
    else:
        LOG.debug('finding genes (connected components)...')
        _sttime = time.time()
        me, se = UT.mese(ae)
        genes = mcore_allcomponents3(sj, me, np, depth=depth)
        # genes = [set([_id's]),...]
        # SE genes
        genes += [set([x]) for x in se['_id']]
        if cachename:
            pickle.dump(genes, open(cachename,'wb'))
        LOG.debug(' time: {0:.3f}s'.format(time.time()-_sttime))
    
    ### WRITE EXONS W/ GENE number
    LOG.debug('assigning gidx...')
    _sttime = time.time()
    # ae['_gidx'] = 0
    # ae.index = ae['_id']
    # for i, ids in enumerate(genes):
    #     ae.ix[ids, '_gidx'] = i+1
    i2g = {}
    for i, ids in enumerate(genes):
        for x in ids:
            i2g[x] = i+1
    ae['_gidx'] = [i2g[x] for x in ae['_id']]

    ## set sj _gidx, use acceptor=>_gidx map (exon a_id, sj a_id)
    a2g = dict(UT.izipcols(ae, ['a_id','_gidx']))
    d2g = dict(UT.izipcols(ae, ['d_id','_gidx']))
    sj['_gidx'] = [a2g.get(x,d2g.get(y,-1)) for x,y in UT.izipcols(sj,['a_id','d_id'])]

    
    # This shouldn't happen
    nidx = ae['_gidx']==0
    if N.sum(nidx)>0:
        LOG.warning('###### WARNING!!!!!! exons with no gene assignment:{0}'.format(N.sum(nidx)))
        #ae.loc[nidx, '_gidx'] = N.arange(len(ae),len(ae)+N.sum(nidx))        
        
    return genes

def find_genes4(sj, ae, filepre, cachename=None, np=1, override=False, depth=500, separatese=True):
    """ 
    Adds _gidx column to ae
    Connection: 1) by junctions, 2) by overlap in the same strand

    Returns genes [set([_id,..]), ...]
    """
    if '_id' not in ae.columns:
        LOG.info('setting ex _id...')
        UT.set_ids(ae)
    if '_id' not in sj.columns:
        LOG.info('setting sj _id...')
        UT.set_ids(sj)
    if 'cat' not in ae.columns:
        UT.set_exon_category(sj,ae)
    if 'a_id' not in ae.columns:
        UT.set_ad_info(sj,ae)

    ### FIND GENES
    if cachename and os.path.exists(cachename) and not override:
        LOG.info('loading cached genes (connected components)...')
        genes = pickle.load(open(cachename, 'rb'))
    else:
        LOG.info('finding genes (connected components)...')
        _sttime = time.time()
        if separatese:
            me, se = UT.mese(ae)
            genes = mcore_allcomponents4(sj, me, filepre, np, depth=depth)
            # SE genes
            genes += [set([x]) for x in se['_id']]
        else:
            genes = mcore_allcomponents4(sj, ae, filepre, np, depth=depth)
        # version 4 graph: uses overlaps in addition to junctions to connect
        # genes = [set([_id's]),...]
        if cachename:
            pickle.dump(genes, open(cachename,'wb'))
        LOG.info(' time: {0:.3f}s'.format(time.time()-_sttime))
    
    ### WRITE EXONS W/ GENE number
    LOG.info('assigning gidx...')
    _sttime = time.time()
    i2g = {} # eid => _gidx
    i2gn = {} # eidt => gname
    g2gn = {}
    i2s = dict(UT.izipcols(ae, ['_id','strand'])) # eid => strand
    #i2c = dict(UT.izipcols(ae, ['_id','cat'])) # eid => category
    s2n = {'+':'P','-':'N','.':''}
    c2n = {'s':'S','i':'G','5':'G','3':'G'}
    for i, ids in enumerate(genes):
        gid = i+1
        strand = s2n[i2s[list(ids)[0]]]
        cat = 'S' if len(ids)==1 else 'G'
        if strand=='N': # negative strand
            gid = -gid
        gname = '{0}{1}{2}'.format(strand,cat,abs(gid))
        g2gn[gid] = gname
        for x in ids:
            i2g[x] = gid
            i2gn[x] = gname

    ae['_gidx'] = [i2g[x] for x in ae['_id']]
    ae['gname'] = [i2gn[x] for x in ae['_id']]

    ## set sj _gidx, use acceptor=>_gidx map (exon a_id, sj a_id)
    a2g = dict(UT.izipcols(ae, ['a_id','_gidx']))
    d2g = dict(UT.izipcols(ae, ['d_id','_gidx']))
    sj['_gidx'] = [a2g.get(x,d2g.get(y,0)) for x,y in UT.izipcols(sj,['a_id','d_id'])]
    sj['gname'] = [g2gn.get(x,'') for x in sj['_gidx']]
    
    # This shouldn't happen
    nidx = ae['_gidx']==0
    if N.sum(nidx)>0:
        LOG.warning('###### WARNING!!!!!! exons with no gene assignment:{0}'.format(N.sum(nidx)))
        #ae.loc[nidx, '_gidx'] = N.arange(len(ae),len(ae)+N.sum(nidx))        
        
    return genes


