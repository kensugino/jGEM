"""

.. module:: assembler2
    :synopsis: assemble genes from RNASeq data (reads and junction coverages (bigwig) and junction paths)
    jGEM version 2 assembler

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
# system imports
import subprocess
import multiprocessing
import gzip
import os
import time
import shutil
from functools import reduce
from operator import iadd, iand
from collections import Counter
from itertools import repeat
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)
import json

# 3rd party imports
import pandas as PD
import numpy as N
import matplotlib.pyplot as P


# LocalAssembler imports
from collections import Counter
from matplotlib.collections import BrokenBarHCollection
from functools import partial, reduce
from operator import iadd
import bisect
from scipy.optimize import nnls


# library imports
from jgem import utils as UT
from jgem import bigwig as BW
from jgem import bedtools as BT
from jgem import gtfgffbed as GGB
from jgem import taskqueue as TQ

####### BigWigs ########################################################################
     

class SjExBigWigs(object):
    
    def __init__(self, bwpre, sjbwpre=None, mixunstranded=True):
        if sjbwpre is None:
            sjbwpre = bwpre
        if type(bwpre)!=type([]):
            bwpre = [bwpre]
        if type(sjbwpre)!=type([]):
            sjbwpre = [sjbwpre]
        self.bwpre = bwpre
        self.sjbwpre = sjbwpre
        self.mixunstranded = mixunstranded
        S2S = {'+':'.p','-':'.n','.':'.u','r+':'.rp','r-':'.rn','r.':'.ru'}
        bwp = {
            'ex': {s:[b+'.ex{0}.bw'.format(S2S[s]) for b in bwpre] for s in S2S},
            'sj': {s:[b+'.sj{0}.bw'.format(S2S[s]) for b in sjbwpre] for s in S2S},
        }
        self.bwpaths = bwpaths = {'ex':{},'sj':{}}
        if mixunstranded:
            bwpaths['ex']['+'] = {'p':bwp['ex']['+']+bwp['ex']['.'],}
            bwpaths['ex']['-'] = {'p':bwp['ex']['-']+bwp['ex']['.'],}
            bwpaths['ex']['.'] = {'p':bwp['ex']['.'],}        
            bwpaths['ex']['a'] = {'p':bwp['ex']['+']+bwp['ex']['-']+bwp['ex']['.'],}
            bwpaths['sj']['+'] = {'p':bwp['sj']['+']+bwp['sj']['.'],}
            bwpaths['sj']['-'] = {'p':bwp['sj']['-']+bwp['sj']['.'],}
            bwpaths['sj']['.'] = {'p':bwp['sj']['.'],}
            bwpaths['sj']['a'] = {'p':bwp['sj']['+']+bwp['sj']['-']+bwp['sj']['.'],}
        else:
            bwpaths['ex']['+'] = {'p':bwp['ex']['+'],}
            bwpaths['ex']['-'] = {'p':bwp['ex']['-'],}
            bwpaths['ex']['.'] = {'p':bwp['ex']['.'],}        
            bwpaths['ex']['a'] = {'p':bwp['ex']['+']+bwp['ex']['-'],}
            bwpaths['sj']['+'] = {'p':bwp['sj']['+'],}
            bwpaths['sj']['-'] = {'p':bwp['sj']['-'],}
            bwpaths['sj']['.'] = {'p':bwp['sj']['.'],}
            bwpaths['sj']['a'] = {'p':bwp['sj']['+']+bwp['sj']['-'],}

        # bwpaths['ex']['+'] = {'p':[bwp['ex']['+'],bwp['ex']['.']],
        #                       'n':[bwp['ex']['r+'],bwp['ex']['r.']]}
        # bwpaths['ex']['-'] = {'p':[bwp['ex']['-'],bwp['ex']['.']],
        #                       'n':[bwp['ex']['r-'],bwp['ex']['r.']]}
        # bwpaths['ex']['.'] = {'p':[bwp['ex']['.']],
        #                       'n':[bwp['ex']['r.']]}
        # bwpaths['sj']['+'] = {'p':[bwp['sj']['+'],bwp['sj']['.']],
        #                       'n':[bwp['sj']['r+'],bwp['sj']['r.']]}
        # bwpaths['sj']['-'] = {'p':[bwp['sj']['-'],bwp['sj']['.']],
        #                       'n':[bwp['sj']['r-'],bwp['sj']['r.']]}
        # bwpaths['sj']['.'] = {'p':[bwp['sj']['.']],
        #                       'n':[bwp['sj']['r.']]}
        
        self.make_bws()
    
    def make_bws(self):
        bwp = self.bwpaths
        self.bws = bws = {}
        for k in ['ex','sj']: 
            # bws[k] = {s: BW.MultiBigWigs(plus=bwp[k][s]['p'],
            #                          minus=bwp[k][s]['n']) for s in ['+','-','.']}
            bws[k] = {s: BW.MultiBigWigs(plus=bwp[k][s]['p']) for s in ['+','-','.','a']}
            for s in bws[k]:
                bws[k][s].make_bws()
        
    def __enter__(self):
        for k in ['ex','sj']:
            for s in self.bws[k]:
                self.bws[k][s].__enter__()

    def __exit__(self, exc_type, exc_value, traceback):
        for k in ['ex','sj']:
            for s in self.bws[k]:
                self.bws[k][s].__exit__(exc_type, exc_value, traceback)


####### Classifiers ####################################################################

class LogisticClassifier(object):
    
    def __init__(self, json, dstcol):
        self.json = json
        self.b0 = json['intercept'] #b0
        self.bs = N.array(json['coef']) #bs
        self.cols = json['cols']
        self.dstcol = dstcol
        
    def classify(self, df):
        X = df[self.cols].values
        d = self.b0 + N.dot(X, self.bs)
        e = N.exp(-d)
        p = 1./(1.+e)
        df[self.dstcol] = p>0.5
        

# for intergenic        
itg_p = dict(coef = N.array([ -0.40,  -4.72, 0.86]),
           intercept = 11.8638183,
           cols = ['lemax','lgap','llen'],
           th = 0.05,
           zoom = 1)
INTG = LogisticClassifier(json=itg_p, dstcol='exon')

e53_p = dict(coef = N.array([2.51, -0.77]),
             intercept= -2.7,
             cols = ['sdiff','smean'],
             sdiffth= 1,
             zoom = 1)
E53C = LogisticClassifier(json=e53_p, dstcol='e53')

e53m_p = dict(coef = N.array([1.542, -0.368]),
              intercept = -0.329,
              cols = ['sdiff','smean'],
             sdiffth= 0.5,
             zoom = 1)
E53CM = LogisticClassifier(json=e53m_p, dstcol='e53')


class EdgeFinder(object):
    
    def __init__(self, json):
        self.json=json
        a_lsin,a_lgap = json['coef']
        b0 = json['intercept']
        th = json['th']
        self.zoom = zoom = json['zoom']
        self.a_lsin = a_lsin
        self.a_lgap = a_lgap
        self.b0 = b0
        self.c0 = -b0/a_lgap
        self.c1 = -a_lsin/a_lgap
        self.th = th
        self.maxsize = json['maxsize']
        
    def find(self, sja, exa, direction):
        # sja1, exa1 : pos0=>pos+1(<), pos-1=>pos0(>)
        c0,c1,th = self.c0,self.c1,self.th
        # def _find_gap_from_idx(idx):
        #     if len(idx)==0:
        #         return
        #     i0 = idx[0] # current gap start
        #     i1 = i0 # current gap end
        #     for i2 in idx[1:]: # next point <= th
        #         if (i1+1)!=i2: # not continuous
        #             yield (i0, i1-i0+1) # position and gapsize
        #             i0 = i2 # new start
        #         i1 = i2
        #     if i0!=i1:
        #         yield (i0, i1-i0+1)
        zoom = self.zoom
        def _find_gap_from_idx(idx):
            if len(idx)==0:
                return []
            dif = idx[1:]-idx[:-1] # if continuous dif==1
            idx2 = N.nonzero(dif>1)[0] # non contiguous   
            idxst = N.array([idx[0]]+list(idx[idx2+1]))
            idxed = N.array(list(idx[idx2])+[idx[-1]])
            gsize = idxed - idxst + 1
            return zip(idxst, gsize)
        if direction=='>':
            # lsin = abs(sja[1]-sja[0])
            ein = N.mean(exa[1:11])
            lein = N.log2(zoom*ein+1)
            gapth = 2**(c0+c1*lein)-1
            #print('gapth={0:.2f}, lsin={1:.2f}'.format(gapth, lsin))
            # pos => pos0, find position where lgap > gapth
            idx = N.nonzero(exa[1:]<=th*ein)[0]
            #print(idx)
            epos = len(exa)-1 # all the way to the end
            for x in _find_gap_from_idx(idx):
                if x[1]>gapth: # found
                    epos = x[0] # start offset pos
                    break
            epos = min(epos, self.maxsize)
        else:
            # lsin = abs(sja[-1]-sja[-2])
            ein = N.mean(exa[-12:-1])
            lein = N.log2(zoom*ein+1)
            gapth = 2**(c0+c1*lein)-1
            #print('gapth={0:.2f}, lsin={1:.2f}'.format(gapth, lsin))
            # pos0 <= pos, going opposite way
            idx = N.nonzero(exa[:-1][::-1]<=th*ein)[0]
            epos = -len(exa)
            for x in _find_gap_from_idx(idx):
                if x[1]>gapth:
                    epos = -x[0]+1
                    break
            epos = max(epos, -self.maxsize)            
        return epos
        
e5_p  = dict(coef=[-0.285,-0.81], intercept=5.6, th=0, zoom=1, maxsize=3000)
EF5 = EdgeFinder(e5_p)
e3_p = dict(coef=[-0.25,-0.51], intercept=4.6, th=0, zoom=1, maxsize=25000) # -0.25, -0.5, 4.5
EF3 = EdgeFinder(e3_p) 


####### Gene Graph ###########################################

class PathNumUpperLimit(Exception):
    pass

def overlap(x,y):
    if x[-1]<y[0]:
        return False
    try:
        l = y.index(x[-1])+1
        return x[-l:]==y[:l]
    except:
        return False

def find_all_connections(paths):
    # sid, jid
    conn = {} # sid => [sid,...] overlap
    vals = paths[['sid','jid']].values
    n = len(vals)
    for i in range(n):
        sid1,jid1 = vals[i]
        for j in range(n):
            if i==j:
                continue
            sid2,jid2 = vals[j]
            if overlap(jid1, jid2):
                conn.setdefault(sid1, []).append(sid2)
    return conn

def merge_pathcode(a,b):
    # a: d1|a1,...,dk|ak,...,dn|an
    # b:           dk|ak,...,dn|an,...,dm|am
    alast = a.split(',')[-1]
    pos = b.find(alast)
    #print(a, b[pos:], b[pos+len(alast)+1:])
    if pos<0:
        raise RuntimeError('non overlapping pathcodes {0} and {1}'.format(a,b))
    c = b[pos+len(alast)+1:]
    if len(c)>0:
        return '{0},{1}'.format(a, c)
    return a

def assign_jids(sjs, strand):
    if '+' in strand:
        sjs.sort_values(['tst','ted'], inplace=True)
    else:
        sjs.sort_values(['ted','tst'], inplace=True, ascending=False)
    sjs['sid'] = N.arange(1,len(sjs)+1)

    sju0 = set()
    for n in sjs['name']:
        for x in n.split(','):
            if '|' in x:
                sju0.add(x)
    tmp = sorted([[[int(y) for y in x.split('|')],x] for x in sju0])
    if '-' in strand:
        tmp = tmp[::-1]
    sju = {}
    usj = {} # reverse map
    for i,x in enumerate(tmp):
        sju[x[1]]=i+1 # junction name:=>numbering 1=>N from 5' to 3'
        usj[i+1] = x[1] # numbering:=>junction name

    sjs['jid'] = [[sju[x] for x in y.split(',') if '|' in x] for y in sjs['name']]
    sjs['jc'] = [','.join([str(x) for x in y]) for y in sjs['jid']]
    s2jc = UT.df2dict(sjs, 'sid', 'jc')

    # jid => strand
    j2strand = {}
    for jids,strand in sjs[['jid','strand']].values:
        for j in jids:
            j2strand.setdefault(j, []).append(strand)
    j2s = {j:Counter(j2strand[j]).most_common()[0][0] for j in j2strand}

    return s2jc, sju, usj, j2s

def assign_jids1(sjs, strand):
    if '+' in strand:
        sjs.sort_values(['tst','ted'], inplace=True)
    else:
        sjs.sort_values(['ted','tst'], inplace=True)
    sjs['sid'] = N.arange(1,len(sjs)+1)

    sju0 = set()
    for n in sjs['name']:
        for x in n.split(','):
            if '|' in x:
                sju0.add(x)
    tmp = sorted([[[int(y) for y in x.split('|')],x] for x in sju0])
    if '-' in strand:
        tmp = tmp[::-1]
    sju = {}
    for i,x in enumerate(tmp):
        sju[x[1]]=i+1 # junction name:=>numbering 1=>N from 5' to 3'

    sjs['jid'] = [[sju[x] for x in y.split(',') if '|' in x] for y in sjs['name']]
    sjs['jc'] = [','.join([str(x) for x in y]) for y in sjs['jid']]

    return sju

def assign_jids2(sjs, strand, sju):
    if '+' in strand:
        sjs.sort_values(['tst','ted'], inplace=True)
    else:
        sjs.sort_values(['ted','tst'], inplace=True)
    sjs['sid'] = N.arange(1,len(sjs)+1)
    sjs['jid'] = [[sju[x] for x in y.split(',') if '|' in x] for y in sjs['name']]
    sjs['jc'] = [','.join([str(x) for x in y]) for y in sjs['jid']]

class GeneGraph(object):
    
    def __init__(self, exons, sjpaths, strand, upperpathnum=5000):
        self.exs = exons#.copy()
        self.sjs = sjpaths#.copy()
        self.strand = strand
        self.upperpathnum0 = upperpathnum
        self.upperpathnum = max(len(sjpaths), upperpathnum)

    def prep_gstree(self, sjs):
        # decompose into unitary junctions, make GSTree
        self.gs_sid2jc, self.gs_sju, self.gs_usj, self.gs_j2s = assign_jids(sjs, self.strand)
        #recs = sjs[['sid','jid']].values
        # self.gstree = gstree = GSTree(recs)
        # self.p2p = gstree.find_all_connections()
        self.p2p = find_all_connections(sjs)
        self.gs_sjs = sjs
        self.gs_sid2pc = UT.df2dict(sjs, 'sid', 'name')


    def prep_ggraph(self, sjs):
        exs = self.exs
        exs['eid'] = N.arange(1,len(exs)+1)
        strand = self.strand
        self.sid2jc, self.sju, self.usj, self.j2s = assign_jids(sjs, strand)
        if strand=='+':
            sjs['apos'] = sjs['ted']
            sjs['dpos'] = sjs['tst']
            exs['apos'] = exs['st']
            exs['dpos'] = exs['ed']
            #st,ed,tst,ted='st','ed','tst','ted'
            #sjs.sort_values(['dpos','apos'], ascending=True, inplace=True)
        else:
            sjs['apos'] = sjs['tst']
            sjs['dpos'] = sjs['ted']
            exs['apos'] = exs['ed']
            exs['dpos'] = exs['st']
            #st,ed,tst,ted='ed','st','ted','tst'
            #sjs.sort_values(['dpos','apos'], ascending=False, inplace=True)
        etbl = exs[['apos','eid','dpos']]
        stbl = sjs[['apos','sid','dpos']]
        stbl1 = stbl.rename(columns={'dpos':'dpos1','sid':'sid1'})
        stbl2 = stbl.rename(columns={'sid':'sid2','apos':'apos2'})
        j1 = PD.merge(stbl1, etbl, how='outer', on='apos', sort=False)
        j2 = PD.merge(j1, stbl2, how='outer', on='dpos', sort=False)
        self.j2 = j2
        self.s2s = j2.groupby('sid1')['sid2'].apply(lambda x: [int(y) for y in set(x) if not N.isnan(y)])    
        self.sid2pc = UT.df2dict(sjs, 'sid', 'name')
        
    # def sj2ex(self, sjrec):
    #     return self.exons[self.exon['st']==sjrec['ted']]
    # def ex2sj(self, exrec):
    #     return self.sjpaths[self.sjpaths['tst']==exrec['ed']]
    # def sj2sj(self, sid):
    #     return self.sjs.set_index('sid').ix[self.s2s.ix[sid]]
    
    def find_a_tree_sj_ovl(self, sid, visited, allus):
        # return all paths in a list of **jids** 
        # allus: path components all unstranded 
        if sid in visited:
            return visited[sid]
        children =  self.p2p.get(sid, [])
        pc = self.gs_sid2jc[sid]
        j2s = self.gs_j2s
        unstranded = all(['.' in j2s[int(x)] for x in pc.split(',')])
        if len(children)>0:
            rslt = set()
            for y in children:
                for x in self.find_a_tree_sj_ovl(y, visited, allus):
                    npc = merge_pathcode(pc, x) #'{0}.{1}'.format(pc,x)
                    rslt.add(npc)
                    if x in allus and unstranded:
                        allus[npc]=True
        else:
            rslt = set([pc])
            if unstranded:
                allus[pc] = True
        visited[sid] = rslt
        return rslt

    def find_a_tree_ex(self, sid, visited, allus):
        # return all unique paths in as a list of pathcodes
        # allus: jc => True when jc components all unstranded use j2s dic
        if sid in visited:
            return visited[sid] # cached
        children = self.s2s.ix[sid]
        pc = self.sid2jc[sid]
        j2s = self.j2s
        unstranded = all(['.' in j2s[int(x)] for x in pc.split(',')])
        if len(children)>0:
            rslt = set()
            for y in children:
                for x in self.find_a_tree_ex(y, visited, allus):
                    npc = '{0},{1}'.format(pc,x)
                    rslt.add(npc)
                    if x in allus and unstranded:
                        allus[npc]=True
        else:
            rslt = set([pc])
            if unstranded:
                allus[pc] = True
        visited[sid] = rslt
        if len(rslt)>self.upperpathnum:
            raise PathNumUpperLimit
        return rslt
    
    def find_all_paths0(self):
        self.prep_gstree(self.sjs)
        newpaths = self.find_all_paths_sj_ovl0(self.sjs)
        if newpaths is not None:
            cols = ['name','tst','ted','strand']
            paths = PD.concat([self.sjs[cols], newpaths[cols]], ignore_index=True)
            self.mergedpaths = paths = paths.sort_values(['tst','ted'])
            paths['sid'] = N.arange(1,len(paths)+1)
        else:
            paths = self.sjs
        self.prep_ggraph(paths)
        return self.find_all_paths_ex(paths)

    def find_all_paths(self):
        self.prep_ggraph(self.sjs)
        paths = self.find_all_paths_ex(self.sjs)
        if len(paths)>1:
            self.prep_gstree(paths)
            paths = self.find_all_paths_sj_ovl(paths)
        return paths

    def find_all_paths_sj_ovl(self, sjs):
        notused = sjs
        allpaths = set()
        self.gs_visited = visited = {}
        self.gs_allus = allus = {}  # pc element all unstranded
        while len(notused)>0:
            if self.strand=='+':
                rootpos = notused['tst'].min()
                rootsids = notused[notused['tst']==rootpos]['sid'].values
                tst,ted = 'tst','ted'
            else:
                rootpos = notused['ted'].max()
                rootsids = notused[notused['ted']==rootpos]['sid'].values
                tst,ted = 'ted','tst'
            for rid in rootsids:
                #print('processing rid {0}'.format(rid))
                x = self.find_a_tree_sj_ovl(int(rid), visited, allus) # set of paths
                allpaths.update(x)
                #print('  subtree #{0}'.format(len(x)))
            notused = sjs[~sjs['sid'].isin(visited)]
        
        # make for each path (tst,ted,pathcode) 
        sjsi = sjs.set_index('sid')
        usj = self.gs_usj
        recs = {}
        self.allpaths_sj = allpaths
        for jc in allpaths:
            # sids = pc2sids[name]
            # sjsub = sjsi.ix[sids]
            # sjsub = sjsi.ix[list(path)]
            # jids = sorted(set(reduce(iadd, sjsub['jid'].values, [])))
            # name = ','.join([usj[x] for x in jids])
            s = '.'+self.strand if allus.get(jc,False) else self.strand
            name = ','.join([usj[int(x)] for x in jc.split(',')])
            if name not in recs:
                # cnt = Counter(sjsub['strand'].values)
                # pos = list(sjsub['tst'])+list(sjsub['ted'])
                tmp = name.split(',')
                dpos = int(tmp[0].split('|')[0])
                apos = int(tmp[-1].split('|')[1])
                pos = [dpos, apos]
                tst = N.min(pos)
                ted = N.max(pos)
                recs[name] = dict(
                    tst = tst,
                    ted = ted,
                    strand = s, #cnt.most_common()[0][0],
                )
        if len(recs)>0:
            df = PD.DataFrame(recs).T
            df.index.name = 'name'
            df['sid'] = N.arange(len(df))
            df = df.reset_index().sort_values(['tst','ted'])
        else:
            df = None
        self.pathsdf_sj = df
        return df

   
    def _find_all_paths_ex_1(self, sjs):
        # root?
        notused = sjs
        allpaths = set()
        self.visited = visited = {}
        self.allus = allus = {}
        while len(notused)>0:
            if self.strand=='+':
                rootpos = notused['tst'].min()
                rootsids = notused[notused['tst']==rootpos]['sid'].values
                tst,ted = 'tst','ted'
            else:
                rootpos = notused['ted'].max()
                rootsids = notused[notused['ted']==rootpos]['sid'].values
                tst,ted = 'ted','tst'
            for rid in rootsids:
                #print('processing rid {0}'.format(rid))
                x = self.find_a_tree_ex(int(rid), visited, allus) # set of paths
                allpaths.update(x)
                if len(allpaths)>self.upperpathnum:
                    raise PathNumUpperLimit
            notused = sjs[~sjs['sid'].isin(visited)]
        
        # make for each path (tst,ted,pathcode) 
        #print('path finding finished, #paths={0}'.format(len(allpaths)))
        self.allpaths_ex = allpaths
        return allpaths

    def find_all_paths_ex(self, sjs):
        sjrth = 0.002
        uth = sjs['sc1'].min()
        while True:
            try:
                allpaths = self._find_all_paths_ex_1(sjs)
                if len(allpaths)>1000:
                    chrom = sjs.iloc[0]['chr']
                    st = sjs['st'].min()
                    ed = sjs['ed'].max()
                    LOG.warning('#path>1000 ({0}) at {1}:{2}-{3}'.format(len(allpaths),chrom,st,ed))
                break
            except PathNumUpperLimit:
                LOG.warning('Too many paths (>{0}). Possible repeats. Increasing stringency.'.format(self.upperpathnum))
                chrom = sjs.iloc[0]['chr']
                stmin = sjs['st'].min()
                edmax = sjs['ed'].max()
                LOG.debug('location: {0}:{1}-{2}'.format(chrom,stmin,edmax))

                uth += 0.4
                # sc2min = sjs['sc2'].min()
                mcnt = sjs['sc2']-sjs['sc1'] # multi mappers
                mth = max(0, mcnt.max()/2.)
                sjrth += 0.0005
                n0 = len(sjs)
                # ucnt threshold increases by 1, mcnt threshold decrease by 10, sjratio by 0.001
                sjs = sjs[(sjs['sc1']>uth)&(mcnt<=mth)&(sjs['sjratio']>sjrth)].copy()
                n1 = len(sjs)
                LOG.debug('#sjs:{0}=>{1}, uth:{2}, mth:{3}, sjrth:{4}'.format(n0,n1,uth,mth,sjrth))
                self.prep_ggraph(sjs)


        allus = self.allus
        sjsi = sjs.set_index('sid')
        usj = self.usj
        recs = {}
        for jc in allpaths: # path = pathcode 
            # sids = pc2sids[name]
            # sjsub = sjsi.ix[sids]
            # name = ','.join(sjsub['name'].values)
            s = '.'+self.strand if allus.get(jc,False) else self.strand
            name = ','.join([usj[int(x)] for x in jc.split(',')])
            if name not in recs:
                tmp = name.split(',')
                dpos = int(tmp[0].split('|')[0])
                apos = int(tmp[-1].split('|')[1])
                pos = [dpos, apos]
                # cnt = Counter(sjsub['strand'].values)
                # pos = list(sjsub['tst'])+list(sjsub['ted'])
                tst = N.min(pos)
                ted = N.max(pos)
                recs[name] = dict(
                    tst = tst,
                    ted = ted,
                    strand = s, #cnt.most_common()[0][0],
                )
        df = PD.DataFrame(recs).T
        df.index.name = 'name'
        self.pathsdf_ex = df = df.reset_index()
        return df

####### Colors  ########################################################################
# move to plotutil

import matplotlib.colors as C
import matplotlib.cm as CM

class Colors(object):
    
    def __init__(self, mapname, vmax, vmin=0, nl=32):
        self.mn = mapname
        self.vmin = vmin
        self.vmax = vmax
        self.d = d = 1./nl
        if mapname=='C':
            self.rgba = [(1.-x,1.,1.,1.) for x in N.arange(0,1+d,d)]
        elif mapname=='M':
            self.rgba = [(1.,1.-x,1.,1.) for x in N.arange(0,1+d,d)]
        elif mapname=='Y':
            self.rgba = [(1.,1.,1.-x,1.) for x in N.arange(0,1+d,d)]
        elif mapname=='R':
            self.rgba = [(1.,1.-x,1.-x,1.) for x in N.arange(0,1+d,d)]
        elif mapname=='G':
            self.rgba = [(1.-x,1.,1.-x,1.) for x in N.arange(0,1+d,d)]
        elif mapname=='B':
            self.rgba = [(1.-x,1.-x,1.,1.) for x in N.arange(0,1+d,d)]
        else:
            cm = P.get_cmap(mapname)
            cnorm = C.Normalize(vmin=0,vmax=1.)
            self.sm = sm = CM.ScalarMappable(norm=cnorm,cmap=cm)
            self.rgba = [sm.to_rgba(x) for x in N.arange(0,1+d,d)]
            
    def to_rgba(self, v):
        d = self.d
        if self.mn in ['R','G','B','C','M','Y']:
            vn = max(0., (float(v)-self.vmin)/self.vmax)
            vn = min(1., vn)
            vni = int(vn/d)
            return self.rgba[vni]
        return self.sm.to_rgba(float(v))
    
    def RGB(self, v):
        # to RGB string (e.g. 255,0,0)
        rgba = [str(int(255*x)) for x in self.to_rgba(v)][:3]
        return ','.join(rgba)
        
####### Local Assembler ###############################################################


def detect_exons(sjpaths, offset, sja, exa, classifier=INTG):
    x = N.log2(sja+1)
    xd = (x[1:]-x[:-1])
    # use sjpaths to get donor/acceptor positions
    idxp = N.nonzero(xd>4)[0]+1 # donor(+ strand), acceptor(- strand)
    idxn = N.nonzero(xd<-4)[0]+1 # acceptor(+ strand), donor(- strand)
    tmp = [[x,'p',0] for x in idxp]+[[x,'n',0] for x in idxn]
    gaps0 = find_np_pairs(tmp, maxsize=10000)
    # tst: donor (+ strand), acceptor (- strand) => idxp
    # ted: => idxn
    idxp = set(sjpaths['tst'].values-offset)
    idxn = set(sjpaths['ted'].values-offset)
    tmp = [[x,'p',0] for x in idxp]+[[x,'n',0] for x in idxn]
    # gaps = find_np_pairs(tmp, xd)
    gaps1 = find_np_pairs(tmp)
    gaps = sorted(set(gaps0+gaps1))
    zoom = classifier.json['zoom']
    covfactor = classifier.json['th']
    def _gen_params():
        #find_maxgap = cyas2.find_maxgap
        for st,ed in gaps:
            # lemax, lgap, llen
            exsub = exa[st:ed]
            if len(exsub)==0:
                print(st,ed)
            emax = exsub.max()
            th = emax*covfactor
            lemax = N.log2(zoom*emax+1)
            lgap = N.log10(find_maxgap2(exsub, th)+1)
            llen = N.log10(ed-st+1)
            sdmax = max(xd[st-1],xd[ed-1])
            mp = float(N.sum(exsub>th))/len(exsub)
            yield (st,ed,lemax,lgap,llen,sdmax,mp)
    cols = ['ost','oed','lemax','lgap','llen','sdmax','mp']
    df = PD.DataFrame([x for x in _gen_params()], columns=cols)
    classifier.classify(df)
    return df
    
def find_np_pairs(tmp, maxsize=20000):
    # tmp=[ [pos, 'n' or 'p', 0], ...]
    # find n=>p pairs
    tmp = sorted(tmp)
    n = len(tmp)
    def _gen_sted():
        for i in range(n):
            x = tmp[i]
            if x[1]=='n':
                for j in range(i+1,n):
                    y = tmp[j]
                    if (y[0]-x[0]>maxsize):
                        break
                    if (y[1]=='p')&(y[0]>x[0]):
                        x[2],y[2]=1,1
                        yield (x[0],y[0])
                        break
        for j in range(n):
            y = tmp[j]
            if (y[1]=='p')&(y[2]==0): # unused
                for i in range(j-1,-1,-1):
                    x = tmp[i]
                    if (y[0]-x[0]>maxsize):
                        break
                    if (x[1]=='n')&(x[0]<y[0]):
                        x[2],y[2]=1,1
                        yield (x[0],y[0])
                        break
    return [x for x in _gen_sted()]
    
def find_genespan(st0, ed0, gaps, th=1):
    gaps2 = gaps[(gaps['exon']==False)&(gaps['sdmax']>th)]
    tmp = [[st0,'n',0],[ed0,'p',0]] + \
          [[x,'p',0] for x in set(gaps2['st'].values)] + \
          [[x,'n',0] for x in set(gaps2['ed'].values)]
    df = PD.DataFrame(find_np_pairs(tmp), columns=['st','ed'])
    df.sort_values(['st','ed'],inplace=True)
    return df
    
def detect_53(sja, exa, strand, classifier=E53C):
    zoom = classifier.json['zoom']
    if strand=='+':
        x = N.log2(zoom*sja+1)
    else:
        x = N.log2(zoom*sja[::-1]+1)
    xd = (x[1:]-x[:-1])
    xm = (x[1:]+x[:-1])/2.
    # idxp = N.nonzero(xd>0)[0] # source
    # idxn = N.nonzero(xd<0)[0] # sink
    sdiffth = classifier.json['sdiffth']
    idxp = N.nonzero(xd>sdiffth)[0] # source
    idxn = N.nonzero(xd<-sdiffth)[0] # sink
    # pos, sdiff, smean, kind
    n = len(x)
    if strand=='+':
        recs = [(i+1, N.abs(xd[i]), xm[i], '5') for i in idxp]+\
               [(i+1, N.abs(xd[i]), xm[i], '3') for i in idxn]
    else:
        recs = [(n-i-1, N.abs(xd[i]), xm[i], '5') for i in idxp]+\
               [(n-i-1, N.abs(xd[i]), xm[i], '3') for i in idxn]
        # position is off by one for negative strand (since xd is calculated in + strand)
    df = PD.DataFrame(recs, columns=['pos','sdiff','smean','kind'])
    classifier.classify(df)
    df = df[df['e53']==True].copy()
    # # genespans
    # idxp = df[df['kind']=='5']['pos']
    # idxn = df[df['kind']=='3']['pos']
    # if strand=='+':
    #     tmp = [[x,'n',0] for x in idxp]+[[x,'p',0] for x in idxn]
    # else:
    #     tmp = [[x,'p',0] for x in idxp]+[[x,'n',0] for x in idxn]
    # gspans = unionregion(find_np_pairs(tmp))
    return df#, gspans
    
# (~12sec ==> use cythonized ~ 3sec ==> use vectorized below ~ 1 sec
# def find_maxgap(arr, covfactor):
#     emax = arr.max()
#     th = emax*covfactor
#     emin = arr.min()
#     if (emin>th):
#         return 0
#     idx = N.nonzero(arr<=th)
#     if len(idx[0])==0:
#         return 0
#     cmax = 1
#     cst = idx[0][0]
#     ced = cst
#     for i in idx[0][1:]:
#         if i==ced+1: # continuous
#             ced = i
#         else:
#             cmax = max(cmax, ced-cst+1)
#             cst = ced = i
#     cmax = max(cmax, ced-cst+1)
#     return cmax    


def find_maxgap2(arr, th):
    idx = N.nonzero(arr<=th)[0]
    if len(idx)==0:
        return 0
    dif = idx[1:]-idx[:-1] # if continuous dif==1
    idx2 = N.nonzero(dif>1)[0] # non contiguous   
    idxst = N.array([idx[0]]+list(idx[idx2+1]))
    idxed = N.array(list(idx[idx2])+[idx[-1]])
    # gap start idx2[x]+1 ==> gap end idx2[x+1]
    gsize = idxed - idxst + 1
    return N.max(gsize)

    
def unionregion(df, sfld='st', efld='ed'):
    if type(df)==type([]):
        sted = sorted(df)
    else:
        sted = df.sort_values([sfld,efld])[[sfld,efld]].values
    def _gen():
        rec0 = sted[0]
        for rec1 in sted[1:]:
            if rec0[1] < rec1[0]: # ed1<st0 new interval
                yield rec0
                rec0 = rec1
            else: # overlapping/contiguous
                rec0[1] = max(rec0[1], rec1[1]) # update the end
        yield rec0
    recs = [x for x in _gen()]
    return recs
    
def fill_gap(sja, sj, exons, strand, offset):
    sjac = sja.copy()
    # need to fill exons inside the path as well
    if strand =='+':
        steds = sorted(set([tuple([int(z) for z in x.split(',')])  \
                                for y in sj['name'] 
                                for x in y.split('|')[1:-1]]))
    else:
        steds = sorted(set([tuple([int(z) for z in x.split(',')][::-1])  \
                                for y in sj['name'] 
                                for x in y.split('|')[1:-1]]))

    for st,ed in steds: # this way both 5' and 3' positions will be correct
        ost = st-offset
        oed = ed-offset
        sjac[ost:oed] = min(sjac[ost-1],sjac[oed+1])

    if len(exons)==0:
        return sjac
    #gaps0 = gaps[gaps['exon']==True]
    gaps1 = unionregion(exons, 'ost', 'oed')
    # if strand=='+': # >>> fill from right
    #     for st,ed in gaps1:
    #         sjac[st:ed] = sjac[st-1]
    # else: # <<< fill from left
    #     for st,ed in gaps1:
    #         sjac[st:ed] = sjac[ed+1]
    # return sjac
    # above only 5' positions are correct (3' are offset by exon length)
    for st,ed in gaps1: # this way both 5' and 3' positions will be correct
        sjac[st:ed] = min(sjac[st-1],sjac[ed+1])
    return sjac
    
def reversepathcode(pc):
    return ','.join(['|'.join(x.split('|')[::-1]) for x in pc.split(',')][::-1])

# def fix_unstranded_pathcode(sj, strand):
#     if strand=='-':
#         idx = sj['strand']=='.'
#         sj.loc[idx, 'name'] = [reversepathcode(x) for x in sj[idx]['name']]

def find_set(paths, sjs, chrom, st, ed, strand, covth=0.1):
    # idxs0 = (sjs['tst']>=st)&(sjs['ted']<=ed)&(sjs['strand'].isin([strand,'.']))&(sjs['chr']==chrom)
    idxs1 = (sjs['tst']>=st)&(sjs['ted']<=ed)&(sjs['strand'].isin(STRS[strand]))&(sjs['chr']==chrom)
    if N.sum(idxs1)==0:
        return []
    # sj = sjs[idxs0].copy()
    sj = sjs[idxs1].copy()
    idxp = (paths['tst']>=st)&(paths['ted']<=ed)&(paths['strand'].isin(STRS[strand]))&(paths['chr']==chrom)
    t = paths[idxp].copy()
    #fix_unstranded_pathcode(sj, strand)
    #fix_unstranded_pathcode(t, strand)
    sju = assign_jids1(t, strand)
    try:
        assign_jids2(sj, strand, sju)
        # assign_jids2(t, strand, sju)
    except:
        print(st,ed,strand)
        raise
    z = N.zeros((len(t),len(sj)))
    for sid1,jc1 in t[['sid','jc']].values:
        jc1b = ','+jc1+','
        for sid2,jc2 in sj[['sid','jc']].values:
            z[sid1-1][sid2-1] = ','+jc2+',' in jc1b
    th = t['tcov'].max()*covth
    sids0 = t[t['tcov']>=th]['sid'].values # initial set
    z0 = z[sids0-1].sum(axis=0) # sjpath usage
    if all(z0>0): # already every sjpaths covered
        return [t.set_index('sid').ix[sids0]]
    # go through the rest in the order of largest tcov until everything is covered
    sids1 = t[t['tcov']<th].sort_values('tcov', ascending=False)['sid'].values
    sids0 = list(sids0)
    cscore = N.sum(z0>0)
    tscore = z.shape[1]
    for sid in sids1:
        z0 = z0+z[sid-1]
        score = N.sum(z0>0)
        if score>cscore: # contribute
            sids0.append(sid)
            if score==tscore:
                return [t.set_index('sid').ix[sids0]]
            cscore = score
    return [t.set_index('sid').ix[sids0]]

def set_ad_pos(df, which='sj'):
    if which == 'sj':
        idx = df['strand'].isin(['-','.-'])
    else:
        idx = df['strand'].isin(['+','.+'])
    df.loc[idx, 'dpos'] = df[idx]['ed']
    df.loc[~idx,'dpos'] = df[~idx]['st']
    df.loc[idx, 'apos'] = df[idx]['st']
    df.loc[~idx,'apos'] = df[~idx]['ed']
    df['dpos'] = df['dpos'].astype(int)
    df['apos'] = df['apos'].astype(int)

def _pc(st,ed,strand,sep):
    if strand in ['+','.+','.']:
        return '{0}{2}{1}'.format(int(st),int(ed),sep)
    return '{1}{2}{0}'.format(int(st),int(ed),sep)

AFLD = {'ex':{'+':'st','-':'ed','.':'st'},
        'sj':{'+':'ed','-':'st','.':'ed'}}
DFLD = {'ex':{'+':'ed','-':'st','.':'ed'},
        'sj':{'+':'st','-':'ed','.':'st'}}
STRS = {'+':['+','.+'],
        '-':['-','.-'],
        '.':['.+','.-']}
EXDFCOLS = ['chr','st','ed','strand','name','kind','ecov']
SJDFCOLS = ['chr','st','ed','strand','name','kind','tcnt'  ]#,'donor','acceptor','dp','ap']
PATHCOLS = ['chr','st','ed','name','strand','tst','ted','tcov0','tcov1','tcov', 'tcov0a','tcov0b','tcov0c']

class LocalAssembler(object):
    
    def __init__(self, bwpre, chrom, st, ed, dstpre, 
                 sjbwpre=None,
                 refcode='gen9',
                 sjpaths=None, 
                 discardunstranded=False,
                 uth=1, 
                 mth=3, 
                 sjratioth=2e-3, 
                 usjratioth=1e-2,
                 #covfactor=0.05, 
                 covth=0.1,
                 upperpathnum=3000, # if num of paths larger than this increase stringency for sjs
                 pathcheckth=300, # above this num of sjs check sc1(ucnt)==0 if >50% remove
                 ):
        self.bname = '{0}:{1}-{2}'.format(chrom,st,ed)
        self.bwpre = bwpre
        self.sjbwpre = sjbwpre
        self.dstpre = dstpre
        self.discardunstranded = discardunstranded
        self.refcode = refcode # for classifier
        self.chrom = chrom
        self.st = int(st)
        self.ed = int(ed)
        self.uth = uth
        self.mth = mth
        self.sjratioth = sjratioth
        self.usjratioth = usjratioth
        #self.covfactor = covfactor
        self.covth = covth
        self.upperpathnum = upperpathnum
        self.pathcheckth = pathcheckth
        mixunstranded = not discardunstranded
        self.sjexbw = sjexbw = SjExBigWigs(bwpre, sjbwpre, mixunstranded=mixunstranded)
        self.arrs = arrs = {}
        with sjexbw: # get bw arrays
            for k in ['ex','sj']:
                arrs[k] = {}
                for s in ['+','-']:
                    arrs[k][s] = sjexbw.bws[k][s].get(chrom, st, ed)
        self._sjpaths=sjpaths
        self.load_classifiers()
        # self.load_and_filter_sjpaths(sjpaths) # => move to process

    def __str__(self):
        return 'LocalAssembler({0}:{1}-{2}, bwpre:{3}, dstpre:{4}'.format(self.chrom,self.st,self.ed,self.bwpre,self.dstpre)
    def _log(self, msg, level='debug'):
        _l = getattr(LOG, level)
        _l('{0}: {1}'.format(self.bname, msg))

    def loginfo(self, msg):
        self._log(msg, 'info')

    def logdebug(self, msg):
        self._log(msg, 'debug')
        
    def load_classifiers(self):
        refcode = self.refcode
        # [TODO] how to specify classifiers when assembling replicates?
        if (refcode is None) or (not UT.isstring(self.bwpre)): # load default
            self.intg=INTG
            self.e53c=E53C 
            self.ef5=EF5 
            self.ef3=EF3
            return

        # pathpre  from bwpre, sjbwpre
        # pathpre: bwpre+'.{refcode}'
        pathpre = self.bwpre+'.'+refcode
        path = pathpre+'.exonparams.json'
        if os.path.exists(path):
            with open(path,'r') as fp:
                self.exonparams = ep = json.load(fp)
            self.intg = LogisticClassifier(json=ep, dstcol='exon')        
        else:
            LOG.warning('{0} does not exists, reverting to default'.format(path))
            self.intg = INTG

        path = pathpre+'.gap5params.json'
        if os.path.exists(path):        
            with open(path,'r') as fp:
                self.gap5params = g5p = json.load(fp)        
            self.ef5 = EdgeFinder(g5p)
        else:
            LOG.warning('{0} does not exists, reverting to default'.format(path))
            self.ef5 = EF5 
        path = pathpre+'.gap3params.json'   
        if os.path.exists(path):        
            with open(path,'r') as fp:
                self.gap3params = g3p = json.load(fp)        
            self.ef3 = EdgeFinder(g3p)
        else:
            LOG.warning('{0} does not exists, reverting to default'.format(path))
            self.ef3 = EF3 

        if self.sjbwpre is not None:
            pathpre = self.sjbwpre+'.'+refcode
        path = pathpre+'.e53params.json'
        if os.path.exists(path):                
            with open(path,'r') as fp:
                self.e53params = e5p = json.load(fp)    
            self.e53c = LogisticClassifier(json=e5p, dstcol='e53')
        else:
            LOG.warning('{0} does not exists, reverting to default'.format(path))
            self.e53c = E53C

    def _read_sjpaths(self):
        sjpaths0 = self._sjpaths
        if sjpaths0 is not None:
            return sjpaths0

        chrom,st,ed = self.chrom,self.st,self.ed
        if UT.isstring(self.bwpre):
            chrfiltered = self.bwpre+'.sjpath.{0}.filtered.bed.gz'.format(chrom)
            if os.path.exists(chrfiltered):
                sj = GGB.read_bed(chrfiltered)
            else:
                chrpath = self.bwpre+'.sjpath.{0}.bed.gz'.format(chrom) # separated chrom file exists?
                if os.path.exists(chrpath):
                    sj = GGB.read_bed(chrpath)
                else:
                    sj = GGB.read_bed(self.bwpre+'.sjpath.bed.gz')
            idx0 = (sj['chr']==chrom)&(sj['tst']>=st)&(sj['ted']<=ed)        
            sj0 = sj[idx0].copy()
            # merged sjpath has 53exon in pathcode => remove
            if len(sj0)>0:
                name0 = sj0.iloc[0]['name']
                if len(name0.split('|'))<len(name0.split(',')):
                    sj0['name'] = [','.join(x.split(',')[1:-1]) for x in sj0['name']]
            return sj0

        # list of bwpres, load and merge
        sjps0 = [GGB.read_bed(b+'.sjpath.bed.gz') for b in self.bwpre]
        sjps = []
        for sj in sjps0:
            idx0 = (sj['chr']==chrom)&(sj['tst']>=st)&(sj['ted']<=ed)        
            sj0 = sj[idx0].copy()
            if len(sj0)>0:
                name0 = sj0.iloc[0]['name']
                if len(name0.split('|'))<len(name0.split(',')):
                    sj0['name'] = [','.join(x.split(',')[1:-1]) for x in sj0['name']]            
            sjps.append(sj0)
        sjp = PD.concat(sjps, ignore_index=True)
        sjg = sjp.groupby(['chr','name'])
        sj = sjg.first()
        # chr,st,ed,name,sc1,strand,tst,ted,sc2,#exons,esizes,estarts
        sj['st'] = sjg['st'].min().astype(int)
        sj['ed'] = sjg['ed'].max().astype(int)
        sj['sc1'] = sjg['sc1'].sum()
        sj['sc2'] = sjg['sc2'].sum()
        self._sjps = sjps 
        sj = sj.reset_index()
        return sj

    def load_and_filter_sjpaths(self):
        sjpaths0 = self._read_sjpaths()
        bwpre = self.bwpre
        chrom,st,ed = self.chrom,self.st,self.ed
        uth,mth = self.uth,self.mth
        sjratioth = self.sjratioth
        usjratioth = self.usjratioth

        #idx0 = (sjpaths0['chr']==chrom)&(sjpaths0['tst']>=st)&(sjpaths0['ted']<=ed)
        #s0 = sjpaths0[idx0]
        self._sjpaths0 = s0 = sjpaths0
        idxu = s0['strand']=='.' # unstranded => duplicates into '.+','.-'
        if self.discardunstranded:
            s1 = s0[~idxu].copy()
        else:
            if N.sum(idxu)>0:
                # if unstranded junction share donor/acceptor sites with stranded one then assign that strand
                sj_pn = s0[~idxu]
                if len(sj_pn)>0:
                    tst2str = UT.df2dict(sj_pn, 'tst', 'strand')
                    ted2str = UT.df2dict(sj_pn, 'ted', 'strand')
                    sj_u = s0[idxu].copy()
                    sj_u['strand'] = [tst2str.get(x,y) for x,y in sj_u[['tst','strand']].values]
                    sj_u['strand'] = [ted2str.get(x,y) for x,y in sj_u[['ted','strand']].values]
                    idx_n = sj_u['strand']=='-'
                    sj_u.loc[idx_n, 'name'] = [reversepathcode(x) for x in sj_u[idx_n]['name']]
                    idxu2 = sj_u['strand']=='.'
                    sj_upn = sj_u[~idxu2]
                    self.logdebug('{0} unstranded assigned +/- with matching ends'.format(len(sj_upn)))
                    sj_up = sj_u[idxu2].copy()
                    sj_un = sj_u[idxu2].copy()
                else:
                    sj_upn = None
                    sj_up = s0[idxu].copy()
                    sj_un = s0[idxu].copy()

                sj_up['strand'] = '.+'
                sj_un['name'] = [reversepathcode(x) for x in sj_un['name']]
                sj_un['strand'] = '.-'
                s1 = PD.concat([sj_pn, sj_upn, sj_up, sj_un], ignore_index=True)
            else:
                s1 = s0
            
        self.sjpaths0 = s1 # original unfiltered (just position restricted)
        sc1 = self.sjpaths0['sc1']
        if self.sjpaths0.dtypes['sc2']=='O': # merged sjpath
            self.sjpaths0['sc2'] = sc1
        sc2 = self.sjpaths0['sc2']
        idx1 = (sc1>=uth)|(sc2-sc1>=mth)
        self.sjpaths = sjpaths = self.sjpaths0[idx1].copy()
        # max ratio to cov (sj+ex) > sjratioth
        with self.sjexbw:
            sjaa = self.sjexbw.bws['sj']['a'].get(chrom, st, ed)
            exaa = self.sjexbw.bws['ex']['a'].get(chrom, st, ed)
        a = sjaa+exaa # all of the coverages
        o = int(self.st)
        # sjpaths['minscov'] = [N.min(a[s-o:e-o]) for s,e in sjpaths[['tst','ted']].values]]
        sjpaths['sjratio'] = [x/N.min(a[int(s-o):int(e-o)]) for x,s,e in sjpaths[['sc2','tst','ted']].values]
        # .values => dtype float matrix => s,e float
        n0 = len(sjpaths)
        idxpn = (sjpaths['strand'].isin(['+','-']))&(sjpaths['sjratio']>sjratioth)
        idxu = (sjpaths['strand'].isin(['.+','.-']))&(sjpaths['sjratio']>usjratioth)
        self.sjpaths = sjpaths[idxpn|idxu].copy()
        n1 = len(self.sjpaths)
        self.loginfo('sjratio filter: {0}=>{1}'.format(n0,n1))

        
    def find_exons(self):
        #covfactor=self.covfactor
        arrs = self.arrs
        self.filled  = {}
        self.exons  = {}
        self.gspans = {}
        self.gaps = {}
        sjs = self.sjpaths
        for s in ['+','-']:
            sja = arrs['sj'][s]
            exa = arrs['ex'][s]
            sj = sjs[sjs['strand'].isin(STRS[s])]
            df = detect_exons(sj, self.st, sja, exa, classifier=self.intg)
            self.exons[s] = df[df['exon']==True].copy()            
            self.gaps[s] = df
            self.filled[s] = fill_gap(sja, sj, self.exons[s], s, self.st)
            self.gspans[s] = self._get_spans(s)

    
    def _get_spans(self, strand):
        sj0 = self.sjpaths
        ex0 = self.exons[strand]
        o = self.st
        sj1a = sj0[sj0['strand'].isin(STRS[strand])]
        sj1 = [list(x-o) for x in sj1a[['st','ed']].values]
        ex1 = [list(x) for x in ex0[['ost','oed']].values]
        arr = sj1+ex1
        if len(arr)==0:
            return []
        return UT.union_contiguous_intervals(arr)
            
    def find_53(self):
        self.e53pos = {}
        for s in ['+','-']:
            sja = self.filled[s]
            exa = self.arrs['ex'][s]            
            self.e53pos[s] = detect_53(sja, exa, s, classifier=self.e53c)
        
    def find_paths(self, ost, oed, strand):
        sj,ex = self._get_sub_sjex(ost,oed,strand)
        if len(sj)==0:
            return None
        # fix_unstranded_pathcode(sj, strand)
        self.gg = gg = GeneGraph(ex, sj, strand, upperpathnum=self.upperpathnum)
        df = gg.find_all_paths()
        # fix_unstranded_pathcode(sj, strand)            
        return df    
        
    def _get_sub_sjex(self, ost, oed, strand):
        #st,ed = gspan.iloc[0]
        gaps = self.exons[strand]
        idx = (gaps['ost']>=ost)&(gaps['oed']<=oed)
        ex = gaps[idx].copy()
        offset = self.st
        ex['st'] = ex['ost']+offset
        ex['ed'] = ex['oed']+offset
        st = ost+offset
        ed = oed+offset
        sjpaths = self.sjpaths
        idx0 = (sjpaths['chr']==self.chrom)&(sjpaths['strand'].isin(STRS[strand]))
        idx1 = (sjpaths['tst']>=st)&(sjpaths['ted']<=ed)
        idx2 = (sjpaths['tst']>=st)&(sjpaths['tst']<=ed)
        idx3 = (sjpaths['ted']>=st)&(sjpaths['ted']<=ed)
        idx4 = (sjpaths['chr']==self.chrom)&(sjpaths['strand']==strand)
        #idx5 = sjpaths['sc1']>th
        idx = (idx0&idx1)|((idx2|idx3)&idx4) #&idx5)
        sj = sjpaths[idx]
        n0 = len(sj)
        n1 = self.pathcheckth
        if n0>n1:
            n2 = N.sum(sj['sc1']==0)
            if n2>100:
                LOG.warning('num sj ({0}>{1}) removed non unique junctions({2})'.format(n0,n1,n2))
                idx5 = sj['sc1']>0
                sj = sj[idx5]
        sj = sj.copy()
        return sj, ex
               
    def find_all_paths(self):
        pathdfs = []
        for strand in ['+','-']:
            #spans = self._get_spans(strand)
            #for st,ed in spans:
            for ost,oed in self.gspans[strand]:
                pathdfs.append(self.find_paths(ost,oed,strand))
        df = PD.concat(pathdfs, ignore_index=True)
        #n0 = len(df)
        #df = df.groupby('name').first().reset_index()
        n1 = len(df)
        # '.' part of others?
        df1 = df[df['strand'].isin(['+','-'])]
        df0 = df[df['strand'].isin(['.+','.-'])]
        allpathcode = '$'.join(df1['name'].values)
        idx = N.array([reversepathcode(x) not in allpathcode for x in df0['name']])
        df2 = df0[idx]
        df = PD.concat([df1,df2],ignore_index=True)
        n2 = len(df)
        self.logdebug('n1(unique pc)={0}, n2(unstrand dup removed)={1}'.format(n1,n2))
        df['chr'] = self.chrom
        self.allpaths = df
        return df
        
    def decompose_paths(self):
        # allpaths => junctions
        ap = self.allpaths
        dfs = []
        dfe = []
        chrom = self.chrom
        def _sgen():
            sted = set()
            for strand in ['+','.+', '-','.-']:
                for p in ap[ap['strand']==strand]['name'].values:
                    for x in p.split(','):
                        st,ed = [int(y) for y in x.split('|')]
                        if st>ed:
                            st,ed = ed,st
                        if (st,ed,strand[-1]) not in sted:
                            yield (chrom,st,ed,strand,x,'j')
                            sted.add((st,ed,strand))
        def _egen():
            sted = set()
            for strand in ['+','.+', '-','.-']:
                for p in ap[ap['strand']==strand]['name'].values:
                    for x in p.split('|')[1:-1]:
                        st,ed = [int(y) for y in x.split(',')]
                        if st>ed:
                            st,ed = ed,st
                        if (st,ed,strand[-1]) not in sted:
                            yield (chrom,st,ed,strand,x,'i')
                            sted.add((st,ed,strand[-1]))
        cols = ['chr','st','ed','strand','name','kind']
        sjdf = PD.DataFrame([x for x in _sgen()], columns=cols)
        exdfi = PD.DataFrame([x for x in _egen()], columns=cols)
        set_ad_pos(sjdf, 'sj')
        set_ad_pos(exdfi, 'ex')
        self.sjdf = sjdf
        self.exdfi = exdfi

        
    def _find_pos(self, strand, direction):
        ap = self.allpaths
        e53 = self.e53pos
        ef5 = self.ef5
        ef3 = self.ef3
        sjdf = self.sjdf
        gspans = self.gspans
        st0 = self.st
        strands = {'+':['+','.+'],'-':['-','.-']}
        sj1 = sjdf[sjdf['strand'].isin(strands[strand])]
        ps = ap[ap['strand'].isin(strands[strand])]
        es = e53[strand]# in local coord
        # direction < (left most) ('+': 5', '-':3' )
        if direction=='<':
            pos_l = set([x-st0 for x in ps['tst'].values])
            sj_l = set([x-st0 for x in sj1['st'].values])
            kind = '5' if strand=='+' else '3'
            es_l = set(es[es['kind']==kind]['pos'].values)
            pos_l.update(es_l.intersection(sj_l))
            return pos_l
        # direction > (right most) ('+':3', '-':5' )
        pos_r = set([x-st0 for x in ps['ted'].values])
        sj_r = set([x-st0 for x in sj1['ed'].values])
        kind = '3' if strand=='+' else '5'
        es_r = set(es[es['kind']==kind]['pos'].values)
        pos_r.update(es_r.intersection(sj_r))
        return pos_r

    def _find_pos0(self, pos, strand, direction):
        gs = N.array(self.gspans[strand])
        if direction=='<': # look for last gspan end
            tgt = gs[:,1]
            i0 = bisect.bisect_left(tgt, pos)
            if i0==0:
                return 1
            return tgt[i0-1]
        tgt = gs[:,0] # look for next gspan start
        i0 = bisect.bisect_right(tgt, pos)
        if i0==len(tgt):
            return len(self.arrs['sj'][strand])-1
        return tgt[i0]

    def _subtract_exons(self, pos, pos0, sja, exa, exs, direction):
        # direction <
        # find exons between pos0=>pos, subtract sja corresponding to further
        # edge of the exon (at position st)
        if direction=='<':
            sja1 = sja[pos0-1:pos+1].copy() # index offset by 1
            exa1 = exa[pos0:pos+1].copy()
            ex = exs[(exs['ost']>=pos0)&(exs['oed']<=pos)]
            for st,ed in ex[['ost','oed']].values:
                st0 = st-pos0
                ed0 = ed-pos0
                exa1[st0:ed0] = exa1[st0:ed0]-(sja1[st0]-sja1[st0+2]) # st0-1 but array itself is offset by 1
            return sja1[1:], exa1 # same length
        # direction >
        # find exons between pos=>pos0, subtract sja corresponding to further
        # edge of the exon (at position ed)
        sja1 = sja[pos-1:pos0+2].copy() # same index but get one more on right
        exa1 = exa[pos-1:pos0].copy()
        ex = exs[(exs['ost']>=pos)&(exs['oed']<pos0)] # don't include = for ed
        for st,ed in ex[['ost','oed']].values:
            st0 = st-pos+1
            ed0 = ed-pos+1
            exa1[st0:ed0] = exa1[st0:ed0]-(sja1[ed0+2]-sja1[ed0])
        return sja1[:len(exa1)], exa1 # same length
        
    def find_53exons(self):
        # edges of all paths + 53 detected above if corresponding node exists in paths
        ef5 = self.ef5
        ef3 = self.ef3
                
        EF = {'<':{'+':ef5, '-':ef3},
              '>':{'+':ef3, '-':ef5}}
        KIND = {'<':{'+':'5','-':'3'},
                '>':{'+':'3','-':'5'}}
        #self._arr = {}
        self._poss = {}
        def _gen(strand):
            sja = self.arrs['sj'][strand]
            exa = self.arrs['ex'][strand]
            exs = self.exons[strand]
            for direction in ['<','>']:
                poss = self._find_pos(strand, direction)
                self._poss[(strand,direction)] = poss
                k = KIND[direction][strand]
                #print('poss:', poss)
                for pos in poss:
                    pos0 = self._find_pos0(pos, strand, direction)
                    sja1, exa1 = self._subtract_exons(pos, pos0, sja, exa, exs, direction)
                    #self._arr[(strand,direction,pos,pos0)]=(sja1,exa1)
                    epos = pos+EF[direction][strand].find(sja1, exa1, direction)
                    #print('  s:{2},d{3},:pos0:{0},pos:{1},epos:{4}'.format(pos0,pos,strand,direction,epos))
                    if pos < epos:
                        yield (pos, epos, k)
                    elif pos > epos:
                        yield (epos, pos, k)
                    else: # pos == epos
                        self.logdebug("didn't find 53exon: strand:{0}, pos:{1}, pos0:{2}, direction:{3}".format(strand,pos,pos0,direction))
        dfs = {}
        for s in ['+','-']:
            dfs[s] = PD.DataFrame([x for x in _gen(s)], columns=['ost','oed','kind'])
        self.exons53 = dfs
        return dfs
        
    def organize(self):
        # dataframes for 
        # 1) exons 
        # 2) junctions
        # 3) unused sjpaths
        # 4) allpaths with 5-3 exons attached
        #    - generate paths with internal 5-3 exons 
        #    - paths grouped by 5-3 end (assign id)
        # 1 ex
        cols0 = ['ost','oed','strand','kind']
        cols = ['chr','st','ed','strand','kind','name', 'apos','dpos']
        # 53 exons
        exs = []
        for s in ['+','-']:
            e53 = self.exons53[s]
            e53['strand'] = s
            exs += [e53[cols0]]
        exdf53 = PD.concat(exs, ignore_index=True).copy()
        exdf53['chr'] = self.chrom
        exdf53['st'] = (exdf53['ost']+self.st).astype(int)
        exdf53['ed'] = (exdf53['oed']+self.st).astype(int)
        exdf53['name'] = [_pc(st,ed,s,',') for st,ed,s in exdf53[['st','ed','strand']].values]
        set_ad_pos(exdf53, 'ex')
        # internal exons <= construct from paths
        exdfi = self.exdfi
        if len(exdfi)>0:
            exdf = PD.concat([exdfi[cols], exdf53[cols]], ignore_index=True)
            exdf['st'] = exdf['st'].astype(int)
            exdf['ed'] = exdf['ed'].astype(int)
        else:
            exdf = exdf53
        # 2 sj
        sjdf = self.sjdf
        # 4 passes, chr, name, strand, tst, ted => add st, ed (as 5,3 edge)
        # first go through all paths, record used 53
        usedids = set()
        paths = []
        pcols = ['chr','st','ed','name','strand','tst','ted']
        ap = self.allpaths
        def _one(s,sfld,efld):
            tstfld,tedfld='t'+sfld,'t'+efld
            apsub = ap[ap['strand'].isin(s)]
            e5 = exdf[(exdf['strand'].isin(s))&(exdf['kind']=='5')].reset_index()[['st','ed','index']]
            e3 = exdf[(exdf['strand'].isin(s))&(exdf['kind']=='3')].reset_index()[['st','ed','index']]
            e5p = PD.merge(e5,apsub,left_on=efld, right_on=tstfld) # inner join
            usedids.update(e5p['index'].values)
            e5p = e5p[[sfld,'chr','name','strand','tst','ted']]
            e3.rename(columns={sfld:sfld+'3'},inplace=True)
            e5p3 = PD.merge(e5p, e3, left_on=tedfld,right_on=sfld+'3')
            usedids.update(e5p3['index'].values)
            e5p3['name'] = ['{0},{1},{2}'.format(s,n,e) for s,n,e in e5p3[[sfld,'name',efld]].values]
            return e5p3[pcols]
        e53p = _one(['+','.+'], 'st', 'ed')
        e53n = _one(['-','.-'], 'ed', 'st')
        # go through unused 53, find paths with matching acceptor/donor make new paths
        uu= exdf.ix[set(exdf.index.values).difference(usedids)]
        def _two(s, e, sfld, efld):
            tstfld,tedfld='t'+sfld,'t'+efld
            u5 = uu[(uu['strand'].isin(s)&(uu['kind']=='5'))]
            u3 = uu[(uu['strand'].isin(s)&(uu['kind']=='3'))]
            # 5' find paths in e that has donor matching u5[efld] attach => d
            sidx = pcols.index(sfld) # 1 or 2
            eidx = pcols.index(efld) # 2 or 1
            tsidx = pcols.index(tstfld) # 5 or 6
            teidx = pcols.index(tedfld) # 6 or 5
            recs = []
            names = set()
            for rec in e[pcols].values:# pcols
                recs.append(rec)
                path = rec[3] # name 
                names.add(path)
                for st, ed in u5[[sfld,efld]].values:
                    donor = ',{0}|'.format(ed)
                    m = path.find(donor)
                    name = str(st)+path[m:]
                    if (m>0) and (name not in names):
                        rec1 = rec.copy() # modify, sfld, tstfld, name
                        rec1[3] = name
                        rec1[sidx] = st
                        rec1[tsidx] = ed
                        recs.append(rec1)
                        names.add(name)
            # 3' find paths in e + d that has acceptor matching u3[sfld] => f
            recs2 = []
            names = set()
            for rec in recs:
                recs2.append(rec)
                path = rec[3]
                names.add(path)
                for st, ed in u3[[sfld,efld]].values:
                    acceptor='|{0},'.format(st)
                    m = path.find(acceptor)
                    name = path[:m]+'|{0},{1}'.format(st,ed)
                    if (m>0) and (name not in names):
                        rec1 = rec.copy()
                        rec1[3] = name
                        rec1[eidx] = ed
                        rec1[teidx] = st
                        recs2.append(rec1)
                        names.add(name)
            return recs2
        recs = _two(['+','.+'], e53p, 'st', 'ed')+ _two(['-','.-'], e53n, 'ed', 'st')
        self.paths = paths = PD.DataFrame(recs, columns=pcols).sort_values(['st','ed'])        
        # groupby st, ed and assign ids? ==> they can be immediately group by [st,ed] so no need
        self.exdf = exdf[cols]
        self.sjdf = sjdf[cols]
        # 3 unused sj
        sjpaths0 = self.sjpaths0
        allpathcode = '$'.join(paths['name'].values)
        idxused = N.array([y in allpathcode for y in sjpaths0['name']])
        self.unusedsj = sjpaths0[~idxused]         
        self.usedsj = sjpaths0[idxused]
        
    def calculate_scovs(self):
        sj = self.sjdf
        sj0 = self.sjpaths0
        sj0mat = sj0[['sc2','name']].values
        sj['tcnt'] = [N.sum([x for x,p in sj0mat if y in p]) for y in sj['name']]
        idx = sj['tcnt']==0
        tmp0 = ['{1}|{0}'.format(*y.split('|')) for y in sj[idx]['name']]
        tmp1 = [N.sum([x for x,p in sj0mat if y in p]) for y in tmp0]
        sj.loc[idx, 'tcnt'] = tmp1
        # sj['donor'] = ['{0}|'.format(y.split('|')[0]) for y in sj['name']]
        #sj['acceptor'] = ['|{0}'.format(y.split('|')[1]) for y in sj['name']]

    def calculate_branchp(self, jids, eids):
        sj0 = self.sjdf
        sj = sj0.set_index('name').ix[jids].reset_index()
        dsump = sj.groupby('dpos')['tcnt'].sum().astype(float)
        jdp = sj['tcnt'].values/dsump.ix[sj['dpos'].values].values
        j2p = dict(zip(sj['name'].values, jdp))
        # exon groupby acceptor
        ex0 = self.exdf
        ex = ex0.set_index('name').ix[eids].reset_index().copy()
        idxz = ex['ecov']==0
        ex.loc[idxz, 'ecov'] = 1e-6
        asump = ex.groupby('apos')['ecov'].sum().astype(float)
        edp = ex['ecov'].values/asump.ix[ex['apos'].values].values
        if (N.sum(N.isnan(edp))>0):
            self._cblcls = locals()
            raise
        e2p = dict(zip(ex['name'].values, edp))
        return j2p, e2p

        
    def tcov_by_nnls(self, s, e, strand):
        o = self.st
        p = self.paths
        idx = (p['tst']>=s)&(p['ted']<=e)&(p['strand'].isin(STRS[strand]))
        ps = p[idx]
        if len(ps)==0:
            return None
        pg = ps.groupby(['tst','ted']).first().reset_index()[['chr','tst','ted','strand','name']].sort_values(['tst','ted'])
        pg['strand'] = strand
        ne = len(pg)
        exa = self.arrs['ex'][strand]
        sja = self.arrs['sj'][strand]
        ne2ecov = self._ne2ecov
        def cov0(s,e):
            # return N.sum(sja[s-o:e-o]+exa[s-o:e-o])/(e-s)
            return N.mean(sja[s-o:e-o])
        def cov1s(s):
            s0 = max(0, s-o-10)
            s1 = max(s0+1,s-o)
            return N.mean(exa[s0:s1])
        def cov1e(e):
            return N.mean(exa[e-o:e-o+10])
        def cov2s(s):
            s0 = max(0, s-o-1)
            return sja[s-o]-sja[s0]
        def cov2e(e):
            e0 = max(0, e-o-1)
            return sja[e-o]-sja[e0]
        # if strand in ['+','.+']:
        #     def cov(s,e):
        #         return sja[s-o]
        # else:
        #     def cov(s,e):
        #         return sja[e-o-1]

        # cov0
        if ne>1:
            pg.rename(columns={'tst':'st','ted':'ed'}, inplace=True)
            pg['eid'] = N.arange(len(pg))
            ci = UT.chopintervals(pg, idcol='eid')
            ci['cov'] = [cov0(s,e) for s,e in ci[['st','ed']].values]
            ci['name1'] = ci['name'].astype(str).apply(lambda x: [int(y) for y in x.split(',')])    
            nc = len(ci)
            mat = N.zeros((nc,ne))
            for i,n1 in enumerate(ci['name1'].values):# fill in rows
                N.put(mat[i], N.array(n1), 1)
            ecov,err = nnls(mat, ci['cov'].values)
            pg['tcov0'] = ecov
            pg.rename(columns={'st':'tst','ed':'ted'}, inplace=True)
        else:
            s,e = pg.iloc[0][['tst','ted']]
            pg['tcov0'] = cov0(s,e)
        keys = [tuple(x) for x in p[idx][['tst','ted']].values]
        p.loc[idx, 'tcov0a'] = pg.set_index(['tst','ted']).ix[keys]['tcov0'].values
        # cov1, cov2
        if ne>1:
            sts = sorted(set(pg['tst'].values))
            eds = sorted(set(pg['ted'].values))
            nst,ned = len(sts),len(eds)
            mat = N.array([(pg['tst']==x).values for x in sts]+[(pg['ted']==x).values for x in eds], dtype=float)
            ### use ecov
            # tmp = [x.split('|') for x in pg['name']]
            # pg['5ename'] = [x[0] for x in tmp]
            # pg['3ename'] = [x[-1] for x in tmp]
            # pg['5ecov'] = [ne2ecov[x] for x in pg['5ename']]
            # pg['3ecov'] = [ne2ecov[x] for x in pg['3ename']]
            # if strand in ['+','.+']:
            #     stecovs = list(pg.groupby('tst')['5ecov'].mean().ix[sts])
            #     edecovs = list(pg.groupby('ted')['3ecov'].mean().ix[eds])
            # else:
            #     stecovs = list(pg.groupby('tst')['3ecov'].mean().ix[sts])
            #     edecovs = list(pg.groupby('ted')['5ecov'].mean().ix[eds])
            # c = N.array(stecovs+edecovs)
            ### end use ecov
            c = N.array([cov1s(x) for x in sts]+[cov1e(x) for x in eds])
            # enforce flux conservation: scale up 5'
            stsum = N.sum(c[:nst])
            edsum = N.sum(c[nst:])
            if strand in ['+','.+']:
                c[:nst] = (edsum/(stsum+1e-6))*c[:nst]
            else:
                c[nst:] = (stsum/(edsum+1e-6))*c[nst:]
            ecov,err = nnls(mat, c)
            pg['tcov1'] = ecov

            mat = N.array([(pg['tst']==x).values for x in sts]+[-1*(pg['ted']==x).values for x in eds], dtype=float)
            c = N.array([cov2s(x) for x in sts]+[cov2e(x) for x in eds])
            # enforce flux conservation: scale up 5'
            stsum = N.sum(c[:nst])
            edsum = N.sum(c[nst:])
            if strand in ['+','.+']:
                c[:nst] = ((-1*edsum)/(stsum+1e-6))*c[:nst]
            else:
                c[nst:] = ((-1*stsum)/(edsum+1e-6))*c[nst:]
            ecov,err = nnls(mat, c)
            pg['tcov2'] = ecov
        else:
            s,e = pg.iloc[0][['tst','ted']]
            pg['tcov1'] = (cov1s(s)+cov1e(e))/2.
            pg['tcov2'] = (cov2s(s)-cov2e(e))/2.

        p.loc[idx, 'tcov0b'] = pg.set_index(['tst','ted']).ix[keys]['tcov1'].values
        p.loc[idx, 'tcov0c'] = pg.set_index(['tst','ted']).ix[keys]['tcov2'].values
        pg['tcov'] = pg[['tcov0','tcov1','tcov2']].mean(axis=1)
        pg.loc[pg['tcov']<0,'tcov'] = 0 # shouldn't really happen
        p.loc[idx, 'tcov0'] = pg.set_index(['tst','ted']).ix[keys]['tcov'].values
        return pg[['chr','tst','ted','strand', 'tcov']]
        
    def tcov_by_utr(self, tst, ted, strand, tcov0):
        p = self.paths
        #cols = ['chr','st','ed','name','strand','tst','ted']
        idx = (p['tst']==tst)&(p['ted']==ted)&(p['strand'].isin(STRS[strand]))
        ps = p[idx]
        if len(ps)==0:
            return None
        pg = ps.groupby(['st','ed']).first().reset_index()[['chr','st','ed','strand']].sort_values(['st','ed'])
        pg['strand'] = strand
        # left and right edges
        lpos = sorted(set(pg['st'].values))
        rpos = sorted(set(pg['ed'].values))
        exa = self.arrs['ex'][strand]
        def cov(s,e):
            o = self.st
            return N.mean(exa[s-o:e-o])
        if len(lpos)>1: # left most up to tst [(lpos0,lpos1),(lpos1,lops2),...,(lpos(n-1),tst)]
            lcov0 = {} # cumulative sum
            for s,e in zip(lpos, lpos[1:]+[tst]):
                lcov0[s] = cov(s,e)
            lcov = {} # individual
            for s,sprev in zip(lpos[1:],lpos[:-1]):
                lcov[s] = max(0, lcov0[s]-lcov0[sprev])
            lcov[lpos[0]] = lcov0[lpos[0]]
            lsum = float(N.sum([x for x in lcov.values()]))
            if lsum==0:
                lsum += 1.
            lp = {s:lcov[s]/lsum for s in lpos}
        else:
            lp = {lpos[0]: 1.}
        if len(rpos)>1: 
            rcov0 = {} # cumulative sum
            for s,e in zip([ted]+rpos[:-1], rpos):
                rcov0[e] = cov(s,e)
            rcov = {} # individual
            for e,enext in zip(rpos[:-1],rpos[1:]):
                rcov[e] = max(0, rcov0[e]-rcov0[enext])
            rcov[rpos[-1]] = rcov0[rpos[-1]]
            rsum = float(N.sum([x for x in rcov.values()]))
            if rsum==0:
                rsum += 1.
            rp = {e:rcov[e]/rsum for e in rpos}
        else:
            rp = {rpos[0]: 1.}
        if len(lpos)>1:
            if len(rpos)>1:
                pg['tcov'] = [tcov0*lp[s]*rp[e] for s,e in pg[['st','ed']].values]
            else:
                pg['tcov'] = [tcov0*lp[s] for s in pg['st']]
        else:
            if len(rpos)>1:
                pg['tcov'] = [tcov0*rp[e] for e in pg['ed']]
            else:
                pg['tcov'] = tcov0

        keys = [tuple(x) for x in p[idx][['st','ed']].values]
        p.loc[idx, 'tcov1'] = pg.set_index(['st','ed']).ix[keys]['tcov'].values
        return pg[['chr','st','ed','strand','tcov']]
        
    def tcov_by_branchp(self, st, ed, tst, ted, strand, tcov0):
        p = self.paths
        idx = (p['st']==st)&(p['ed']==ed)&(p['strand'].isin(STRS[strand]))&\
              (p['tst']==tst)&(p['ted']==ted)
        if N.sum(idx)==0:
            return
        if N.sum(idx)>1:
            # calculate branchp within this group
            jids = set()
            eids = set()
            for n in p[idx]['name']:
                jids.update(n.split(',')[1:-1])
                eids.update(n.split('|'))
            j2p, e2p = self.calculate_branchp(jids, eids)
            def _prob(y):
                epath = y.split('|')[1:-1]
                jpath = y.split(',')[1:-1]
                return N.prod([e2p[x] for x in epath])*N.prod([j2p[x] for x in jpath])
            p.loc[idx,'tcov'] = [tcov0*_prob(y) for y in p[idx]['name']]
        else:
            p.loc[idx,'tcov'] = tcov0
        
        
    def estimate_abundance(self):
        # 1) 5-3 group by NNLS
        # 2) UTR difference
        # 3) within 5-3 group by tree branch prob
        paths = self.paths
        sjs = self.usedsj
        th = self.covth
        trimmed = []
        for s in ['+','-']:
            ps = paths[paths['strand'].isin(STRS[s])]
            if len(ps)==0:
                continue
            for chrom,st,ed in UT.union_contiguous(ps[['chr','st','ed']],returndf=False):
                pg1 = self.tcov_by_nnls(st,ed,s)
                if pg1 is not None:
                    for chrom,tst1,ted1,strand1,tcov1 in pg1.values:
                        pg2 = self.tcov_by_utr(tst1,ted1,strand1,tcov1)
                        if pg2 is not None:
                            for chrom,st2,ed2,strand2,tcov2 in pg2.values:
                                self.tcov_by_branchp(st2,ed2,tst1,ted1,strand2,tcov2)
                # trim paths
                trimmed += find_set(paths, sjs, chrom, st, ed, s, th)
        if len(trimmed)>0:
            self.tpaths = PD.concat(trimmed, ignore_index=True)
        else:
            self.tpaths = UT.make_empty_df(paths.columns)
        
    def calculate_ecovs(self):
        paths = self.paths
        ex = self.exdf
        o = self.st
        if len(ex)==0:
            return
        ex['ecov'] = N.nan
        for strand in ['+','-']:
            ps = paths[paths['strand'].isin(STRS[strand])]
            if len(ps)==0:
                continue
            for chrom,st,ed in UT.union_contiguous(ps[['chr','st','ed']],returndf=False):
                idx = (ex['st']>=st)&(ex['ed']<=ed)&(ex['strand'].isin(STRS[strand]))
                es = ex[idx].copy().sort_values(['st','ed'])
                es['tmpeid'] = N.arange(len(es))
                ne = len(es)
                exa = self.arrs['ex'][strand]
                #sja = self.arrs['sj'][strand]
                def cov(s,e):
                    return N.mean(exa[s-o:e-o])
                if ne>1:
                    ci = UT.chopintervals(es, idcol='tmpeid', sort=False)
                    ci['cov'] = [cov(s,e) for s,e in ci[['st','ed']].values]
                    ci['name1'] = ci['name'].astype(str).apply(lambda x: [int(y) for y in x.split(',')])    
                    nc = len(ci)
                    mat = N.zeros((nc,ne))
                    for i,n1 in enumerate(ci['name1'].values):# fill in rows
                        N.put(mat[i], N.array(n1), 1)
                    ecov,err = nnls(mat, ci['cov'].values)
                    ex.loc[idx,'ecov'] = ecov
                elif ne==1:
                    s,e = es.iloc[0][['st','ed']]
                    ex.loc[idx,'ecov'] = cov(s,e)
        # idxnan = ex['ecov'].isnull()
        # ex.loc[idxnan,'ecov'] = [cov(s,e) for s,e in ex[idxnan][['st','ed']].values]
        self._ne2ecov = UT.df2dict(ex, 'name', 'ecov')
        
    def write(self, cmax=9):
        pre = self.dstpre+'.{0}_{1}_{2}'.format(self.chrom,self.st,self.ed)
        # 1) exon, junctions, allpaths => csv (no header <= to concatenate bundles)
        ecols = EXDFCOLS #['chr','st','ed','strand','name','kind','ecov']
        UT.write_pandas(self.exdf[ecols], pre+'.exdf.txt.gz', '')
        scols = SJDFCOLS #['chr','st','ed','strand','name','kind','tcnt'  ]#,'donor','acceptor','dp','ap']
        UT.write_pandas(self.sjdf[scols], pre+'.sjdf.txt.gz', '')
        pcols = PATHCOLS #['chr','st','ed','name','strand','tst','ted','tcov0','tcov1','tcov']
        if len(self.tpaths)>0:
            UT.write_pandas(self.tpaths[pcols], pre+'.paths.txt.gz', '')
            # 2) unused sjpaths => bed12
            GGB.write_bed(self.unusedsj, pre+'.unused.sjpath.bed.gz', ncols=12)
            # 3) allpaths => gtf or bed12 tcov => sc2 rgb color
            self.bed12 = path2bed12(self.tpaths.copy(), cmax, 'tcov')
            GGB.write_bed(self.bed12, pre+'.paths.bed.gz',ncols=12)
            self.tspan = path2tspan(self.tpaths.copy(), cmax, 'tcov0')
            GGB.write_bed(self.tspan, pre+'.tspans.bed.gz',ncols=12)
        
    def process(self):
        self.load_and_filter_sjpaths()
        if len(self.sjpaths)==0:
            return None
        # self.logdebug('finding exons...')
        self.find_exons()

        # self.logdebug('finding 53 positions...')
        self.find_53()

        # self.logdebug('finding paths...')
        self.find_all_paths()
        self.decompose_paths()

        # self.logdebug('finding 53 exons...')
        self.find_53exons()

        # self.logdebug('making dataframes...')
        self.organize()

        # self.logdebug('calculating covrages...')
        self.calculate_ecovs()
        self.calculate_scovs()
        self.estimate_abundance()

        #self.extract_se_candidates()
        # self.logdebug('writing results...')
        self.write()
        self.loginfo('finished assembling, {0} paths found'.format(len(self.tpaths)))
        if len(self.tpaths)>0:
            return self.bname
        return None
        
    #def extract_se_candidates(self):
    #    pass
        
    def draw_covs(self, st, ed, strand, win=10000, ax=None):
        offset = self.st
        s0 = st-win-offset
        e0 = ed+win-offset
        sjap0 = self.arrs['sj'][strand][s0:e0]
        exap0 = self.arrs['ex'][strand][s0:e0]
        sjap1 = self.filled[strand][s0:e0]
        if ax is None:
            fig,ax = P.subplots(1,1,figsize=(15,3))
        y0 = N.log2(sjap1+1)
        ax.plot(y0, 'r-', alpha=0.8)
        #ax.plot(N.log2(sjap0+1), 'r--')
        ipx = set(N.nonzero(exap0>0)[0])
        n = len(exap0)
        ipx.update([x+1 for x in ipx if x<n-1])
        ipx.update([x-1 for x in ipx if x>1])
        ipx = sorted(ipx)
        ax.fill_between(ipx, 0, N.log2(exap0[ipx]+1), facecolor='m', alpha=0.3)
        #axr.axhline(0,linestyle='--',color='grey')
        #x = N.log2(sjap1+1)
        #if strand=='+':
        #    xd = x[1:]-x[:-1]
        #else:
        #    xd = x[:-1]-x[1:]
        #ax.plot(xd, 'b-')
        # gspan
        gspan = self.gspans[strand]
        h0 = 15
        for i, (s1,e1) in enumerate(gspan):
            if (e1>s0)&(s1<e0):
                gx1,gx2 = s1-s0,e1-s0
                ax.plot([gx1,gx2],[h0+2,h0+2], 'c')
                gx0 = max(min((gx1+gx2)/2., e0-s0), 0)
                ax.text(gx0, h0-2, '{0}'.format(i))
        # 53
        e53p = self.e53pos[strand]
        t5 = e53p[(e53p['kind']=='5')&(e53p['pos']>s0)&(e53p['pos']<e0)]
        t3 = e53p[(e53p['kind']=='3')&(e53p['pos']>s0)&(e53p['pos']<e0)]
        i5p = N.array(t5['pos'].values-s0, dtype=N.int64)
        i3p = N.array(t3['pos'].values-s0, dtype=N.int64)
        if len(i5p)>0:
            ax.plot(i5p, y0[i5p], 'm^')
        if len(i3p)>0:
            ax.plot(i3p, y0[i3p], 'mv')
        # exons
        ex = self.exons[strand]
        ex = ex[(ex['ost']>s0)&(ex['oed']<e0)]
        ymid = h0+2
        h = 4
        yrange = (ymid-h/2., h)
        xranges = [(x-s0,y-x) for x,y in ex[['ost','oed']].values]
        cargs = dict(facecolor='k', edgecolor='k')#, linewidth=0.2)
        bbhc = BrokenBarHCollection(xranges, yrange, **cargs)
        ax.add_collection(bbhc)
        
        ax.set_xlim(0,e0-s0)
        ax.set_ylim(-2,h0+6)
        ax.set_yticks([0,5,10])
        ax.set_xticks([])
        txt = '{0}:{1}-{2}:{3}'.format(self.chrom, st-win, ed+win, strand)
        print(txt)
        ax.text(0,h0+4, txt)
        ax.set_frame_on(False)        
        return ax
        

    def draw_path(self, pathdf, st, ed, strand, covfld='tcov',win=10000, ax=None, delta=500):
        offste = self.st
        st0 = st-win
        ed0 = ed+win
        idx = (((pathdf['tst']>=st0)&(pathdf['tst']<=ed0))|\
              ((pathdf['ted']>=st0)&(pathdf['ted']<=ed0)))&\
              (pathdf['strand'].isin(STRS[strand]))&(pathdf['chr']==self.chrom)
        df = pathdf[idx].sort_values(['tst','ted']).copy()
        esiz = 100
        h = 2
        cnt = 0
        cted = 0
        minypos = 0
        lss = {'+':'-','-':'-','.+':'--','.-':'--'}
        cbs = Colors('gray_r',1.,0.)
        cls = {'+':Colors('R',1.,0.),'-':Colors('B',1.,0.),
               '.+':Colors('gray_r',1.,0.),'.-':Colors('gray_r',1.,0.)}
        if covfld not in df.columns:
            df[covfld] = 1.
        df['ltcov'] = N.log2(df[covfld]+2)
        df['tcovn'] = df['ltcov']/df['ltcov'].max()
        for pc, tst, ted, s, tcov in df[['name','tst','ted','strand','tcovn']].values:
            if cted+delta>tst:
                cnt +=1
            else:
                cnt = 0
            cted = max(ted, cted)
            ymid = -cnt*(h+1)
            minypos = min(ymid, minypos)
            cb = cbs.to_rgba(tcov)
            cl = cls[s].to_rgba(tcov)
            ls = lss[s]
            cargs = dict(facecolor=cb, edgecolor=cb)
            ax.plot([tst-st0,ted-st0],[ymid,ymid],ls=ls, color=cl)
            yrange = (ymid-h/2., h)
            tmp = pc.split(',')
            try:# without 5',3' exons
                # pathcode = dpos0|apos1,dpos1|apos2,...dpos(n-1)|aposn
                tst = int(tmp[0].split('|')[0])
                ted = int(tmp[-1].split('|')[1])
                if tst<ted:
                    tmppc = str(tst-esiz)+','+pc+','+str(ted+esiz)
                else:
                    tmppc = str(tst+esiz)+','+pc+','+str(ted-esiz)                
            except:# with 5',3' exons
                # pathcode = 5pos,dpos0|apos1,dpos1|apos2,...dpos(n-1)|aposn,3pos
                tst = int(tmp[1].split('|')[0])
                ted = int(tmp[-2].split('|')[1])
                tmppc = pc
            if tst<ted:
                exons = [[int(x) for x in y.split(',')] for y in tmppc.split('|')]
            else:
                exons = [[int(x) for x in y.split(',')] for y in tmppc.split('|')]
                exons = [x[::-1] for x in exons[::-1]]
            xranges = [(x-st0,y-x) for x,y in exons]
            bbhc = BrokenBarHCollection(xranges, yrange, **cargs)
            ax.add_collection(bbhc)
        ax.set_ylim(minypos-5, 5)
        ax.set_xlim(0,ed0-st0)
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_frame_on(False)
        return ax

def path2bed12(paths, cmax=9, covfld='tcov'):
    bed = paths
    # strand .+,.-
    idxun = bed['strand']=='.-'
    bed.loc[idxun, 'name']= [reversepathcode(x) for x in bed[idxun]['name']]
    bed = bed.groupby('name').first().reset_index()
    idxu = bed['strand'].isin(['.+','.-'])
    bed.loc[idxu, 'strand']='.'
    # #exons, esizes, estarts
    bed['exonsp'] = [[[int(z) for z in y.split(',')] for y in x.split('|')] for x in bed['name']]
    bed['exonsn'] = [[y[::-1] for y in x][::-1] for x in bed['exonsp']]
    idxp = bed['strand']!='-'
    bed.loc[idxp, 'exons'] = bed[idxp]['exonsp']
    bed.loc[~idxp, 'exons'] = bed[~idxp]['exonsn']
    bed['#exons'] = [len(x) for x in bed['exons']]
    estarts = [[str(y[0]-x[0][0]) for y in x] for x in bed['exons']]
    esizes = [[str(y[1]-y[0]) for y in x] for x in bed['exons']]
    bed['esizes'] = [','.join(x)+',' for x in esizes]
    bed['estarts'] = [','.join(x)+',' for x in estarts]
    # sc1, sc2
    bed['ltcov'] = N.log2(bed[covfld]+2)
    # bed['sc1'] = N.ceil(bed['ltcov']*100).astype(int)
    bed['sc1'] = bed[covfld]
    sm = {'+':Colors('R', cmax),
          '-':Colors('B', cmax),
          '.':Colors('G', cmax)}
    bed['sc2'] = [sm[s].RGB(x) for x,s in bed[['ltcov','strand']].values]
    bed.sort_values(['chr','st','ed'], inplace=True)
    return bed

def path2tspan(paths, cmax=9, covfld='tcov0'):
    bed = paths
    # strand .+,.-
    idxun = bed['strand']=='.-'
    bed.loc[idxun, 'name']= [reversepathcode(x) for x in bed[idxun]['name']]
    bed = bed.groupby('name').first().reset_index()
    idxu = bed['strand'].isin(['.+','.-'])
    bed.loc[idxu, 'strand']='.'
    # #exons, esizes, estarts
    bed['exonsp'] = [[[int(z) for z in y.split(',')] for y in x.split('|')] for x in bed['name']]
    bed['exonsn'] = [[y[::-1] for y in x][::-1] for x in bed['exonsp']]
    idxp = bed['strand']!='-'
    bed.loc[idxp, 'exons'] = bed[idxp]['exonsp']
    bed.loc[~idxp, 'exons'] = bed[~idxp]['exonsn']

    # group same (tst,ted) ##############
    bg = bed.groupby(['tst','ted'])
    bedg = bg.first()
    bedg['exons'] = bg['exons'].apply(lambda g: sorted(set([tuple(x) for y in g for x in y])))
    bedg['st'] = bg['st'].min()
    bedg['ed'] = bg['ed'].max()
    bed = bedg.reset_index()
    ######################################

    bed['#exons'] = [len(x) for x in bed['exons']]
    estarts = [[str(y[0]-x[0][0]) for y in x] for x in bed['exons']]
    esizes = [[str(y[1]-y[0]) for y in x] for x in bed['exons']]
    bed['esizes'] = [','.join(x)+',' for x in esizes]
    bed['estarts'] = [','.join(x)+',' for x in estarts]
    # sc1, sc2
    bed['ltcov'] = N.log2(bed[covfld]+2)
    # bed['sc1'] = N.ceil(bed['ltcov']*100).astype(int)
    bed['sc1'] = bed[covfld]
    sm = {'+':Colors('R', cmax),
          '-':Colors('B', cmax),
          '.':Colors('G', cmax)}
    bed['sc2'] = [sm[s].RGB(x) for x,s in bed[['ltcov','strand']].values]
    bed.sort_values(['chr','st','ed'], inplace=True)
    return bed    

def sjpaths2tspan(sjpaths, cmax=9, strip53=False, sc2color=True):
    # bed12
    bed = sjpaths
    # #exons, esizes, estarts
    bed['exonsp'] = [[[int(z) for z in y.split(',')] for y in x.split('|')] for x in bed['name']]
    bed['exonsn'] = [[y[::-1] for y in x][::-1] for x in bed['exonsp']]
    idxp = bed['strand']!='-'
    bed.loc[idxp, 'exons'] = bed[idxp]['exonsp']
    bed.loc[~idxp, 'exons'] = bed[~idxp]['exonsn']

    bg = sjpaths.groupby(['tst','ted'])
    bedg = bg.first()
    bedg['exons'] = bg['exons'].apply(lambda g: sorted(set([tuple(x) for y in g for x in y])))
    bedg['st'] = bg['st'].min()
    bedg['ed'] = bg['ed'].max()
    bedg['sc1'] = bg['sc1'].sum()
    bed = bedg.reset_index()
    
    bed['#exons'] = [len(x) for x in bed['exons']]
    estarts = [[str(y[0]-x[0][0]) for y in x] for x in bed['exons']]
    esizes = [[str(y[1]-y[0]) for y in x] for x in bed['exons']]
    bed['esizes'] = [','.join(x)+',' for x in esizes]
    bed['estarts'] = [','.join(x)+',' for x in estarts]
    if sc2color:
        bed['ltcov'] = N.log2(bed['sc1']+2)
        sm = {'+':Colors('R', cmax),
              '-':Colors('B', cmax),
              '.':Colors('G', cmax)}
        bed['sc2'] = [sm[s].RGB(x) for x,s in bed[['ltcov','strand']].values]
    else:
        bed['sc2'] = bed['sc1']
    if strip53:
        bed['name'] = [','.join(x.split(',')[1:-1]) for x in bed['name']]
    bed.sort_values(['chr','st','ed'], inplace=True)
    return bed


####### Bundle Finder ###############################################################
    
def find_gaps(bwpre, chrom, csize, gsizeth=5e5, minbundlesize=10e6, sjbwpre=None, sjth=0):
    sjexbw = SjExBigWigs(bwpre, sjbwpre, mixunstranded=False)
    sts = []
    eds = []
    bsize = 2*minbundlesize
    bnum = int(N.ceil(csize/float(bsize)))
    with sjexbw:
        for i in range(bnum):
            st = i*bsize
            ed = min((i+1)*bsize, csize)
            arr = sjexbw.bws['sj']['+'].get(chrom,st,ed)
            arr += sjexbw.bws['sj']['-'].get(chrom,st,ed)
            idx = N.nonzero(arr<=sjth)[0]
            if len(idx)==0:
                continue
            dif = idx[1:]-idx[:-1] # if continuous dif==1
            idx2 = N.nonzero(dif>1)[0] # non contiguous    
            # gap start idx2[x]+1 ==> gap end idx2[x+1]
            gsize = idx[idx2[1:]]-idx[idx2[:-1]+1]
            idx3 = N.nonzero(gsize>gsizeth)[0]
            gst = idx[idx2[idx3]+1]+st
            ged = idx[idx2[idx3+1]]+st
            sts += list(gst)
            eds += list(ged)
    return sts,eds

def find_bundles(bwpre, genome, dstpre, chrom=None, mingap=5e5, minbundlesize=10e6, sjbwpre=None, sjth=0):
    bundles = []
    if chrom is None:
        chroms = UT.chroms(genome) # whole genome
        fpath = dstpre+'.bundles.txt.gz'
        if os.path.exists(fpath):
            df = UT.read_pandas(fpath)
            return df[['chr','st','ed']].values
    else:
        chroms = [chrom]
        fpath = dstpre+'.{0}.bundles.txt.gz'.format(chrom)
        if os.path.exists(fpath):
            df = UT.read_pandas(fpath)
            return df[['chr','st','ed']].values
    chromsizes = UT.df2dict(UT.chromdf(genome), 'chr', 'size')
    for chrom in chroms:
        print('checking {0}...'.format(chrom))
        csize = chromsizes[chrom]
        sts,eds = find_gaps(bwpre, chrom, csize, mingap, minbundlesize,sjbwpre,sjth)
        st = 0
        if len(sts)==0:
            bundles.append((chrom,0,csize))
        else:
            for gs,ge in zip(sts,eds):
                mid = int((gs+ge)/2.)
                if mid-st>minbundlesize:
                    bundles.append((chrom,st,mid))
                    st = mid
            if ge<csize:
                bundles.append((chrom,st,csize))
    df = PD.DataFrame(bundles, columns=['chr','st','ed'])
    UT.write_pandas(df, fpath, 'h')
    return bundles


########## convenience 

def drawspan(la,st,ed,strand,win=10000, figsize=(15,6), df2=None, df3=None, delta=500, 
    df2cov='sc2', df3cov='tcov', maxdisp=None):
    fig, axr = P.subplots(3,1,figsize=figsize,sharex=True)
    P.subplots_adjust(hspace=0)
    la.draw_covs(st,ed,strand, ax=axr[0], win=win)
    if df2 is not None:
        if maxdisp is not None:
            df2 = df2.sort_values(df2cov,ascending=False).iloc[:maxdisp]
        la.draw_path(df2, st,ed,strand, ax=axr[1], covfld=df2cov, win=win, delta=delta)
    if df3 is not None:
        if maxdisp is not None:
            df3 = df3.sort_values(df2cov,ascending=False).iloc[:maxdisp]
        la.draw_path(df3, st,ed,strand, ax=axr[2], win=win, delta=delta, covfld=df3cov)

def drawspan2(la,st,ed,win=10000, figsize=(15,6), df2=None, df3=None, delta=500,
    df2cov='sc2', df3cov='tcov', maxdisp=None):
    o = la.st
    fig, axr = P.subplots(6,1,figsize=figsize,sharex=True)
    P.subplots_adjust(hspace=0)
    strand = '+'
    la.draw_covs(st,ed,strand, ax=axr[0], win=win)
    if df2 is not None:
        if maxdisp is not None:
            df2 = df2.sort_values(df2cov,ascending=False).iloc[:maxdisp]
        la.draw_path(df2, st,ed,strand, ax=axr[1], covfld=df2cov, win=win, delta=delta)
    if df3 is not None:
        if maxdisp is not None:
            df3 = df3.sort_values(df2cov,ascending=False).iloc[:maxdisp]
        la.draw_path(df3, st,ed,strand, ax=axr[2], win=win, delta=delta, covfld=df3cov)
    strand = '-'
    la.draw_covs(st,ed,strand, ax=axr[3], win=win)
    if df2 is not None:
        if maxdisp is not None:
            df2 = df2.sort_values(df2cov,ascending=False).iloc[:maxdisp]
        la.draw_path(df2, st,ed,strand, ax=axr[4], covfld=df2cov, win=win, delta=delta)
    if df3 is not None:
        if maxdisp is not None:
            df3 = df3.sort_values(df2cov,ascending=False).iloc[:maxdisp]
        la.draw_path(df3, st,ed,strand, ax=axr[5], win=win, delta=delta, covfld=df3cov)



######### Chrom Assembler


def bundle_assembler(bwpre, chrom, st, ed, dstpre, upperpathnum, sjbwpre=None, refcode='gen9'):
    bname = bundle2bname((chrom,st,ed))
    bsuf = '.{0}_{1}_{2}'.format(chrom,st,ed)
    csuf = '.{0}'.format(chrom)
    LOG.info('assembling bunle {0}'.format(bname))
    sufs = ['.exdf.txt.gz',
            '.sjdf.txt.gz',
            '.paths.txt.gz',
            '.paths.bed.gz',
            '.tspans.bed.gz',
            '.unused.sjpath.bed.gz']
    done = []
    for x in sufs:
        done.append(os.path.exists(dstpre+bsuf+x) | \
                    os.path.exists(dstpre+csuf+x) | \
                    os.path.exists(dstpre+x) )
    if all(done):
        return bname
    la = LocalAssembler(bwpre, chrom, st, ed, dstpre, upperpathnum=upperpathnum, sjbwpre=sjbwpre, refcode=refcode)
    return la.process()

def bname2bundle(bname):
    # bname = 'chrom:st-ed'
    chrom, sted = bname.split(':')
    st,ed = [int(x) for x in sted.split('-')]
    return chrom,st,ed

def bundle2bname(b):
    return '{0}:{1}-{2}'.format(*b)

def chrom_assembler(bwpre, dstpre, genome, chrom, mingap=5e5, minbundlesize=10e6, np=2):
    bundles = find_bundles(bwpre, genome, dstpre, chrom, mingap, minbundlesize)
    print('{1}: #bundles: {0}'.format(len(bundles), chrom))
    server = TQ.Server(np=np)
    with server:
        for c,st,ed in bundles:
            tname = 'bundle_assembler.{0}:{1}-{2}'.format(c,st,ed)
            args = (bwpre, c, st, ed, dstpre)
            task = TQ.Task(tname, bundle_assembler, args)
            server.add_task(task)
    rslts = {}
    n = len(bundles)
    for i in range(n):
        try:
            name,ans = server.get_result(True, 10) # block, r: (value, name)
            bname = name.split('.')[1]
            rslts[bname2bundle(bname)] = ans
        except multiprocessing.Queue.Empty:
            pass
    bundles1 = []
    for x in bundles:
        if rslts.get(x,None) is None:
            print('no results from bundle {0}'.format(bundle2bname(x)))
        else:
            bundles1.append(x)
    concatenate_bundles(bundles1, chrom, dstpre)
    return '{0}.{1}'.format(dstpre,chrom)


def concatenate_bundles(bundles, bundlestatus, chrom, dstpre):
    # concat results
    sufs = ['exdf.txt.gz', 
           'sjdf.txt.gz',
           'paths.txt.gz',
           'paths.bed.gz',
           'tspans.bed.gz',
           'unused.sjpath.bed.gz']
    files = []
    for suf in sufs:
        dstpath = '{0}.{1}.{2}'.format(dstpre, chrom, suf)
        dstpath2 = '{0}.{1}'.format(dstpre, suf)
        if not os.path.exists(dstpath2):
            if not os.path.exists(dstpath):
                with open(dstpath, 'wb') as dst:
                    for chrom, st, ed in bundles:
                        bname = bundle2bname((chrom,st,ed))
                        srcpath = '{0}.{1}_{2}_{3}.{4}'.format(dstpre, chrom, st, ed, suf)
                        if not os.path.exists(srcpath):
                            if bundlestatus[bname] is None:
                                continue
                            else:
                                raise RuntimeError('{0} does not exists'.format(srcpath))
                        files.append(srcpath)
                        with open(srcpath, 'rb') as src:
                            shutil.copyfileobj(src, dst)
        else:
            files+=['{0}.{1}_{2}_{3}.{4}'.format(dstpre, chrom, st, ed, suf) for chrom,st,ed in bundles]
    # cleanup
    for f in files:
        if os.path.exists(f):
            os.unlink(f)

def concatenate_chroms(chroms, dstpre):
    # concat results
    sufs = ['exdf.txt.gz', 
           'sjdf.txt.gz',
           'paths.txt.gz',
           'paths.bed.gz',
           'tspans.bed.gz',
           'unused.sjpath.bed.gz']
    files = []
    for suf in sufs:
        dstpath = '{0}.{1}'.format(dstpre, suf)
        if not os.path.exists(dstpath):
            with open(dstpath, 'wb') as dst:
                for chrom in chroms:
                    srcpath = '{0}.{1}.{2}'.format(dstpre, chrom, suf)
                    if os.path.exists(srcpath):
                        files.append(srcpath)
                        with open(srcpath, 'rb') as src:
                            shutil.copyfileobj(src, dst)
        else:
            files+=['{0}.{1}.{2}'.format(dstpre, chrom, suf) for chrom in chroms]
    # cleanup
    for f in files:
        if os.path.exists(f):
            os.unlink(f)

def write_stats(dstpre, seinfo):
    dic = {}
    dic.update(seinfo)
    exdf = UT.read_pandas(dstpre+'.exdf.txt.gz', names=EXDFCOLS)
    dic['num_me_exons'] = len(exdf)
    sjdf = UT.read_pandas(dstpre+'.sjdf.txt.gz', names=SJDFCOLS)
    dic['num_junctions'] = len(sjdf)
    paths = UT.read_pandas(dstpre+'.paths.txt.gz', names=PATHCOLS)
    dic['num_paths'] = len(paths)
    unused = GGB.read_bed(dstpre+'.unused.sjpath.bed.gz')
    dic['num_unused_junctions'] = len(unused)
    fname = dstpre+'.stats.txt'
    name = dstpre.split('/')[-1]
    df = PD.DataFrame(dic, index=[name])
    LOG.info('{0}:{1}'.format(name, dic))
    UT.write_pandas(df, fname, 'ih')
    


def sample_assembler(bwpre, dstpre, genome, mingap=5e5, minbundlesize=10e6, np0=2, np1=2, chroms=None):
    server = TQ.Server(np=np0)
    if chroms is None:
        chroms = UT.chroms(genome)
    with server:
        for chrom in chroms:
            tname = 'chrom_assembler.{0}'.format(chrom)
            args = (bwpre, dstpre, genome, chrom, mingap, minbundlesize, np1)
            task = TQ.Task(tname, chrom_assembler, args)
            server.add_task(task)
    rslts = {}
    n = len(chroms)
    for i in range(n):
        try:
            name,ans = server.get_result(True, 10) # block, r: (value, name)
            c = name.split('.')[1]
            rslts[c] = ans
        except multiprocessing.Queue.Empty:
            pass
    chroms1 = []
    for x in chroms:
        if rslts.get(x,None) is None:
            print('no results from {0}'.format(x))
        else:
            chroms1.append(x)    
    concatenate_chroms(chroms1, dstpre)
    return dstpre


def find_SE_chrom(bwpre, dstpre, genome, chrom, exstrand='+', minsizeth=200):
    # find SE candidates and calculate ecovs
    try:
        exdf = UT.read_pandas(dstpre+'.{0}.exdf.txt.gz'.format(chrom), names=EXDFCOLS)
    except:
        exdf = UT.read_pandas(dstpre+'.exdf.txt.gz', names=EXDFCOLS)
    exdf = exdf[exdf['chr']==chrom]
    sjexbw = SjExBigWigs(bwpre)
    chromdf = UT.chromdf(genome).set_index('chr')
    csize = chromdf.ix[chrom]['size']
    with sjexbw:
        exa = sjexbw.bws['ex'][exstrand].get(chrom,0,csize)
    # mask existing exons
    for st,ed in exdf[['st','ed']].values:
        exa[st:ed] = 0
    idx = N.nonzero(exa>0)[0]
    # find continuous segments
    dif = idx[1:]-idx[:-1] # if continuous dif==1
    idx2 = N.nonzero(dif>1)[0] # non contiguous   
    idxst = N.array([idx[0]]+list(idx[idx2+1]))
    idxed = N.array(list(idx[idx2])+[idx[-1]])
    # gap start idx2[x]+1 ==> gap end idx2[x+1]
    gsize = idxed - idxst + 1
    idx3 = N.nonzero(gsize>minsizeth)[0]
    st = idxst[idx3]
    ed = st + gsize[idx3]
    df = PD.DataFrame({'st':st, 'ed':ed}, index=N.arange(len(st)))
    df['ecov'] = [N.mean(exa[x:y]) for x,y in df[['st','ed']].values]
    df['len'] = df['ed']-df['st']
    df['chr'] = chrom
    dstpath = dstpre+'.se.{0}.{1}.txt.gz'.format(chrom,exstrand)
    UT.write_pandas(df[['chr','st','ed','ecov','len']], dstpath, '')
    return dstpath

def find_SE(dstpre, chroms, exstrand='+', sestrand='.', 
    mincovth=5, minsizeth=200, minsep=1000, cmax=9, mergedist=200):
    # concatenate
    dstpath = dstpre+'.se0.txt.gz'
    if not os.path.exists(dstpath):
        with open(dstpath,'wb')  as dst:
            for chrom in chroms:
                srcpath = dstpre+'.se.{0}.{1}.txt.gz'.format(chrom,exstrand)
                with open(srcpath,'rb') as src:
                    shutil.copyfileobj(src,dst)
                os.unlink(srcpath)
    # concatenate 
    if not os.path.exists(dstpre+'.exdf.txt.gz'):
        concatenate_chroms(chroms, dstpre)
    secols = ['chr','st','ed','ecov','len']
    sedf = UT.read_pandas(dstpath, names=secols)
    exdf = UT.read_pandas(dstpre+'.exdf.txt.gz', names=EXDFCOLS) 
    # th = find_threshold(exdf['ecov'].values, sedf['ecov'].values, mincovth, dstpre)
    paths = UT.read_pandas(dstpre+'.paths.txt.gz', names=PATHCOLS)
    th = find_threshold(paths['tcov'].values, sedf['ecov'].values, mincovth, dstpre)
    se0 = sedf[(sedf['ecov']>th)&(sedf['len']>minsizeth)].copy()  # use FPR 1%

    LOG.info('SE covth={0:.2f}, len(se0)={1}'.format(th, len(se0)))
    # se0['strand'] = sestrand
    # se0['name'] = [_pc(st,ed,sestrand,',') for st,ed in se0[['st','ed']].values ]
    # se0['kind'] = 's'
    # c = dstpre+'.se1.txt.gz'
    # UT.write_pandas(se0, c, 'h')
    # return (c, th, len(sedf), len(se0))
    # check irets in bundle workers?
    # find intervals at least minsep separated from exons
    a = dstpre+'.se0.bed'
    b = dstpre+'.espans.bed'
    # cols = ['chr','st','ed','ecov','strand']
    se0 = se0[['chr','st','ed','ecov']].sort_values(['chr','st','ed'])
    UT.write_pandas(se0, a, '')
    exdf['st0'] = exdf['st']-minsep
    exdf['ed0'] = exdf['ed']+minsep
    UT.write_pandas(exdf[['chr','st0','ed0']],b,'')
    c0 = dstpre+'.sedf0.txt.gz'
    BT.bedtoolintersect(a,b,c0,v=True)
    # merge nearby 
    c = dstpre+'.sedf.txt.gz'
    BT.bedtoolmerge(c0,c,d=mergedist,c=4,o='mean')
    se1a = UT.read_pandas(c0,names=['chr','st','ed','ecov'])
    se1 = UT.read_pandas(c, names=['chr','st','ed','ecov'])
    se1['strand'] = sestrand
    se1['name'] = [_pc(st,ed,sestrand,',') for st,ed in se1[['st','ed']].values ]
    se1['kind'] = 's'
    UT.write_pandas(se1[EXDFCOLS], c, '')

    cbed = dstpre+'.se.bed.gz'
    GGB.write_bed(se1, cbed, ncols=3)
    LOG.info('#SE = {0} (before merge {1})'.format(len(se1), len(se1a)))
    os.unlink(a)
    os.unlink(b)
    os.unlink(c0)
    # merge exdf & se ==> update .exdf.txt.gz?
    # update .paths.txt.gz, .paths.bed.gz?
    bed = GGB.read_bed(dstpre+'.paths.bed.gz') # chr,st,ed,name,sc1,strand,sc2,tst,ted,#exons,esizes,estarts

    # BED12
    se1['ltcov'] = N.log2(se1['ecov']+2)
    se1['sc1'] = N.ceil(se1['ltcov']*100).astype(int)
    sm = {'+':Colors('R', cmax),
          '-':Colors('B', cmax),
          '.':Colors('C', cmax)}
    se1['sc2'] = [sm[s].RGB(x) for x,s in se1[['ltcov','strand']].values]
    se1['tst'] = se1['st']
    se1['ted'] = se1['ed']
    se1['#exons'] = 1
    se1['len'] = se1['ed']-se1['st']
    se1['esizes'] = se1['len'].astype(str)+','
    se1['estarts'] = '0,'
    bed1 = PD.concat([bed, se1], ignore_index=True)
    GGB.write_bed(bed1, dstpre+'.paths.withse.bed.gz', ncols=12)
    dic = dict(
        secovth=th, 
        num_se_candidates=len(sedf), 
        num_se_by_covth_and_size=len(se0), 
        num_se_not_near_exons=len(se1a),
        num_se_merge_nearby=len(se1))
    return dic

def find_threshold(x0,x1,minth,dstpre,fdrth=0.5, fprth=0.01):
    x0 = x0[(~N.isnan(x0))&(x0>0)]  # why exdf contains NaN?
    x1 = x1[(~N.isnan(x1))&(x1>0)]
    x0 = N.log2(x0+1)
    x1 = N.log2(x1+1)
    xmax = min(x0.max(), x1.max())
    xmin = max(x0.min(), x1.min())
    delta = (xmax-xmin)/100.
    bins=N.arange(xmin,xmax,delta)
    h0,b0 = N.histogram(x0, bins=bins)
    h1,b1 = N.histogram(x1, bins=bins)
    # def find_best_xmid(h0,h1):
    deltas = []
    for i in range(len(h0)-25):
        scale = float(N.sum(h1[i:]))/N.sum(h0[i:])
        h0s = h0*scale
        delta = N.mean(N.abs(h0s[i:]-h1[i:]))
        deltas.append([delta, i, h0s, scale])
    delta,posi,h0s,scale =sorted(deltas)[0]
    cntf = h0s
    cnto = h1
    cp = N.sum(cntf)   # total positive
    cn = N.sum(cnto)-cp # total negative (observed - positive)
    if cn<0: # due to noise when almost perfect
        th_fpr = minth
        th_fdr = minth
    else:
        fn = N.cumsum(cntf) # false negative (condition positive but < th)
        tp = cp - fn # true positive (condition positive and >= th)
        tn = N.cumsum(cnto) - fn  # true negative
        tn[tn<0]=N.nan
        fp = cn - tn
        fp[fp<0]=N.nan
        #tpr = tp/cp
        fpr = fp/cn
        p = N.sum(cnto)-N.cumsum(cnto) # declared positive
        fdr = fp/p
        fpr[N.isnan(fpr)]=0
        idx = N.nonzero(fdr<=fdrth)[0]
        if len(idx)==0: # not found
            th_fdr = minth
        else:
            th_fdr = 2**(bins[N.min(idx)])-1
        idx0 = N.nonzero(fpr<=fprth)[0]
        if len(idx0)==0: # not found
            th_fpr = minth
        else:
            th_fpr = 2**(bins[N.min(idx0)])-1
    fname = dstpre+'.secovth.pdf'
    title = dstpre.split('/')[-1]
    plot_se_th(b0,h0s,h1,th_fpr,th_fdr,title,fname)
    return th_fpr


def plot_se_th(b0,h0s,h1,th0,th,title,fname=None):
    fig,ax = P.subplots(1,1)
    w = b0[1]-b0[0]
    ax.bar(b0[:-1], h1, width=w, alpha=0.9,color='c', label='single-exon', lw=0)
    ax.bar(b0[:-1], h0s, width=w, alpha=0.5, color='r', label='multi-exon', lw=0)
    ax.set_yscale('log')
    ax.axvline(N.log2(th+1), color='r', linestyle='--', label='FDR 50%')
    ax.axvline(N.log2(th0+1), color='b', linestyle='--', label='FPR 1%')
    ax.set_xlabel('log2(cov+1)')
    ax.set_ylabel('count')
    ax.set_title(title)
    ax.legend()
    if fname is not None:
        fig.savefig(fname)
    

def find_threshold0(x0,x1,minth):
    x0 = x0[~N.isnan(x0)]  # why exdf contains NaN?
    x1 = x1[~N.isnan(x1)]
    x0 = N.log2(x0+1)
    x1 = N.log2(x1+1)
    xmax = min(x0.max(), x1.max())
    xmin = max(x0.min(), x1.min())
    delta = (xmax-xmin)/25.
    bins=N.arange(xmin,xmax,delta)
    h0,b0 = N.histogram(x0, bins=bins)
    h1,b1 = N.histogram(x1, bins=bins)
    xmid = 0.6*(xmax-xmin)+xmin
    scale = float(N.sum(x1>xmid))/N.sum(x0>xmid)
    h0s = scale*h0
    d0 = N.abs(h1-h0s)/h0s
    d0s = smooth(d0,5)
    d0sd = d0s[1:]-d0s[:-1]
    d0sdd = d0sd[1:]-d0sd[:-1]
    for i in range(4,len(d0sdd)):
        if d0sdd[i-1]*d0sdd[i]<0:
            break
    th = b0[1:][i]
    if th<N.log2(minth+1):
        return minth
    if th>0.8*xmax:
        LOG.warning('find_threshold: threshold too large {0} returning 0.8*xmax {1}'.format(th,0.7*xmid))
        return 2**(0.8*xmax)-1
    return 2**th-1


def smooth( v, wsize):
    swin = N.ones(wsize)
    v0 = N.concatenate([swin*v[0], v, swin*v[-1]])
    return N.convolve(v0, swin/float(wsize), 'same')[wsize:-wsize]

LAPARAMS = dict(
    uth=1, 
    mth=3, 
    sjratioth=2e-3, 
    usjratioth=1e-2,
    covfactor=0.05, 
    covth=0.1,
    upperpathnum=5000,
)
SEPARAMS = dict(
    semincovth=5,
    seminsizeth=50,
)
BUNDLEPARAMS = dict(
    mingap=5e5, 
    minbundlesize=30e6, 
)

class SampleAssembler(object):

    def __init__(self, bwpre, dstpre, genome, 
                sjbwpre=None,
                refcode='gen9',
                sjth=0,
                mingap=1e5, 
                minbundlesize=20e6, 
                np=4, 
                chroms=None, 
                maxwaittime=600,
                semincovth=5,
                seminsizeth=50,
                upperpathnum=3000,
                ):
        self.bwpre = bwpre
        self.sjbwpre = sjbwpre
        self.refcode = refcode
        self.sjth = sjth
        self.dstpre = dstpre
        self.genome = genome
        self.mingap = mingap
        self.minbundlesize = minbundlesize
        self.np = np
        self.chroms = chroms
        if self.chroms is None:
            self.chroms = UT.chroms(genome)
        self.semincovth = semincovth
        self.seminsizeth = seminsizeth
        self.upperpathnum = upperpathnum
        self.sestrand = '+'
        self.maxwaittime = maxwaittime # if worker doesn't return within this time limit something is wrong

    def run(self):
        self.server = server = TQ.Server(np=self.np)
        self.bundles = bundles = {} # chr => [(chr,st,ed),...]
        self.bundlestatus = bundlestatus = {} # chrom => bundle (chr,st,ed) => done status
        self.chromstatus = chromstatus = {}
        self.find_se_chrom_status = find_se_chrom_status = {}

        with server:
            for chrom in self.chroms:
                tname = 'find_bundle.{0}'.format(chrom)
                args = (self.bwpre, self.genome, self.dstpre, chrom, self.mingap, 
                        self.minbundlesize, self.sjbwpre, self.sjth)
                task = TQ.Task(tname,find_bundles, args)
                server.add_task(task)
            while server.check_error(self.maxwaittime): # loop
                try:
                    name, rslt = server.get_result(timeout=5) # block until result come in
                except TQ.Empty:
                    name, rslt = None, None
                if name is not None:
                    if name.startswith('find_bundle.'):
                        print('{0}:{1}'.format(name, len(rslt)))
                        chrom = name.split('.')[1]
                        bundles[chrom] = rslt
                        for c,st,ed in rslt:
                            # print('put task##bundle_assembler {0}:{1}-{2}'.format(chrom,st,ed))
                            tname = 'bundle_assembler.{0}:{1}-{2}'.format(c,st,ed)
                            # bwpre, chrom, st, ed, dstpre, upperpathnum
                            args = (self.bwpre, c, st, ed, self.dstpre, self.upperpathnum, self.sjbwpre, self.refcode)
                            task = TQ.Task(tname, bundle_assembler, args)
                            server.add_task(task)
                    if name.startswith('bundle_assembler.'):
                        bname = name.split('.')[1]
                        chrom = bname.split(':')[0]
                        bundlestatus.setdefault(chrom,{})[bname] = rslt
                        if len(bundlestatus[chrom])==len(bundles[chrom]): # all done
                            # print('put task##concatenate_bundles {0}'.format(chrom))
                            tname = 'concatenate_bundles.{0}'.format(chrom)
                            args = (bundles[chrom], bundlestatus[chrom], chrom, self.dstpre)
                            task = TQ.Task(tname, concatenate_bundles, args)
                            server.add_task(task)
                    if name.startswith('concatenate_bundles.'):
                        chrom = name.split('.')[1]
                        chromstatus[chrom] = rslt
                        # start SE finder for chrom
                        # print('put task##find_SE_chrom {0}'.format(chrom))
                        tname = 'find_SE_chrom.{0}'.format(chrom)
                        # bwpre, dstpre, genome, chrom, exstrand='+', minsizeth
                        args = (self.bwpre, self.dstpre, self.genome, chrom, '+',self.seminsizeth)
                        task = TQ.Task(tname, find_SE_chrom, args)
                        server.add_task(task)
                    if name.startswith('find_SE_chrom.'):
                        chrom = name.split('.')[1]
                        find_se_chrom_status[chrom] = rslt
                        if len(find_se_chrom_status)==len(self.chroms):
                            # print('start SE finder')
                            tname = 'find_SE'
                            # dstpre, chroms, exstrand='+', sestrand='.', mincovth=5, minsizeth
                            args = (self.dstpre, self.chroms, '+', '.', self.semincovth, self.seminsizeth)
                            task = TQ.Task(tname, find_SE, args)
                            server.add_task(task)
                    if name== 'find_SE':
                        tname = 'write_stats'
                        args = (self.dstpre, rslt)
                        task = TQ.Task(tname, write_stats, args)
                        server.add_task(task)
                    if name=='write_stats':
                        break
                        # tname = 'concatenate_chroms'
                        # args = (self.chroms, self.dstpre)
                        # task = TQ.Task(tname, concatenate_chroms, args)
                        # server.add_task(task)
                    # if name== 'concatenate_chroms':
                    #     break
            print('Exit Loop')
        print('Done')
