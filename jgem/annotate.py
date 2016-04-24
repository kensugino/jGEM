"""Copyright (c) 2015-2016 Ken Sugino

.. module:: compare
    :synopsis: compare to a reference and annotate known genes

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
# system imports
import gzip
import os
import subprocess
from collections import Counter
from operator import iadd

# 3rd party libraries
import pandas as PD
import numpy as N
import matplotlib.pyplot as P

# library imports
from jgem import utils as UT
from jgem import gtfgffbed as GGB
from jgem import bedtools as BT
from jgem import bigwig as BW


class Collector(object):
    """
    calculate 
    1. gcov
    2. ecov
    3. jcnt
    for all samples
    """
    def __init__(self, code, snames=None):
        if snames is None:
            snames = SI.sortdf(SI.df1)
        self.paths = paths = MD.paths[code]
        self.expath = paths['ex']
        self.sjpath = paths['sj']
        self.snames = snames # all sample
        self.code = code
                
        #         self.allsjcntname = MD.MRGDIR+'allsamples.{0}.sjcnt.txt.gz'.format(code)
        #         self.allecovname = MD.MRGDIR+'allsamples.{0}.ecov.txt.gz'.format(code)
        #         self.allgcov1name = MD.MRGDIR+'allsamples.{0}.gcov1.txt.gz'.format(code)
        #         self.allgcov2name = MD.MRGDIR+'allsamples.{0}.gcov2.txt.gz'.format(code)
        #         self.allecovname = MD.MRGDIR+'allsamples.{0}.uecov.txt.gz'.format(code)
        #         self.allgcov1name = MD.MRGDIR+'allsamples.{0}.ugcov1.txt.gz'.format(code)
        #         self.allgcov2name = MD.MRGDIR+'allsamples.{0}.ugcov2.txt.gz'.format(code)
        # get all sample path from MD
        self.aspaths = MD.ASPATHS[code]
        
    def collectsj(self):
        self.msj = msj = UT.read_pandas(self.paths['sj'])
        msj['locus'] = GGB.calc_locus_strand(msj)
        for i, s in enumerate(self.snames):
            print(i,s)
            fn = FN.FileNames(s)
            sj = UT.read_pandas(fn.sjname2())
            sj['locus'] = GGB.calc_locus_strand(sj)
            l2u = dict(UT.izipcols(sj, ['locus','sc1']))
            msj[s] = [l2u.get(x,0) for x in msj['locus']]        
        UT.save_tsv_nidx_whead(msj, self.aspaths['sjcnt'])
        
    def collectecov(self,unique=False):
        self.mex = mex = UT.read_pandas(self.paths['ex'])
        w = 'uecov' if unique else 'ecov'
        for i,s in enumerate(self.snames):
            print(i,s)
            fn = FN.FileNames(s)
            if unique:
                path = os.path.join(fn.sjexdir,fn.sname+'.{0}.unique.ecov.txt.gz'.format(self.code))
            else:
                path = os.path.join(fn.sjexdir,fn.sname+'.{0}.ecov.txt.gz'.format(self.code))
            ecov = UT.read_pandas(path)
            mex[s] = ecov.set_index('eid').ix[mex['_id'].values]['ecov'].values
        UT.save_tsv_nidx_whead(mex, self.aspaths[w])
        
    def collectgcov(self,unique=False):
        mex = UT.read_pandas(self.paths['ex'])
        tmp = mex.groupby('_gidx')
        gcov1 = tmp[['chr']].first().reset_index().sort_values('_gidx')
        #gcov1['st'] = tmp['st'].min()
        #gcov1['ed'] = tmp['ed'].max()
        gcov2 = gcov1.copy()
        gids = gcov1['_gidx'].values
        for i,s in enumerate(self.snames):
            print(i,s)
            fn = FN.FileNames(s)
            if unique:
                path = os.path.join(fn.sjexdir,fn.sname+'.{0}.unique.gcov.txt.gz'.format(self.code))
            else:
                path = os.path.join(fn.sjexdir,fn.sname+'.{0}.gcov.txt.gz'.format(self.code))
            cov = UT.read_pandas(path)
            covi = cov.set_index('_gidx')
            gcov1[s] = covi['gcov'].ix[gids].values
            gcov2[s] = covi['gcov2'].ix[gids].values
        gcov1['len'] = covi['len'].ix[gids].values
        gcov2['len'] = covi['len'].ix[gids].values
        gcov1['st'] = covi['st'].ix[gids].values
        gcov2['st'] = covi['st'].ix[gids].values
        gcov1['ed'] = covi['ed'].ix[gids].values
        gcov2['ed'] = covi['ed'].ix[gids].values
        w = 'ugcov' if unique else 'gcov'
        UT.save_tsv_nidx_whead(gcov1, self.aspaths[w+'1'])
        UT.save_tsv_nidx_whead(gcov2, self.aspaths[w+'2'])
        self.gcov1 = gcov1
        self.gcov2 = gcov2
        
    def collectgcov_length(self,unique=True, length=2000):
        mex = UT.read_pandas(self.paths['ex'])
        tmp = mex.groupby('_gidx')
        gcov1 = tmp[['chr']].first().reset_index().sort_values('_gidx')
        gids = gcov1['_gidx'].values
        for i,s in enumerate(self.snames):
            print(i,s)
            fn = FN.FileNames(s)
            if unique:
                path = os.path.join(fn.sjexdir,fn.sname+'.{0}.unique.{1}.gcov.txt.gz'.format(self.code,length))
            else:
                path = os.path.join(fn.sjexdir,fn.sname+'.{0}.{1}.gcov.txt.gz'.format(self.code,length))
            cov = UT.read_pandas(path)
            covi = cov.set_index('_gidx')
            gcov1[s] = covi['gcov'].ix[gids].values
        gcov1['len'] = covi['len'].ix[gids].values
        gcov1['st'] = covi['st'].ix[gids].values
        gcov1['ed'] = covi['ed'].ix[gids].values
        w = 'ugcov{0}'.format(length) if unique else 'gcov{0}'.format(length)
        UT.save_tsv_nidx_whead(gcov1, self.aspaths[w])
        self.gcov1 = gcov1
        
    def collectgcov1000(self,unique=True):
        mex = UT.read_pandas(self.paths['ex'])
        tmp = mex.groupby('_gidx')
        gcov1 = tmp[['chr']].first().reset_index().sort_values('_gidx')
        gids = gcov1['_gidx'].values
        for i,s in enumerate(self.snames):
            print(i,s)
            fn = FN.FileNames(s)
            if unique:
                path = os.path.join(fn.sjexdir,fn.sname+'.{0}.unique.1000.gcov.txt.gz'.format(self.code))
            else:
                path = os.path.join(fn.sjexdir,fn.sname+'.{0}.1000.gcov.txt.gz'.format(self.code))
            cov = UT.read_pandas(path)
            covi = cov.set_index('_gidx')
            gcov1[s] = covi['gcov'].ix[gids].values
        gcov1['len'] = covi['len'].ix[gids].values
        gcov1['st'] = covi['st'].ix[gids].values
        gcov1['ed'] = covi['ed'].ix[gids].values
        w = 'ugcov1000' if unique else 'gcov1000'
        UT.save_tsv_nidx_whead(gcov1, self.aspaths[w])
        self.gcov1 = gcov1

    def collectall(self):
        self.collectsj()
        self.collectecov()
        self.collectgcov()
        
    def collectunique(self):
        self.collectecov(unique=True)
        self.collectgcov(unique=True)

    def write_wiggletool_cmd(self):
        bws = [FN.FileNames(x).bwname() for x in self.snames]
        cmd = "#!/bin/bash\nwiggletools write /groups/neuroseq/home/suginok/mrgd/allsample.bw mean "
        cmd += ' '.join(bws)
        open(os.path.join(MD.MRGDIR,'make_allsamplebw.sh','w')).write(cmd)
        

class ComparatorPaths(object):
    
    def __init__(self,code,evaldir=None):
        paths = MD.paths[code]
        self.ex = paths['ex']
        self.sj = paths['sj']
        self.code = code
        if evaldir is None:
            evaldir = MD.EVALDIR #'/groups/neuroseq/home/suginok/eval/'
        prefix = os.path.join(evaldir,code)
        self.prefix = prefix
        self.exbed = prefix+'.ex.bed.gz'
        self.gbed = prefix+'.gene.bed.gz'
        
    def exovl(self,code2):
        return self.prefix+'.ex.ovl.{0}.txt.gz'.format(code2)

    def geneovl(self,code2):
        # anti-sense overlap
        return self.prefix+'.gene.ovl.{0}.txt.gz'.format(code2)
    
    def g2g(self,code2):
        return self.prefix+'.g2g.{0}.txt.gz'.format(code2)
    
class Comparator(object):
    """
    compare to ref
    """
    def __init__(self, tgtcode, refcode='gen4', refnamefld='gene_name',
                mergedsj=MD.MRGDIR+'merge3-sj-u0m5.txt.gz'):
        self.tgt = tgtcode
        self.ref = refcode
        self.refnamefld = refnamefld
        self.p1 = p1 = ComparatorPaths(tgtcode)
        self.p2 = p2 = ComparatorPaths(refcode)
        self.mergedsj = mergedsj
        # calc overlap
        c = p1.exovl(p2.code)
        cols = ['chr','st','ed','cat','_id','strand','cov','_gidx','len']
        self.e1 = e1 = UT.read_pandas(p1.ex)
        self.e2 = e2 = UT.read_pandas(p2.ex)
        if 'cov' not in e1.columns:
            if 'ecov' in e1.columns:
                e1['cov'] = e1['ecov']
            else:
                raise ValueException('cov,ecov not in e1 columns')
        if 'cov' not in e2.columns:
            if 'ecov' in e2.columns:
                e2['cov'] = e2['ecov']
            else:
                raise ValueException('cov,ecov not in e1 columns')
        e1['len'] = e1['ed']-e1['st']
        e2['len'] = e2['ed']-e2['st']
        a = UT.save_tsv_nidx_nhead(e1[cols],p1.exbed)
        b = UT.save_tsv_nidx_nhead(e2[cols],p2.exbed)
        c = BT.bedtoolintersect(a,b,c,wao=True)
        ocols = cols + ['b_'+x for x in cols] + ['ovl']
        self.ov = ov = UT.read_pandas(c, names=ocols)

        # gene overlap
        gcols = ['chr','st','ed','strand']
        def _gbed(ex):
            gr = ex.groupby('_gidx')
            g = gr[gcols].first()
            g['st'] = gr['st'].min()
            g['ed'] = gr['ed'].max()
            return g.reset_index()
        g1 = _gbed(e1)
        g2 = _gbed(e2)
        gcols2 = gcols+['_gidx']
        a2 = UT.save_tsv_nidx_nhead(g1[gcols2],p1.gbed)
        b2 = UT.save_tsv_nidx_nhead(g2[gcols2],p2.gbed)
        c2 = p1.geneovl(p2.code)
        c2 = BT.bedtoolintersect(a2,b2,c2,wao=True)
        gocols = gcols2 + ['b_'+x for x in gcols2] + ['ovl']
        self.gov = gov = UT.read_pandas(c2,names=gocols)
        
    def assign_tcode_ex(self):
        ov = self.ov
        # restrict to proper ones (e.g. strand)
        idxchr = ov['chr']==ov['b_chr'] # all items should satisfy this...
        idxstrand = ov['strand']==ov['b_strand'] # strand should be same
        # above will remove all known SE
        idxstrand = idxstrand|(ov['strand']=='.') # 
        self.eo = eo = ov[idxchr&idxstrand].copy() 
        eo['dlen'] = N.abs(eo['len']-eo['b_len'].astype(int))
        eo = eo.sort_values(['_id','dlen'],ascending=True)
        self.eog = eog = eo.groupby('_id',sort=False).first().reset_index()
        # exons: eknown, etcode
        ex = self.e1
        e2k = dict([(x,'k') for x in eog['_id']])
        ex['eknown_'+self.ref] = [e2k.get(x,'u') for x in ex['_id']]
        idxse = ex['cat']=='s'
        idxek = ex['eknown_'+self.ref]=='k'
        ex.loc[idxek&idxse,'etcode_'+self.ref] = 'k.se'
        ex.loc[idxek&(~idxse),'etcode_'+self.ref] = 'k.me'
        ex.loc[(~idxek)&idxse,'etcode_'+self.ref] = 'u.se'
        ex.loc[(~idxek)&(~idxse),'etcode_'+self.ref] = 'u.me'

    def assign_tcode_as(self):
        # whether exon overlap in anti-sense direction
        # assign code to exon-wise and then gene-wise
        ov = self.ov
        # restrict to anti-sense
        idxchr = ov['chr']==ov['b_chr'] # all items should satisfy this...
        idxstrand = ov['strand']!=ov['b_strand'] # opposite strand 
        # above will include all SE
        idxstrand = idxstrand&(ov['strand']!='.') # 
        self.eoas = eoas = ov[idxchr&idxstrand].copy() 
        # take closest in length
        eoas['dlen'] = N.abs(eoas['len']-eoas['b_len'].astype(int))
        eoas = eoas.sort_values(['_id','dlen'],ascending=True)
        self.eoasg = eoasg = eoas.groupby('_id',sort=False).first().reset_index()
        # exons: eknown, etcode
        ex = self.e1
        e2k = dict([(x,'y') for x in eoasg['_id']])
        ex['ex_as_ovl_'+self.ref] = [e2k.get(x,'n') for x in ex['_id'].values]
        gasovl = ex.groupby('_gidx')['ex_as_ovl_'+self.ref].apply(lambda x: 'y' in set(x))
        g2k = dict([(x,'y') for x in gasovl[gasovl].index.values])
        ex['gene_as_ovl_'+self.ref] = [g2k.get(x,'n') for x in ex['_gidx'].values]

    def assign_intergenic(self):
        # whether gene is intergenic or not
        gov = self.gov  
        ex = self.e1
        gov2 = gov[gov['ovl']>0]
        # ['chr','st','ed','strand','_gidx','b_chr','b_st','b_ed','b_strand','b__gidx','ovl']
        ovlgids = set(gov2['_gidx'].values)
        g2o = dict([(x,'n') for x in ovlgids])
        self.e1['intergenic_'+self.ref] = [g2o.get(x,'y') for x in ex['_gidx'].values]
        
    def assign_tcode_gene(self):
        eog = self.eog
        ex = self.e1
        # genes: gknown, gtcode
        g2k = dict([(x,'k') for x in set(eog['_gidx'].values)])
        ex['gknown_'+self.ref] = [g2k.get(x,'u') for x in ex['_gidx']]
        idxse = ex['cat']=='s'
        idxgk = ex['gknown_'+self.ref]=='k'
        ex.loc[idxgk&idxse,'gtcode_'+self.ref] = 'k.se'
        ex.loc[idxgk&(~idxse),'gtcode_'+self.ref] = 'k.me'
        ex.loc[(~idxgk)&idxse,'gtcode_'+self.ref] = 'u.se'
        ex.loc[(~idxgk)&(~idxse),'gtcode_'+self.ref] = 'u.me'
        # ref_gidxs, ref_gidx0
        eogg = eog.groupby('_gidx')['b__gidx']
        g2cnt = {k:Counter(list(v)) for k,v in eogg}
        g2gi = {k:','.join(map(str,v.keys())) for k,v in g2cnt.items()}
        g2gi0 = {k:v.most_common()[0][0] for k,v in g2cnt.items()}
        #self.eogg = eogg = eog.groupby('_gidx')['b__gidx'].apply(lambda x: Counter(list(x))).reset_index()
        #g2gi = dict([(x, ','.join(map(str, y.keys()))) for x,y in UT.izipcols(eogg, ['_gidx','b__gidx'])])
        #g2gi0 = dict([(x, y.most_common()[0][0]) for x,y in UT.izipcols(eogg, ['_gidx','b__gidx'])])
        #ex['ref_gidxs'] = [g2gi.get(x, N.nan) for x in ex['_gidx']]
        #ex['ref_gidx0'] = [g2gi0.get(x, N.nan) for x in ex['_gidx']]
        # 2016-02-16 
        ex['ref_gidxs_'+self.ref] = [g2gi.get(x, N.nan) for x in ex['_gidx']]
        ex['ref_gidx0_'+self.ref] = [g2gi0.get(x, N.nan) for x in ex['_gidx']]
        # ref_symbols, ref_symbol0
        sfld = self.refnamefld
        self.e2g = e2g = self.e2.groupby('_gidx')[sfld].apply(lambda x: list(set(x))).reset_index()
        self.g2s = g2s = dict(UT.izipcols(e2g, ['_gidx',sfld]))
        g2gs = {k:','.join(set(reduce(iadd,[g2s[int(z)] for z in v.keys()]))) for k,v in g2cnt.items()}
        g2gs0 = {k:','.join(g2s[int(v.most_common()[0][0])]) for k,v in g2cnt.items()}
        #g2gs = dict([(x, ','.join(set(reduce(iadd,[g2s[z] for z in y.keys()])))) for x,y 
        #             in UT.izipcols(eogg, ['_gidx','b__gidx'])])
        #g2gs0 = dict([(x, ','.join(g2s[y.most_common()[0][0]])) for x,y 
        #             in UT.izipcols(eogg, ['_gidx','b__gidx'])])
        #ex['ref_syms'] = [g2gs.get(x, N.nan) for x in ex['_gidx']]
        #ex['ref_sym0'] = [g2gs0.get(x, N.nan) for x in ex['_gidx']]
        # 2016-02-16 save ref_syms, ref_sym0 columns from old ones
        ex['ref_syms_'+self.ref] = [g2gs.get(x, N.nan) for x in ex['_gidx']]
        ex['ref_sym0_'+self.ref] = [g2gs0.get(x, N.nan) for x in ex['_gidx']]
        
    def assign_tcode_sj(self):
        self.s1 = s1 = UT.read_pandas(self.p1.sj)
        self.s2 = s2 = UT.read_pandas(self.p2.sj)
        self.msj = msj =GGB.read_bed(self.mergedsj) # sc1=>u, m
        if 'locus' not in s1.columns:
            s1['locus'] = GGB.calc_locus_strand(s1)
        if 'locus' not in s2.columns:
            s2['locus'] = GGB.calc_locus_strand(s2)
        if 'locus' not in msj.columns:
            msj['locus'] = GGB.calc_locus_strand(msj)
        l2c = dict([(x,'k.me') for x in s2['locus']])
        l2u = dict(UT.izipcols(msj, ['locus','sc1']))
        l2m = dict(UT.izipcols(msj, ['locus','tst']))
        
        s1['etcode_'+self.ref] = [l2c.get(x,'u.me') for x in s1['locus']]
        g2c = dict(UT.izipcols(self.e1,['_gidx','gtcode_'+self.ref]))
        s1['gtcode_'+self.ref] = [g2c.get(x,'u.me') for x in s1['_gidx']]
        
        s1['ucnt'] = [l2u.get(x,0) for x in s1['locus']]
        s1['mcnt'] = [l2m.get(x,0) for x in s1['locus']]
        s1['jcnt'] = [x or y for x,y in s1[['ucnt','mcnt']].values]
        
    def assign_tcode(self):
        self.assign_tcode_ex()
        self.assign_intergenic()
        self.assign_tcode_as()
        self.assign_tcode_gene()
        self.assign_tcode_sj()
        UT.save_tsv_nidx_whead(self.e1, self.p1.ex)
        UT.save_tsv_nidx_whead(self.s1, self.p1.sj)
        
   