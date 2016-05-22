"""

.. module:: findparams
    :synopsis: Using reference annotation and data to find parameters for assembler.

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
# system imports
import subprocess
import multiprocessing
import gzip
import os
import time
from functools import reduce
from operator import iadd, iand
from collections import Counter
from itertools import repeat
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

# 3rd party imports
import pandas as PD
import numpy as N
import matplotlib.pyplot as P

# library imports
from jgem import utils as UT
from jgem import bedtools as BT
from jgem import bigwig as BW
from jgem import gtfgffbed as GGB

from jgem.bigwig import BWObj, BWs 


def find_maxgap(arr, emin, emax, th, win, gapmode):
    if (emin>th):
        return 0
    if (emax<=th):
        return win
    idx0 = N.nonzero(arr>th) # first find furthest point
    maxpos = idx0[0][-1] # after this all zero
    idx = N.nonzero(arr<=th)
    if len(idx[0])==0:
        return 0
    if (gapmode!='i')&(idx[0][0]>maxpos):
        return 0
    cmax = 1
    cst = idx[0][0]
    ced = cst
    for i in idx[0][1:]:
        if (i>maxpos)&(gapmode!='i'):
            break        
        if i==ced+1: # continuous
            ced = i
        else:
            cmax = max(cmax, ced-cst+1)
            cst = ced = i
    cmax = max(cmax, ced-cst+1)
    return cmax

def find_firstgap(arr, emin, emax, th, win):
    if (emin>th):
        return 0, win # no gap
    if (emax<=th):
        return win, 0 # max
    idx = N.nonzero(arr<=th) # gap pos
    cst = idx[0][0] # current gap start
    ced = cst # current gap end
    for i in idx[0][1:]:
        if i==ced+1: # continuous
            ced = i
        else: # end of first gap and pos
            break
    return ced-cst+1, cst
      
def average(arr, n):
    end =  n * int(len(arr)/n)
    return N.mean(arr[:end].reshape(-1, n), 1)


COPYCOLS = ['chr','st','ed','_gidx','locus','gene_id','gene_type']

class ParamFinder(object):
    """
    Args:
        refpre: pathprefix to ref (assume .ex.txt.gz, .sj.txt.gz)

    """
    def __init__(self, refpre, bwpre, genome):
        self.refpre = refpre
        self.bwpre = bwpre
        self.genome = genome
        # .ex.p.bw, .ex.n.bw, .ex.u.bw, .sj.p.bw, .sj.n.bw, .sj.u.bw
        S2S = {'+':'.p','-':'.n','.':'.u'}
        self.bwpaths = bwp = {
            'ex': {s:bwpre+'.ex.{0}.bw'.format(S2S[s]) for s in S2S},
            'sj': {s:bwpre+'.sj.{0}.bw'.format(S2S[s]) for s in S2S},
        }
        self.bws = bws = {'ex':{},'sj':{}}
        bws['ex']['.'] = BWs([bwp['ex']['.']])
        bws['ex']['+'] = BWs([bwp['ex']['+'],bwp['ex']['.']]) if os.path.exists(bwp['ex']['+']) else bws['ex']['.']
        bws['ex']['-'] = BWs([bwp['ex']['-'],bwp['ex']['.']]) if os.path.exists(bwp['ex']['-']) else bws['ex']['.']
        bws['sj']['+'] = BWs([bwp['sj']['+'],bwp['sj']['.']])
        bws['sj']['-'] = BWs([bwp['sj']['-'],bwp['sj']['.']])
        bws['sj']['.'] = BWs([bwp['sj']['.']])

        self.extract_nonovl_exons()
        self.extract_exi53()
        self.extract_53_pair() # intergenic

    def extract_nonovl_exons(self):
        # nonovl exons    
        self.ex = ex = UT.read_pandas(self.refpre+'.ex.txt.gz')
        self.sj = sj = UT.read_pandas(self.refpre+'.sj.txt.gz')
        ex['gene_type'] = ex['extra'].str.split(';').str[2].str.split().str[1].str[1:-1]
        cols0 = ['chr','st','ed','_id']
        a = self.refpre+'.ex.bed.gz'
        a = UT.write_pandas(ex[cols0], a, '')
        b = self.refpre+'.sj.bed.gz'
        b = UT.write_pandas(sj[cols0], b, '')
        c1 = self.refpre+'.ex-ovl-sj.txt.gz'
        c2 = self.refpre+'.ex-ovl-ex.txt.gz'
        c1 = BT.bedtoolintersect(a,b,c1,wao=True)
        c2 = BT.bedtoolintersect(a,a,c2,wo=True)

        cols = cols0+['b_'+x for x in cols0]+['ovl']
        sov = UT.read_pandas(c1, names=cols)
        sov['len'] = sov['ed']-sov['st']
        sov['ovlratio'] = sov['ovl']/sov['len']
        sovg = sov.groupby('_id')['ovlratio'].max()
        snonov = sovg[sovg<1.] # not completely covered by junction

        eov = UT.read_pandas(c2, names=cols)
        eovsize = eov.groupby('_id').size()
        enonov = eovsize[eovsize==1] # only overlaps with self

        LOG.info('#non-ex-ovl-ex={0}, #non-sj-ovl-ex={1}'.format(len(enonov), len(snonov)))
        ids = set(enonov.index).intersection(snonov.index)
        LOG.info('#non-ovl-ex={0}'.format(len(ids)))
        self.nov_ex = novex = ex.set_index('_id').ix[ids].sort_values(['chr','st','ed'])
        self.ne_i = novex[novex['cat']=='i']
        self.ne_5 = novex[novex['cat']=='5']
        self.ne_3 = novex[novex['cat']=='3']
        self.ne_s = novex[novex['cat']=='s']

    def extract_exi53(self):
        # internal exons overlapping with either 5 or 3 prime exons?
        cols0 = ['chr','st','ed','_id']
        cols = cols0+['b_'+x for x in cols0]+['ovl']
        ex = self.ex

        exi = ex[ex['cat']=='i'] # internal exons
        ai = self.refpre + '.exi.bed.gz'
        ai = UT.write_pandas(exi[cols0], ai, '')

        e5 = ex[ex['cat']=='5']
        a5 = self.refpre + '.ex5.bed.gz'
        a5 = UT.write_pandas(e5[cols0], a5, '')

        e3 = ex[ex['cat']=='3']
        a3 = self.refpre + '.ex3.bed.gz'
        a3 = UT.write_pandas(e3[cols0], a3, '')

        a5i = self.refpre + '.ex5-ovl-exi.txt.gz'
        a3i = self.refpre + '.ex3-ovl-exi.txt.gz'

        e5i0 = BT.calc_ovlratio(a5,ai,a5i,4,4)
        e3i0 = BT.calc_ovlratio(a3,ai,a3i,4,4)

        self.e5i = e5i = e5i0[e5i0['ovlratio']==1]
        self.e3i = e3i = e3i0[e3i0['ovlratio']==1]

    def calc_flux(self, exdf):
        ebw = self.bws['ex']
        sbw = self.bws['sj']
        recs = []
        cols = ['_id', 'sdelta','ecovavg','ecovmin','ecovmax','sin','sout']
        for strand in ['+','-','.']:
            exdfsub = exdf[exdf['strand']==strand]
            with ebw[strand]:
                with sbw[strand]:
                    for chrom, st, ed, _id in exdf[['chr','st','ed', '_id']].values:
                        ecov = ebw[strand].get(chrom,st-1,ed+1)
                        scov = sbw[strand].get(chrom,st-1,ed+1)
                        if strand=='+':
                            sd = scov[-1]-scov[0]
                            sin,sout = scov[0],scov[-1]
                        else:
                            sd = scov[0]-scov[-1]
                            sin,sout = scov[-1],scov[0]
                        recs.append([_id, sd, ecov.mean(), ecov.min(), ecov.max(),sin,sout])
        df = PD.DataFrame(recs, columns=cols)
        exdfi = exdf.set_index('_id').ix[df['_id'].values]
        for f in COPYCOLS:
            df[f] =exdfi[f].values
        return df

    def extract_35_pair(ex, tmpprefix):
        ex['_apos'] = ex['a_pos'].str.split(':').str[1].astype(int)
        ex['_dpos'] = ex['d_pos'].str.split(':').str[1].astype(int)
        ex.loc[ex['cat']=='3','spos'] = ex['_apos']
        ex.loc[ex['cat']=='5','spos'] = ex['_dpos']
        cols = ['chr','st','ed','name','strand','_gidx1','_gidx2']
        def _find(ecs, chrom, strand):
            asc = strand=='+'
            e53 = ecs[ecs['cat'].isin(['3','5'])].sort_values('spos', ascending=asc)
            esorted = echrstrand.sort_values('_apos', ascending=asc)
            v1 = e53.iloc[:-1][['spos','cat','_gidx','_id']].values
            v2 = e53.iloc[1:][['spos','cat','_gidx','_id']].values
            pairs = []
            if strand=='+':
                for r1,r2 in zip(v1,v2):
                    if (r1[1]=='3')&(r2[1]=='5'): # 
                        name = 'e{1}-e{3}'.format(r1[2],r1[3],r2[2],r2[3])
                        pairs.append((chrom,r1[0],r2[0],name,strand,r1[2],r2[2]))
            else:
                for r1,r2 in zip(v1,v2):
                    if (r1[1]=='3')&(r2[1]=='5'): # 
                        name = 'g{0}.e{1}-g{2}.e{3}'.format(r1[2],r1[3],r2[2],r2[3])
                        pairs.append((chrom,r2[0],r1[0],name,strand,r1[2],r2[2]))

            df = PD.DataFrame(pairs, columns=cols)
            return df
        rslts = []
        for chrom in ex['chr'].unique():
            for strand in ['+','-']:
                echrstrand = ex[(ex['chr']==chrom)&(ex['strand']==strand)]
                rslts.append(_find(echrstrand, chrom, strand))
        df = PD.concat(rslts, ignore_index=True).sort_values(['chr','st','ed'])
        # intersect with internal exons
        a = tmpprefix+'.exi.bed' # ncol 3
        b = tmpprefix+'.53.bed' #ncol 5
        c = tmpprefix+'.i53ovl.txt'
        exi = ex[ex['cat']=='i'].sort_values(['chr','st','ed'])
        UT.write_pandas(exi[['chr','st','ed']], a, '')
        UT.write_pandas(df, b, '')
        c = BT.bedtoolintersect(b, a, c, wao=True)
        cols1 = cols+['b_chr','b_st','b_ed','ovl']
        cdf = UT.read_pandas(c, names=cols1)
        sdf = cdf[cdf['ovl']==0][cols]
        sdf['locus'] = UT.calc_locus(sdf)
        sdf['len'] = sdf['ed']-sdf['st']
        UT.write_pandas(sdf, tmpprefix+'.e53pair.bed.gz')
        return sdf

    def calc_params_mp(self, beddf,  win=600, siz=10, covfactor=1e-3, direction='>', gapmode='53', np=10):
        chroms = UT.chroms(self.genome)
        args = []
        for c in chroms:
            bedc = beddf[beddf['chr']==c]
            args.append((bedc, self.bws, win, siz, covfactor, direction, gapmode))
        rslts = UT.process_mp(calc_params_chr, args, np=np, doreduce=False)
        df = PD.concat(rslts, ignore_index=True)
        return df
        
    def parseplot(self, locus, figsize=(15,6)):
        bws = self.bws
        chrom,tmp,strand = locus.split(':')
        st,ed = [int(x.replace(',','')) for x in tmp.split('-')]
        print(locus,st,ed,strand)
        a1 = bws['ex'][strand].get_as_array(chrom,st,ed)
        b1 = bws['sj'][strand].get_as_array(chrom,st,ed)
        fig,axr = P.subplots(2,1,figsize=figsize)
        axr[0].plot(a1)
        axr[0].plot(b1)
        axr[1].plot(a1)
        axr[1].plot(a1+b1,'r--')
    
    def plotex(self, exrec, win=50, figsize=(15,6)):
        chrom,st,ed,strand = exrec[['chr','st','ed','strand']].values
        print(chrom,st,ed,strand)
        bws = self.bws
        a1 = bws['ex'][strand].get_as_array(chrom,st-win,ed+win)
        b1 = bws['sj'][strand].get_as_array(chrom,st-win,ed+win)
        fig,ax = P.subplots(1,1,figsize=figsize)
        ax.plot(a1)
        ax.plot(b1)
        # ex: max,min,maxl,maxr,exl10,exr10 sj: sjl10, sjr10
        dic = dict(
            exmax = N.max(a1[win:-win]),
            exmin = N.min(a1[win:-win]),
            maxl = N.max(a1[:win]),
            maxr = N.max(a1[-win:]),
            exl10 = N.mean(a1[win:win+10]),
            exr10 = N.mean(a1[-win-10:-win]),
            sjl10 = N.mean(b1[win-10:win]),
            sjr10 = N.mean(b1[-win:-win+10]),
        )
        return dic
    
    def pltscatter(self, eids, exdf, sidf, fld='d_id', ecov='ecov', jcov='jcov'):
        etgt = exdf.set_index('_id').ix[eids]
        #e1d = etgt[etgt[fld]!=0].groupby(fld)['ecov_d_m35'].mean()
        e1d = etgt[etgt[fld]!=0].set_index(fld)[ecov]
        s1d = g4s1.groupby(fld)[jcov].sum()
        print(len(e1d), len(s1d))
        # all exons
        x = N.log2(e1d.values+1)
        y = N.log2(s1d.ix[e1d.index].values+1)
        P.plot(x, y, '.', ms=1);
        x0 = N.linspace(0,20)
        P.plot(x0,x0,'r')
        P.plot(x0,x0+1.5,'r--')
        P.plot(x0,x0-1.5,'r--')
        P.xlabel('ecov mean')
        P.ylabel('jcov sum')
        return x,y

    def analyze_i(self, t, bnum=1000):
        xfld = 'emax'
        yfld = 'gap'
        ymax = 14
        #t = pi
        #bnum = 1000
        y2fld = 'mp'
        y2max = 1.1

        t['len'] = t['ed']-t['st']
        t['mp'] = 1. - t[['gap','len']].min(axis=1)/t['len']
        t['mima'] = t['emin']/t['emax']
        t = t[t['emax']>0]
        
        fig,axr = P.subplots(1,3,figsize=(12,4))
        ax = axr[0]
        x = N.log2(t[xfld].values+1)
        y = N.log2(t['gap'].values+1)
        ax.plot(x,y,'.',ms=5,alpha=0.1)
        ax.plot([0,ymax],[0,ymax],'g--')
        ax.set_ylim([-1,ymax])
        avx,avy = UT.calc_binned(x,y,num=bnum)
        ax.plot(avx,avy,'r.-')
        ax.set_title('gap')

        ax = axr[1]
        x = N.log2(t[xfld].values+1)
        y = N.log2(t['mima'].values+1)
        ax.plot(x,y,'.',ms=5,alpha=0.1)
        ax.plot([0,y2max],[0,y2max],'g--')
        ax.set_ylim([-0.1,y2max])
        avx,avy = UT.calc_binned(x,y,num=bnum)
        ax.plot(avx,avy,'r.-')
        ax.set_title('mima')

        ax = axr[2]
        x = N.log2(t['sIn'].values+1)
        y = N.log2(t['emin'].values+1)
        ax.plot(x,y,'.',ms=5,alpha=0.1)
        ax.plot([0,ymax],[0,ymax],'g--')
        # ax.set_ylim([-0.1,1.1])
        avx,avy = UT.calc_binned(x,y,num=bnum)
        ax.plot(avx,avy,'r.-')
        ax.set_title('emin')
    
    


def calc_params_chr(exdf, bws, win=300, siz=10, covfactor=1e-3, direction='>', gapmode='i'):
    ebw = bws['ex']
    sbw = bws['sj']
    recs = []
    cols = ['emax','emin',
            'emaxIn','eminIn','gapIn','gposIn',
            'emaxOut','eminOut','gapOut','gposOut',
            'eIn','sIn',
            'eOut','sOut',
            'gap']
    for strand in ['+','-','.']:
        exdfsub = exdf[exdf['strand']==strand]
        with ebw[strand]:
            with sbw[strand]:
                for chrom,st,ed in exdfsub[['chr','st','ed']].values:
                    #win = ed-st # same size as exon
                    left = max(0, st-win)
                    if left==0:
                        print('st-win<0:{0}:{1}-{2}'.format(chrom,st,ed))
                    right = ed+win
                    stpos = st-left
                    edpos = ed-left
                    a1 = ebw[strand].get_as_array(chrom,left,right)
                    b1 = sbw[strand].get_as_array(chrom,left,right)
                    a1[N.isnan(a1)] = 0.
                    b1[N.isnan(b1)] = 0.
                     
                    exl10 = N.mean(a1[stpos:stpos+siz])
                    sjl10 = N.mean(b1[stpos-siz:stpos])
                    exr10 = N.mean(a1[edpos-siz:edpos])
                    sjr10 = N.mean(b1[edpos:edpos+siz])
                    exmax = N.max(a1[stpos:edpos])
                    exmin = N.min(a1[stpos:edpos])
                    #gapth = sjl10*covfactor if strand=='+' else sjr10*covfactor
                    gapth = exmax*covfactor
                    if ((direction=='>')&(strand=='+'))|((direction!='>')&(strand=='-')):
                        gap = find_maxgap(a1[stpos:edpos],exmin, exmax, gapth, win, gapmode)
                    else:
                        gap = find_maxgap(a1[stpos:edpos][::-1],exmin, exmax, gapth, win, gapmode)
                    maxl = N.max(a1[:stpos])
                    maxr = N.max(a1[edpos:])
                    minl = N.min(a1[:stpos])
                    minr = N.min(a1[edpos:])
                    gapl,posl = find_firstgap(a1[:stpos][::-1],minl,maxl,gapth,win)
                    gapr,posr = find_firstgap(a1[edpos:],minr,maxr,gapth,win)
                    if strand=='+':
                        recs.append([exmax,exmin,
                                     maxl,minl,gapl,posl, 
                                     maxr,minr,gapr,posr, 
                                     exl10,sjl10,
                                     exr10,sjr10,
                                     gap])
                    else:
                        recs.append([exmax,exmin,
                                     maxr,minr,gapr,posr, 
                                     maxl,minl,gapl,posl, 
                                     exr10,sjr10,
                                     exl10,sjl10,
                                     gap])

    df = PD.DataFrame(recs, columns=cols)
    exdfi = exdf.set_index('_id').ix[df['_id'].values]
    for f in COPYCOLS:
        df[f] =exdfi[f].values
    return df
