"""

.. module:: merge3
    :synopsis: merge assemblies from different cell types
    jGEM version 3 merger

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
from jgem import assembler3 as A3

import jgem.cy.bw as cybw


class PrepBWSJ(object):
    
    def __init__(self, j2pres, genome, dstpre, libsizes=None, np=10):
        self.j2pres = j2pres
        self.libsizes = libsizes # scale = 1e6/libsize
        self.genome = genome
        self.dstpre = dstpre
        self.np = np
        
    def __call__(self):
        # exdf => ex.p, ex.n, ex.u
        # sjdf => sj.p, sj.n, sj.u
        # paths => sjpath.bed
        # divide into tasks (exdf,sjdf,paths) x chroms
        self.server = server = TQ.Server(name='PrepBWSJ', np=self.np)
        self.chroms = chroms = UT.chroms(self.genome)
        csizes = UT.df2dict(UT.chromdf(self.genome), 'chr', 'size')
        self.exstatus = exstatus = {}
        self.sjstatus = sjstatus = {}
        self.pastatus = pastatus = {}
        self.sdstatus = sdstatus = {}
        exdone=False
        sjdone=False
        padone=False
        sddone=False
        with server:
            for chrom in chroms:
                # exdf tasks
                tname = 'prep_exwig_chr.{0}'.format(chrom)
                args = (self.j2pres, self.libsizes, self.dstpre, chrom, csizes[chrom])
                task = TQ.Task(tname, prep_exwig_chr, args)
                server.add_task(task)
                # exdf tasks
                tname = 'prep_sjwig_chr.{0}'.format(chrom)
                args = (self.j2pres, self.libsizes, self.dstpre, chrom, csizes[chrom])
                task = TQ.Task(tname, prep_sjwig_chr, args)
                server.add_task(task)
                # exdf tasks
                tname = 'prep_sjpath_chr.{0}'.format(chrom)
                args = (self.j2pres, self.libsizes, self.dstpre, chrom)
                task = TQ.Task(tname, prep_sjpath_chr, args)
                server.add_task(task)
                tname = 'prep_sjdf_chr.{0}'.format(chrom)
                args = (self.j2pres, self.libsizes, self.dstpre, chrom)
                task = TQ.Task(tname, prep_sjdf_chr, args)
                server.add_task(task)
            while server.check_error():
                try:
                    name, rslt = server.get_result(timeout=5) # block until result come in
                except TQ.Empty:
                    name, rslt = None, None
                if name is not None:
                    if name.startswith('prep_exwig_chr.'):
                        chrom = name.split('.')[1]
                        exstatus[chrom] = rslt
                        if len(exstatus)==len(chroms): # all finished
                            print('$$$$$$$$ putting in prep_exbw $$$$$$$$$$$')
                            tname='prep_exbw'
                            args = (self.dstpre, chroms, self.genome)
                            task = TQ.Task(tname, prep_exbw, args)
                            server.add_task(task)
                    if name.startswith('prep_sjwig_chr.'):
                        chrom = name.split('.')[1]
                        sjstatus[chrom] = rslt
                        if len(sjstatus)==len(chroms): # all finished
                            print('$$$$$$$$ putting in prep_sjbw $$$$$$$$$$$')
                            tname='prep_sjbw'
                            args = (self.dstpre, chroms, self.genome)
                            task = TQ.Task(tname, prep_sjbw, args)
                            server.add_task(task)
                    if name.startswith('prep_sjpath_chr.'):
                        chrom = name.split('.')[1]
                        pastatus[chrom] = rslt
                        if len(pastatus)==len(chroms): # all finished
                            print('$$$$$$$$ putting in prep_sjpath $$$$$$$$$$$')
                            tname='prep_sjpath'
                            args = (self.dstpre, chroms)
                            task = TQ.Task(tname, prep_sjpath, args)
                            server.add_task(task)
                    if name.startswith('prep_sjdf_chr.'):
                        chrom = name.split('.')[1]
                        sdstatus[chrom] = rslt
                        if len(sdstatus)==len(chroms): # all finished
                            print('$$$$$$$$ putting in prep_sjdf $$$$$$$$$$$')
                            tname='prep_sjdf'
                            args = (self.dstpre, chroms)
                            task = TQ.Task(tname, prep_sjdf, args)
                            server.add_task(task)
                    if name=='prep_exbw':
                        print('$$$$$$$$ prep_exbw done $$$$$$$$$$$')
                        exdone=True
                    if name=='prep_sjbw':
                        print('$$$$$$$$ prep_sjbw done $$$$$$$$$$$')
                        sjdone=True
                    if name=='prep_sjpath':
                        print('$$$$$$$$ prep_sjpath done $$$$$$$$$$$')
                        padone=True
                    if name=='prep_sjdf':
                        print('$$$$$$$$ prep_sjdf done $$$$$$$$$$$')
                        sddone=True
                    if exdone&sjdone&padone&sddone:
                        break
            print('Exit Loop')
        print('Done')
                        


def prep_exwig_chr(j2pres, libsizes, dstpre, chrom, csize):
    ss = ['p','n','u']
    s2s = {'p':['+'],'n':['-'],'u':['.+','.-','.']}
    a = {s:N.zeros(csize) for s in ss}
    wigpaths = {s:dstpre+'.ex.{0}.{1}.wig'.format(s,chrom) for s in ss}
    if all([os.path.exists(dstpre+'.ex.{0}.bw'.format(s)) for s in ss]):
        return wigpaths
    if all([os.path.exists(dstpre+'.ex.{0}.wig'.format(s)) for s in ss]):
        return wigpaths
    if all([os.path.exists(wigpaths[s]) for s in ss]):
        return wigpaths
    if libsizes is None:
        n = 1
        scales = N.ones(len(j2pres))
    else:
        n = len(j2pres)
        scales = [1e6/float(x) for x in libsizes]
    for pre,scale in zip(j2pres, scales):
        exdf = UT.read_pandas(pre+'.exdf.txt.gz',names=A3.EXDFCOLS)
        exdf = exdf[exdf['chr']==chrom]
        for s in ss:
            exsub = exdf[exdf['strand'].isin(s2s[s])]
            for st,ed,ecov in exsub[['st','ed','ecov']].values:
                a[s][st:ed] += ecov*scale
        sedf = UT.read_pandas(pre+'.sedf.txt.gz',names=A3.EXDFCOLS)
        sedf = sedf[sedf['chr']==chrom]
        for s in ss:
            sesub = sedf[sedf['strand'].isin(s2s[s])]
            for st,ed,ecov in sesub[['st','ed','ecov']].values:
                a[s][st:ed] += ecov*scale
    for s in ['p','n','u']:
        if libsizes is not None:
            a[s] /= float(n) # average
        cybw.array2wiggle_chr64(a[s], chrom,  wigpaths[s], 'w')
    return wigpaths  

def prep_sjwig_chr(j2pres, libsizes, dstpre, chrom, csize):
    ss = ['p','n','u']
    s2s = {'p':['+'],'n':['-'],'u':['.+','.-']}
    a = {s:N.zeros(csize) for s in ss}
    wigpaths = {s:dstpre+'.sj.{0}.{1}.wig'.format(s,chrom) for s in ss}
    if all([os.path.exists(dstpre+'.sj.{0}.bw'.format(s)) for s in ss]):
        return wigpaths
    if all([os.path.exists(dstpre+'.sj.{0}.wig'.format(s)) for s in ss]):
        return wigpaths
    if all([os.path.exists(wigpaths[s]) for s in ss]):
        return wigpaths
    if libsizes is None:
        n = 1
        scales = N.ones(len(j2pres))
    else:
        n = len(j2pres)
        scales = [1e6/float(x) for x in libsizes]
    for pre,scale in zip(j2pres, scales):
        sjdf = UT.read_pandas(pre+'.sjdf.txt.gz',names=A3.SJDFCOLS)
        sjdf = sjdf[sjdf['chr']==chrom]
        for s in ss:
            sjsub = sjdf[sjdf['strand'].isin(s2s[s])]
            for st,ed,tcnt in sjsub[['st','ed','tcnt']].values:
                a[s][st:ed] += tcnt*scale
    for s in ['p','n','u']:
        if libsizes is not None:
            a[s] /= float(n) # average
        cybw.array2wiggle_chr64(a[s], chrom,  wigpaths[s], 'w')
    return wigpaths    

def prep_sjpath_chr(j2pres, libsizes, dstpre, chrom):
    pc2st = {}
    pc2ed = {}
    pc2tst = {}
    pc2ted = {}
    pc2strand = {}
    pc2tcov = {}
    # pc2tcov0 = {}
    # chr,st,ed,name,sc1(tcov),strand,tst,ted,sc2(),#exons,estarts,esizes
    # cols = ['st','ed','name','strand','tst','ted','tcov0','tcov']
    path = dstpre+'.sjpath.{0}.bed.gz'.format(chrom)
    path0 = dstpre+'.sjpath.bed.gz'
    if os.path.exists(path0):
        return path
    if os.path.exists(path):
        return path
    
    cols = ['st','ed','name','strand','tst','ted','tcov']

    if libsizes is None:
        n = 1
        scales = N.ones(len(j2pres))
    else:
        n = len(j2pres)
        scales = [1e6/float(x) for x in libsizes]
    for pre,scale in zip(j2pres, scales):
        paths = UT.read_pandas(pre+'.paths.txt.gz', names=A3.PATHCOLS)
        paths = paths[paths['chr']==chrom]
        for st,ed,name,s,tst,ted,tcov in paths[cols].values:
            pc = ','.join(name.split(',')[1:-1]) # trim 53exons => intron chain
            if pc=='':
                continue # ignore no junction path
            pc2st[pc] = min(st, pc2st.get(pc,st))
            pc2ed[pc] = max(ed, pc2ed.get(pc,ed))
            pc2tst[pc] = tst
            pc2ted[pc] = ted
            pc2strand[pc] = s
            pc2tcov[pc] = pc2tcov.get(pc,0)+scale*tcov
            #pc2tcov0[pc] = pc2tcov0.get(pc,0)+scale*tcov0
    df = PD.DataFrame({'st':pc2st,'ed':pc2ed,'tst':pc2tst,'ted':pc2ted,
                       'strand':pc2strand,'tcov':pc2tcov})
    df['chr'] = chrom
    df.index.name = 'name'
    df.reset_index(inplace=True)
    # create bed12: parse name => #exons, esizes, estarts
    df['pc'] = df['name'].copy()
    idxp = df['strand'].isin(['+','.+'])
    if libsizes is not None:
        df['tcov'] = df['tcov']/float(n)
    df.loc[idxp,'name'] = ['{0},{1},{2}'.format(s,p,e) for s,p,e in df[idxp][['st','pc','ed']].values]
    df.loc[~idxp,'name'] = ['{2},{1},{0}'.format(s,p,e) for s,p,e in df[~idxp][['st','pc','ed']].values]
    df = df.groupby('pc').first() # get rid of unstranded duplicates
    cmax = 9+N.log2(N.mean(scales))
    bed = A3.path2bed12(df, cmax)
    # reset sc1 to tcov (from log2(tcov+2)*100)
    bed['sc1'] = bed['tcov']
    GGB.write_bed(bed, path, ncols=12)
    return path
    
def prep_sjdf_chr(j2pres, libsizes, dstpre, chrom):
    pc2st = {}
    pc2ed = {}
    pc2strand = {}
    pc2tcnt = {}
    pc2ucnt = {}
    # chr,st,ed,name,sc1(tcov),strand,tst,ted,sc2(),#exons,estarts,esizes
    # cols = ['st','ed','name','strand','tst','ted','tcov0','tcov']
    path = dstpre+'.sjdf.{0}.txt.gz'.format(chrom)
    path0 = dstpre+'.sjdf.txt.gz'
    if os.path.exists(path0):
        return path
    if os.path.exists(path):
        return path
    
    cols = ['st','ed','name','strand','st','ed','tcnt','ucnt']
    # cols = A3.SJDFCOLS

    if libsizes is None:
        n = 1
        scales = N.ones(len(j2pres))
    else:
        n = len(j2pres)
        scales = [1e6/float(x) for x in libsizes]
    for pre,scale in zip(j2pres, scales):
        paths = UT.read_pandas(pre+'.sjdf.txt.gz', names=A3.SJDFCOLS)
        paths = paths[paths['chr']==chrom]
        for st,ed,pc,s,st,ed,tcnt,ucnt in paths[cols].values:
            pc2st[pc] = st
            pc2ed[pc] = ed
            pc2strand[pc] = s
            pc2tcnt[pc] = pc2tcnt.get(pc,0)+scale*tcnt
            pc2ucnt[pc] = pc2ucnt.get(pc,0)+scale*ucnt
    df = PD.DataFrame({'st':pc2st,'ed':pc2ed,'st':pc2st,'ed':pc2ed,
                       'strand':pc2strand,'tcnt':pc2tcnt,'ucnt':pc2ucnt})
    df['chr'] = chrom
    df['kind'] = 'j'
    if libsizes is not None:
        df['tcnt'] = df['tcnt']/float(n)
        df['ucnt'] = df['ucnt']/float(n)
    df.index.name = 'name'
    df.reset_index(inplace=True)
    UT.write_pandas(df[A3.SJDFCOLS], path, '')
    return path

def prep_exbw(dstpre, chroms, genome):
    return _prep_bw(dstpre, chroms, genome, 'ex')

def prep_sjbw(dstpre, chroms, genome):
    return _prep_bw(dstpre, chroms, genome, 'sj')

def _prep_bw(dstpre, chroms, genome, w):
    # concatenate
    ss = ['p','n','u']
    files = []
    bwpaths = {s: dstpre+'.{1}.{0}.bw'.format(s,w) for s in ss}
    if all([os.path.exists(bwpaths[s]) for s in ss]):
        return bwpaths
    for s in ss:
        dstwig = dstpre+'.{1}.{0}.wig'.format(s,w)
        with open(dstwig, 'wb') as dst:
            for c in chroms:
                srcpath = dstpre+'.{2}.{0}.{1}.wig'.format(s,c,w)
                with open(srcpath,'rb') as src:
                    shutil.copyfileobj(src,dst)
                files.append(srcpath)
        files.append(dstwig)         
        print('converting wig to bigwig {0}'.format(dstwig))
        BT.wig2bw(dstwig, UT.chromsizes(genome), bwpaths[s])
    # clean up        
    for f in files:
        os.unlink(f)
    return bwpaths

def prep_sjpath(dstpre, chroms):
    dstpath = dstpre+'.sjpath.bed.gz'
    if os.path.exists(dstpath):
        return dstpath
    files = []
    with open(dstpath, 'wb') as dst:
        for c in chroms:
            srcpath = dstpre+'.sjpath.{0}.bed.gz'.format(c)
            with open(srcpath,'rb') as src:
                shutil.copyfileobj(src,dst)
            files.append(srcpath)
    # for f in files: # keep separate chr files 
    #     os.unlink(f)
    return dstpath

def prep_sjdf(dstpre, chroms):
    dstpath = dstpre+'.sjdf.txt.gz'
    if os.path.exists(dstpath):
        return dstpath
    files = []
    with open(dstpath, 'wb') as dst:
        for c in chroms:
            srcpath = dstpre+'.sjdf.{0}.txt.gz'.format(c)
            with open(srcpath,'rb') as src:
                shutil.copyfileobj(src,dst)
            files.append(srcpath)
    # for f in files: # keep separate chr files 
    #     os.unlink(f)
    return dstpath

SJFILTERPARAMS = dict(
    th_detected=1,
    th_maxcnt=1,
    th_maxoverhang=15,
    th_minedgeexon=15,
    th_sjratio=1e-3,
    filter_unstranded=False,# there are substantial number of high cov unstranded
)
class SJFilter(object):

    def __init__(self, bwsjpre, statspath, genome, np=10, **kw):
        self.bwsjpre = bwsjpre
        self.statspath = statspath
        self.genome = genome
        self.np = np
        self.params = SJFILTERPARAMS.copy()
        self.params.update(kw)

    def __call__(self):
        chroms = UT.chroms(self.genome)
        csizedic = UT.df2dict(UT.chromdf(self.genome), 'chr', 'size')
        args = []
        for c in chroms:
            csize = csizedic[c]
            args.append((self.bwsjpre, self.statspath, c, csize, self.params))
        
        rslts = UT.process_mp(filter_sjpath, args, np=self.np, doreduce=False)
        dstpath = self.bwsjpre+'.sjpath.filtered.bed.gz'
        with open(dstpath,'wb') as dst:
            for c in chroms:
                srcpath =  self.bwsjpre+'.sjpath.{0}.filtered.bed.gz'.format(c)
                with open(srcpath, 'rb') as src:
                    shutil.copyfileobj(src, dst)

        rslts = UT.process_mp(filter_sjdf, args, np=self.np, doreduce=False)
        dstpath = self.bwsjpre+'.sjdf.filtered.txt.gz'
        with open(dstpath,'wb') as dst:
            for c in chroms:
                srcpath =  self.bwsjpre+'.sjdf.{0}.filtered.txt.gz'.format(c)
                with open(srcpath, 'rb') as src:
                    shutil.copyfileobj(src, dst)



def locus2pc(l):
    chrom,sted,strand = l.split(':')
    st,ed = sted.split('-')
    st = str(int(st)-1)
    if strand in ['+','.']:
        return  '|'.join([st,ed])
    return '|'.join([ed,st])

def filter_sjpath(bwsjpre, statspath, chrom, csize, params):
    # read in junction stats
    stats = UT.read_pandas(statspath)
    if 'chr' not in stats:
        stats['chr'] = [x.split(':')[0] for x in stats['locus']]
    if '#detected' in stats:
        stats.rename(columns={'#detected':'detected'}, inplace=True)
    stats = stats[stats['chr']==chrom].copy()
    if 'pc' not in stats:
        stats['pc'] = [locus2pc(x) for x in stats['locus']]
    flds = ['detected','maxcnt','maxoverhang']
    dics = {f: UT.df2dict(stats, 'pc', f) for f in flds}
    # read sjpath
    fpath_chr =  bwsjpre+'.sjpath.{0}.bed.gz'.format(chrom)
    dstpath = bwsjpre+'.sjpath.{0}.filtered.bed.gz'.format(chrom)
    if os.path.exists(fpath_chr):
        sj = GGB.read_bed(fpath_chr)
    else:
        fpath = bwsjpre+'.sjpath.bed.gz'
        sj = GGB.read_bed(fpath)
        sj = sj[sj['chr']==chrom].copy()
    name0 = sj.iloc[0]['name']
    if len(name0.split('|'))<len(name0.split(',')): # exons attached?
        sj['name'] = [','.join(x.split(',')[1:-1]) for x in sj['name']]
    # filter unstranded
    if params['filter_unstranded']:
        sj = sj[sj['strand'].isin(['+','-'])].copy()
    # filter with stats
    for f in flds:
        sj[f] = [N.min([dics[f].get(x,0) for x in y.split(',')]) for y in sj['name']]
        sj = sj[sj[f]>params['th_'+f]].copy() # filter 
    # edge exon size    
    sj['eflen'] = [int(x.split(',')[0]) for x in sj['esizes']]
    sj['ellen'] = [int(x.split(',')[-2]) for x in sj['esizes']]    
    eth = params['th_minedgeexon']
    sj = sj[(sj['eflen']>eth)&(sj['ellen']>eth)].copy()
    # calculate sjratio, sjratio
    if params['filter_unstranded']:
        sjexbw = A3.SjExBigWigs(bwsjpre, mixunstranded=False)
    else:
        sjexbw = A3.SjExBigWigs(bwsjpre, mixunstranded=True)
    with sjexbw:
        sa = sjexbw.bws['sj']['a'].get(chrom,0,csize)
        ea = sjexbw.bws['ex']['a'].get(chrom,0,csize)
    a = sa+ea
    # sj['sjratio'] = [x/N.mean(a[int(s):int(e)]) for x,s,e in sj[['sc1','tst','ted']].values]
    sj['sjratio'] = [x/N.max(a[int(s):int(e)]) for x,s,e in sj[['sc1','tst','ted']].values]
    sj = sj[sj['sjratio']>params['th_sjratio']]
    GGB.write_bed(sj, dstpath, ncols=12)

def filter_sjdf(bwsjpre, statspath, chrom, csize, params):
    # read in junction stats
    stats = UT.read_pandas(statspath)
    if 'chr' not in stats:
        stats['chr'] = [x.split(':')[0] for x in stats['locus']]
    if '#detected' in stats:
        stats.rename(columns={'#detected':'detected'}, inplace=True)
    stats = stats[stats['chr']==chrom].copy()
    if 'pc' not in stats:
        stats['pc'] = [locus2pc(x) for x in stats['locus']]
    flds = ['detected','maxcnt','maxoverhang']
    dics = {f: UT.df2dict(stats, 'pc', f) for f in flds}
    # read sjdf
    fpath_chr =  bwsjpre+'.sjdf.{0}.txt.gz'.format(chrom)
    dstpath = bwsjpre+'.sjdf.{0}.filtered.txt.gz'.format(chrom)
    if os.path.exists(fpath_chr):
        sj = UT.read_pandas(fpath_chr, names=A3.SJDFCOLS)
    else:
        fpath = bwsjpre+'.sjdf.txt.gz'
        sj = UT.read_pandas(fpath, names=A3.SJDFCOLS)
        sj = sj[sj['chr']==chrom].copy()
    # filter unstranded
    if params['filter_unstranded']:
        sj = sj[sj['strand'].isin(['+','-'])].copy()
    # filter with stats
    for f in flds:
        # sj[f] = [N.min([dics[f].get(x,0) for x in y.split(',')]) for y in sj['name']]
        sj[f] = [dics[f].get(y,0) for y in sj['name']]
        sj = sj[sj[f]>params['th_'+f]].copy() # filter 
    # edge exon size    
    # sj['eflen'] = [int(x.split(',')[0]) for x in sj['esizes']]
    # sj['ellen'] = [int(x.split(',')[-2]) for x in sj['esizes']]    
    # eth = params['th_minedgeexon']
    # sj = sj[(sj['eflen']>eth)&(sj['ellen']>eth)].copy()
    # calculate sjratio, sjratio
    if params['filter_unstranded']:
        sjexbw = A3.SjExBigWigs(bwsjpre, mixunstranded=False)
    else:
        sjexbw = A3.SjExBigWigs(bwsjpre, mixunstranded=True)
    with sjexbw:
        sa = sjexbw.bws['sj']['a'].get(chrom,0,csize)
        ea = sjexbw.bws['ex']['a'].get(chrom,0,csize)
    a = sa+ea
    # sj['sjratio'] = [x/N.mean(a[int(s):int(e)]) for x,s,e in sj[['tcnt','st','ed']].values]
    sj['sjratio'] = [x/N.max(a[int(s):int(e)]) for x,s,e in sj[['tcnt','st','ed']].values]
    sj = sj[sj['sjratio']>params['th_sjratio']]
    UT.write_pandas(sj[A3.SJDFCOLS], dstpath, '')


class LocalEstimator(A3.LocalAssembler):

    def __init__(self, modelpre, bwpre, chrom, st, ed, dstpre, tcovth):
        self.modelpre = modelpre
        self.tcovth = tcovth
        A3.LocalAssembler.__init__(self, bwpre, chrom, st, ed, dstpre)
        bed12 = GGB.read_bed(modelpre+'.paths.withse.bed.gz')
        idx = (bed12['chr']==chrom)&(bed12['tst']>=st)&(bed12['ted']<=ed)
        self.paths = bed12[idx].copy()
        sj = GGB.read_bed(bwpre+'.sjpath.bed.gz')
        idx0 = (sj['chr']==chrom)&(sj['tst']>=st)&(sj['ted']<=ed)        
        self.sjpaths0 = sj[idx0].copy()        
        # load exdf, sjdf
        sjdf = UT.read_pandas(modelpre+'.sjdf.txt.gz', names=A3.SJDFCOLS)
        exdf = UT.read_pandas(modelpre+'.exdf.txt.gz', names=A3.EXDFCOLS)
        idx = (sjdf['chr']==chrom)&(sjdf['st']>=st)&(sjdf['ed']<=ed)
        self.sjdf = sjdf[idx].copy()
        idx = (exdf['chr']==chrom)&(exdf['st']>=st)&(exdf['ed']<=ed)
        self.exdf = exdf[idx].copy()
        A3.set_ad_pos(self.sjdf, 'sj')
        A3.set_ad_pos(self.exdf, 'ex')
        self.sjexbw = sjexbw = A3.SjExBigWigs(bwpre, None, mixunstranded=True)
        self.arrs = arrs = {}
        with sjexbw: # get bw arrays
            for k in ['ex','sj']:
                arrs[k] = {}
                for s in ['+','-']:
                    arrs[k][s] = sjexbw.bws[k][s].get(chrom, st, ed)

    def process(self):
        self.calculate_ecovs()
        self.calculate_scovs()
        self.estimate_abundance()
        self.write()

    def calculate_scovs(self):
        sj = self.sjdf
        sj0 = self.sjpaths0
        sj0mat = sj0[['sc1','sc2','name']].values
        tmp = [[(sc1,sc2) for sc1,sc2,p in sj0mat if y in p] for y in sj['name']]
        sj['ucnt'] = [N.sum([x[0] for x in y]) for y in tmp]
        sj['tcnt'] = [N.sum([x[1] for x in y]) for y in tmp]
        # idx = sj['tcnt']==0
        # tmp0 = ['{1}|{0}'.format(*y.split('|')) for y in sj[idx]['name']]
        # tmp1 = [N.sum([x for x,p in sj0mat if y in p]) for y in tmp0]
        # sj.loc[idx, 'tcnt'] = tmp1
        idxz = sj['tcnt']==0
        if N.sum(idxz)>0:
            sj.loc[idxz,'tcnt'] = 1e-6
        self.sjdfi = sj.set_index('name')

    def calculate_ecovs(self):
        ex = self.exdf
        o = self.st
        if len(ex)==0:
            return
        ex['ecov'] = N.nan
        for strand in ['+','-']:
            spans = self._get_spans(strand)
            for st,ed in spans:
                idx = (ex['st']>=st)&(ex['ed']<=ed)&(ex['strand'].isin(STRS[strand]))
                es = ex[idx].copy().sort_values(['st','ed'])
                es['tmpeid'] = N.arange(len(es))
                ne = len(es)
                exa = self.arrs['ex'][strand]
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
                    try:
                        ecov,err = nnls(mat, ci['cov'].values)
                        ex.loc[idx,'ecov'] = ecov
                    except:
                        LOG.warning('!!!!!! Exception in NNLS (calculate_ecov) @{0}:{1}-{2}, setting to mean !!!!!!!!!'.format(self.chrom, st, ed))
                        ex.loc[idx,'ecov'] = cov(st,ed)
                elif ne==1:
                    s,e = es.iloc[0][['st','ed']]
                    ex.loc[idx,'ecov'] = cov(s,e)
        idxz = ex['ecov']==0
        ex.loc[idxz, 'ecov'] = 1e-6
        self.exdfi = ex.set_index('name')

    def calculate_branchp(self, jids, eids):
        sj0 = self.sjdfi
        sj = sj0.ix[jids].reset_index()
        ex0 = self.exdfi
        ex = ex0.ix[eids].reset_index()

        dsump = sj.groupby('dpos')['tcnt'].sum().astype(float)
        jdp = sj['tcnt'].values/(dsump.ix[sj['dpos'].values].values)
        j2p = dict(zip(sj['name'].values, jdp))
        # exon groupby acceptor
        asump = ex.groupby('apos')['ecov'].sum().astype(float)
        eap = ex['ecov'].values/(asump.ix[ex['apos'].values].values)
        e2ap = dict(zip(ex['name'].values, eap))
        dsump = ex.groupby('dpos')['ecov'].sum().astype(float)
        edp = ex['ecov'].values/(dsump.ix[ex['dpos'].values].values)
        e2dp = dict(zip(ex['name'].values, edp))

        return j2p, e2ap, e2dp

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
            try:
                ecov,err = nnls(mat, ci['cov'].values)
                pg['tcov0'] = ecov
            except:
                # too much iteration?
                LOG.warning('!!!!!! Exception in NNLS (tcov_by_nnls) @{0}:{1}-{2}, setting to zero !!!!!!!!!'.format(self.chrom, s, e))
                pg['tcov0a'] = 0
            pg.rename(columns={'st':'tst','ed':'ted'}, inplace=True)
        else:
            s,e = pg.iloc[0][['tst','ted']]
            pg['tcov0a'] = cov0(s,e)
        # cov1, cov2
        if ne>1:
            sts = sorted(set(pg['tst'].values))
            eds = sorted(set(pg['ted'].values))
            nst,ned = len(sts),len(eds)
            mat = N.array([(pg['tst']==x).values for x in sts]+[(pg['ted']==x).values for x in eds], dtype=float)
            c = N.array([cov1s(x) for x in sts]+[cov1e(x) for x in eds])
            # enforce flux conservation: scale up 5'
            stsum = N.sum(c[:nst])
            edsum = N.sum(c[nst:])
            if strand in ['+','.+']:
                c[:nst] = (edsum/(stsum+1e-6))*c[:nst]
            else:
                c[nst:] = (stsum/(edsum+1e-6))*c[nst:]
            ecov,err = nnls(mat, c)
            pg['tcov0b'] = ecov

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
            pg['tcov0c'] = ecov
        else:
            s,e = pg.iloc[0][['tst','ted']]
            pg['tcov0b'] = (cov1s(s)+cov1e(e))/2.
            pg['tcov0c'] = (cov2s(s)-cov2e(e))/2.

        pg['tcov0'] = pg[['tcov0a','tcov0b','tcov0c']].mean(axis=1)
        pg.loc[pg['tcov0']<0,'tcov0'] = 0 # shouldn't really happen
        keys = [tuple(x) for x in p[idx][['tst','ted']].values]
        for f in ['tcov0','tcov0a','tcov0b','tcov0c']:
            p.loc[idx, f] = pg.set_index(['tst','ted']).ix[keys][f].values
        return pg[['chr','tst','ted','strand','tcov0']]
                
    def tcov_by_branchp(self, tst, ted, strand, tcov0):
        p = self.paths
        idx = (p['strand'].isin(STRS[strand]))&(p['tst']==tst)&(p['ted']==ted)
        if N.sum(idx)==0:
            return
        if N.sum(idx)>1:
            # calculate branchp within this group
            jids = set()
            eids = set()
            for n in p[idx]['name']:
                jids.update(n.split(',')[1:-1])
                eids.update(n.split('|'))
            j2p, e2ap, e2dp = self.calculate_branchp(jids, eids)
            def _prob(y):
                epath0 = y.split('|')
                e5 = epath0[0] # use donor p
                epath = epath0[1:] # use acceptor p
                jpath = y.split(',')[1:-1]
                return e2dp[e5]*N.prod([e2ap[x] for x in epath])*N.prod([j2p[x] for x in jpath])
            p.loc[idx,'tcov'] = [tcov0*_prob(y) for y in p[idx]['name']]
        else:
            p.loc[idx,'tcov'] = tcov0

    def estimate_abundance(self):
        # 1) 5-3 group by NNLS
        # 2) within 5-3 group by tree branch prob
        paths = self.paths
        for s in ['+','-']:
            ps = paths[paths['strand'].isin(A3.STRS[s])]
            if len(ps)==0:
                continue
            for chrom,st,ed in UT.union_contiguous(ps[['chr','st','ed']],returndf=False):
                pg = self.tcov_by_nnls(st,ed,s)
                if pg is not None:
                    for chrom,tst,ted,strand,tcov0 in pg.values:
                        self.tcov_by_branchp(tst,ted,strand,tcov0)

    def write(self):
        pre = self.dstpre+'.{0}_{1}_{2}'.format(self.chrom,self.st,self.ed)
        # 1) exon, junctions, allpaths => csv (no header <= to concatenate bundles)
        ecols = A3.EXDFCOLS #['chr','st','ed','strand','name','kind','ecov']
        UT.write_pandas(self.exdf[ecols], pre+'.covs.exdf.txt.gz', '')
        scols = A3.SJDFCOLS #['chr','st','ed','strand','name','kind','tcnt'  ]#,'donor','acceptor','dp','ap']
        UT.write_pandas(self.sjdf[scols], pre+'.covs.sjdf.txt.gz', '')
        pcols = A3.PATHCOLS #['chr','st','ed','name','strand','tst','ted','tcov0','tcov1','tcov']
        UT.write_pandas(self.paths[pcols], pre+'.covs.paths.txt.gz', '')
        # write colored bed12 for tcov > th
        tgt = self.paths[self.paths['tcov']>=self.tcovth].copy()
        self.bed12 = A3.path2bed12(tgt, cmax=9, covfld='tcov')
        GGB.write_bed(self.bed12, pre+'.covs.paths.bed.gz',ncols=12)

def bundle_estimator(modelpre, bwpre, chrom, st, ed, dstpre, tcovth):
    bname = A3.bundle2bname((chrom,st,ed))
    bsuf = '.{0}_{1}_{2}'.format(chrom,st,ed)
    csuf = '.{0}'.format(chrom)
    sufs = ['.covs.exdf.txt.gz',
            '.covs.sjdf.txt.gz',
            '.covs.paths.txt.gz',
            '.covs.paths.bed.gz',
            ]
    done = []
    for x in sufs:
        done.append(os.path.exists(dstpre+bsuf+x) | \
                    os.path.exists(dstpre+csuf+x) | \
                    os.path.exists(dstpre+x) )
    if all(done):
        LOG.info('bunle {0} already done, skipping'.format(bname))
        return bname
    LOG.info('processing bunle {0}'.format(bname))
    la = LocalEstimator(modelpre, bwpre, chrom, st, ed, dstpre, tcovth)
    return la.process()    

def concatenate_bundles(bundles, dstpre):
    # concat results
    sufs = ['covs.exdf.txt.gz', 
           'covs.sjdf.txt.gz',
           'covs.paths.txt.gz',
           'covs.paths.bed.gz',
           ]
    files = []
    for suf in sufs:
        dstpath = '{0}.{1}'.format(dstpre, suf)
        if not os.path.exists(dstpath):
            with open(dstpath, 'wb') as dst:
                for chrom, st, ed in bundles:
                    bname = A3.bundle2bname((chrom,st,ed))
                    srcpath = '{0}.{1}_{2}_{3}.{4}'.format(dstpre, chrom, st, ed, suf)
                    files.append(srcpath)
                    with open(srcpath, 'rb') as src:
                        shutil.copyfileobj(src, dst)
        else:
            files+=['{0}.{1}_{2}_{3}.{4}'.format(dstpre, chrom, st, ed, suf) for chrom,st,ed in bundles]
    # cleanup
    for f in files:
        if os.path.exists(f):
            os.unlink(f)


def estimatecovs(modelpre, bwpre, dstpre, genome, tcovth=1, np=6):
    bed = GGB.read_bed(modelpre+'.paths.withse.bed.gz')
    chroms = bed['chr'].unique()
    csizedic = UT.df2dict(UT.chromdf(genome), 'chr', 'size')
    bundles = []
    args = []
    for chrom in chroms:
        sub = bed[(bed['chr']==chrom)]
        uc = UT.union_contiguous(sub[['chr','st','ed']], returndf=True)
        # total about 30K=> make batch of ~1000
        n = len(uc)
        nb = int(N.ceil(n/1000.))
        for i in range(nb):
            sti = 1000*i
            edi = min(1000*(i+1), len(uc)-1)
            st = max(uc.iloc[sti]['st'] - 100, 0)
            ed = min(uc.iloc[edi]['ed'] + 100, csizedic[chrom])
            args.append([modelpre, bwpre, chrom, st, ed, dstpre, tcovth])
            bundles.append((chrom,st,ed))

    rslts = UT.process_mp(bundle_estimator, args, np=np, doreduce=False)

    concatenate_bundles(bundles, dstpre)



