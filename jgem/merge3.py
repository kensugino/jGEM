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

############# Merge Prep ######################################################

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


############# SJ Filter #######################################################

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
        dstpath = self.bwsjpre+'.filtered.sjpath.bed.gz'
        with open(dstpath,'wb') as dst:
            for c in chroms:
                srcpath =  self.bwsjpre+'.filtered.sjpath.{0}.bed.gz'.format(c)
                with open(srcpath, 'rb') as src:
                    shutil.copyfileobj(src, dst)

        rslts = UT.process_mp(filter_sjdf, args, np=self.np, doreduce=False)
        dstpath = self.bwsjpre+'.filtered.sjdf.txt.gz'
        with open(dstpath,'wb') as dst:
            for c in chroms:
                srcpath =  self.bwsjpre+'.filtered.sjdf.{0}.txt.gz'.format(c)
                with open(srcpath, 'rb') as src:
                    shutil.copyfileobj(src, dst)

        # make sj.bw
        sjfiltered2bw(self.bwsjpre, self.genome, self.np)
        for s in ['p','n','u']:
            src = self.bwsjpre + '.ex.{0}.bw'.format(s)
            dst = self.bwsjpre + '.filtered.ex.{0}.bw'.format(s)
            cmd = ['ln','-s', src, dst]
            subprocess.call(cmd)



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
    dstpath = bwsjpre+'.filtered.sjpath.{0}.bed.gz'.format(chrom)
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
    dstpath = bwsjpre+'.filtered.sjdf.{0}.txt.gz'.format(chrom)
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



def sjfiltered2wig(bwpre, chrom, chromsize):
    a = {'+':N.zeros(chromsize, dtype=N.float64),
         '-':N.zeros(chromsize, dtype=N.float64),
         '.':N.zeros(chromsize, dtype=N.float64)}
    path = bwpre+'.filtered.sjdf.{0}.txt.gz'.format(chrom)
    sjchr = UT.read_pandas(path, names=A3.SJDFCOLS)
    for st,ed,v,strand in sjchr[['st','ed','tcnt','strand']].values:
        a[strand[0]][st:ed] += v
    for strand in a:
        wig = bwpre+'.filtered.sjdf.{0}.{1}.wig'.format(chrom, strand)
        cybw.array2wiggle_chr64(a[strand], chrom, wig)
    return path

def sjfiltered2bw(bwpre, genome, np=12):
    chroms = UT.chroms(genome)
    chromdf = UT.chromdf(genome).sort_values('size',ascending=False)
    chroms = [x for x in chromdf['chr'] if x in chroms]
    chromdic = UT.df2dict(chromdf, 'chr', 'size')
    args = [(bwpre, c, chromdic[c]) for c in chroms]
    rslts = UT.process_mp(sjfiltered2wig, args, np=np, doreduce=False)
    S2N = {'+':'p','-':'n','.':'u'}
    rmfiles = []
    for strand in ['+','-','.']:
        s = S2N[strand]
        wigpath = bwpre+'.filtered.sj.{0}.wig'.format(s)
        with open(wigpath, 'w') as dst:
            for chrom in chroms:
                f = bwpre+'.filtered.sjdf.{0}.{1}.wig'.format(chrom, strand)
                with open(f,'r') as src:
                    shutil.copyfileobj(src, dst)
                rmfiles.append(f)
        bwpath = bwpre+'.filtered.sj.{0}.bw'.format(s)
        BT.wig2bw(wigpath, UT.chromsizes(genome), bwpath)
        rmfiles.append(wigpath)
    for f in rmfiles:
        os.unlink(f)
    

############# Cov Estimator ######################################################

class LocalEstimator(A3.LocalAssembler):

    def __init__(self, modelpre, bwpre, chrom, st, ed, dstpre, tcovth, usegeom=False):
        self.modelpre = modelpre
        self.tcovth = tcovth
        self.usegeom = usegeom
        A3.LocalAssembler.__init__(self, bwpre, chrom, st, ed, dstpre)
        bed12 = GGB.read_bed(modelpre+'.paths.withse.bed.gz')
        assert(all(bed12['tst']<bed12['ted']))
        idx = (bed12['chr']==chrom)&(bed12['tst']>=st)&(bed12['ted']<=ed)
        self.paths = bed12[idx].copy()
        eids = set()
        sids = set()
        for n in self.paths['name']:
            eids.update(n.split('|'))
            sids.update(n.split(',')[1:-1])
        tgt1 = bwpre+'.filtered.{0}.bed.gz'.format(chrom)
        tgt2 = bwpre+'.{0}.bed.gz'.format(chrom)
        tgt3 = bwpre+'.sjpath.bed.gz'
        if os.path.exists(tgt1):
            sj = GGB.read_bed(tgt1)
        elif os.path.exists(tgt2):
            sj = GGB.read_bed(tgt2)
        else:
            sj = GGB.read_bed(tgt3)
        idx0 = (sj['chr']==chrom)&(sj['tst']>=st)&(sj['ted']<=ed)        
        self.sjpaths0 = sj[idx0].copy()        
        # load exdf, sjdf
        sjdf = UT.read_pandas(modelpre+'.sjdf.txt.gz', names=A3.SJDFCOLS)
        sjdf['tst'] = sjdf['st'] # for sjpath compatibility
        sjdf['ted'] = sjdf['ed']
        sjdf['sc1'] = sjdf['ucnt']
        sjdf['sc2'] = sjdf['tcnt']
        sjdf = sjdf[(sjdf['chr']==chrom)&(sjdf['st']>=st)&(sjdf['ed']<=ed)]
        sjdf = sjdf[sjdf['name'].isin(sids)]
        self.sjdf = sjdf.groupby(['chr','st','ed','strand']).first().reset_index()

        exdf = UT.read_pandas(modelpre+'.exdf.txt.gz', names=A3.EXDFCOLS)
        exdf = exdf[(exdf['chr']==chrom)&(exdf['st']>=st)&(exdf['ed']<=ed)]
        exdf = exdf[exdf['name'].isin(eids)]
        if os.path.exists(modelpre+'.sedf.txt.gz'):
            sedf = UT.read_pandas(modelpre+'.sedf.txt.gz', names=A3.EXDFCOLS)
            sedf = sedf[(sedf['chr']==chrom)&(sedf['st']>=st)&(sedf['ed']<=ed)]
            sedf = sedf[sedf['name'].isin(eids)]
            exdf = PD.concat([exdf,sedf],ignore_index=True)
        self.exdf = exdf.groupby(['chr','st','ed','strand','kind']).first().reset_index()
        A3.set_ad_pos(self.sjdf, 'sj')
        A3.set_ad_pos(self.exdf, 'ex')
        # filled
        self.filled = {}
        sjs = self.sjdf
        exs = self.exdf
        exs['ost'] = exs['st']-self.st
        exs['oed'] = exs['ed']-self.st
        for s in ['+','-']:
            sja = self.arrs['sj'][s]
            sj = sjs[sjs['strand'].isin(A3.STRS[s])]
            ex = exs[exs['strand'].isin(A3.STRS[s])]
            self.filled[s] = A3.fill_gap(sja, sj, ex, s, self.st)
        fix_i53completematch(self.exdf) # extend 5'3' exons completely matched internal exons


    def process(self):
        self.calculate_ecovs()
        self.calculate_scovs()
        self.estimate_abundance()
        self.write()
        return 

    def calculate_scovs(self):
        sj = self.sjdf
        sj0 = self.sjpaths0
        sj0mat = sj0[['sc1','sc2','name']].values
        tmp = [[(sc1,sc2) for sc1,sc2,p in sj0mat if y in p] for y in sj['name']]
        sj['ucnt'] = [N.sum([x[0] for x in y]) for y in tmp]
        sj['tcnt'] = [N.sum([x[1] for x in y]) for y in tmp]
        self.sjdfi = sj.set_index('name')

    def calculate_branchp(self, jids, eids):
        sj0 = self.sjdfi
        sj = sj0.ix[jids].reset_index()
        ex0 = self.exdfi
        ex = ex0.ix[eids].reset_index()

        dsump = sj.groupby('dpos')['tcnt'].sum().astype(float)
        tmp = dsump.ix[sj['dpos'].values]
        jdp = sj['tcnt'].values/tmp.values
        idx = N.array(tmp==0, dtype=bool)
        jdp[idx] = 0. 
        j2p = dict(zip(sj['name'].values, jdp))
        
        # exon groupby acceptor
        asump = ex.groupby('apos')['ecov'].sum().astype(float)
        tmp = asump.ix[ex['apos'].values]
        eap = ex['ecov'].values/(tmp.values)
        idx = N.array(tmp==0, dtype=bool)
        eap[idx] = 0. 
        e2ap = dict(zip(ex['name'].values, eap))

        dsump = ex.groupby('dpos')['ecov'].sum().astype(float)
        tmp = dsump.ix[ex['dpos'].values]
        edp = ex['ecov'].values/(tmp.values)
        idx = N.array(tmp==0, dtype=bool)
        edp[idx] = 0. 
        e2dp = dict(zip(ex['name'].values, edp))

        return j2p, e2ap, e2dp

    def tcov_by_nnls(self, s, e, strand):
        o = int(self.st)
        p = self.paths
        idx = (p['tst']>=s)&(p['ted']<=e)&(p['strand'].isin(A3.STRS[strand]))
        ps = p[idx]
        if len(ps)==0:
            return None
        pg = ps.groupby(['tst','ted']).first().reset_index()[['chr','tst','ted','strand','name']].sort_values(['tst','ted'])
        pg['strand'] = strand
        ne = len(pg)
        exa = self.arrs['ex'][strand]
        # sja = self.arrs['sj'][strand]
        sja = self.filled[strand]
        def cov0(s,e):
            # return N.sum(sja[s-o:e-o]+exa[s-o:e-o])/(e-s)
            return N.mean(sja[s-o:e-o])
        # def cov1s(s):
        #     s0 = max(0, int(s)-o-10)
        #     s1 = max(s0+1,int(s)-o)
        #     return N.mean(exa[s0:s1])
        # def cov1e(e):
        #     return N.mean(exa[int(e)-o:int(e)-o+10])
        e_ed2cov = self.eed2cov[strand]
        e_st2cov = self.est2cov[strand]
        def cov1s(s):
            return e_ed2cov.get(s,0)
        def cov1e(e):
            return e_st2cov.get(e,0)
        def cov2s(s): # donor
            # s0 = max(0, s-o-1)
            return max(0, sja[int(s)-o]-sja[int(s)-o-1])
        def cov2e(e): # acceptor
            # e0 = max(0, e-o-1)
            return max(0, sja[int(e)-o-1]-sja[int(e)-o])
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
                pg['tcov0a'] = ecov
            except e:
                # too much iteration?
                LOG.warning('!!!!!! Exception in NNLS (tcov_by_nnls) @{0}:{1}-{2}, setting to zero !!!!!!!!!'.format(self.chrom, s, e))
                pg['tcov0a'] = 0
                # raise e
            pg.rename(columns={'st':'tst','ed':'ted'}, inplace=True)
        else: # this includes single exons
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
            if stsum<1e-9 or edsum<1e-9:
                pg['tcov0b'] = 0
            else:
                c0 = c.copy()
                if strand in ['+','.+']:
                    c[:nst] = (edsum/stsum)*c[:nst]
                else:
                    c[nst:] = (stsum/edsum)*c[nst:]
                try:
                    ecov,err = nnls(mat, c)
                except e:
                    print('s:{0},e:{1},strand:{2}'.format(s,e,strand))
                    print('stsum:', stsum)
                    print('edsum:', edsum)
                    print('nnls error tcov0b', mat, c, c0)
                    print('sts:',sts)
                    print('eds:',eds)
                    print('pg:',pg)
                    pg['tcov0c'] = 0
                    raise e

                pg['tcov0b'] = ecov

            mat = N.array([(pg['tst']==x).values for x in sts]+[(pg['ted']==x).values for x in eds], dtype=float)
            c = N.array([cov2s(x) for x in sts]+[cov2e(x) for x in eds])
            # enforce flux conservation: scale up 5'
            stsum = N.sum(c[:nst])
            edsum = N.sum(c[nst:])
            if stsum<1e-9 or edsum<1e-9:
                pg['tcov0c'] = 0
            else:
                if strand in ['+','.+']:
                    c[:nst] = (edsum/stsum)*c[:nst]
                else:
                    c[nst:] = (stsum/edsum)*c[nst:]
                try:
                    ecov,err = nnls(mat, c)
                except e:
                    print('s:{0},e:{1},strand:{2}'.format(s,e,strand))
                    print('nnls error tcov0c', mat, c)
                    pg['tcov0c'] = 0
                    raise e
                pg['tcov0c'] = ecov
        else:
            s,e = pg.iloc[0][['tst','ted']]
            pg['tcov0b'] = (cov1s(s)+cov1e(e))/2.
            pg['tcov0c'] = (cov2s(s)+cov2e(e))/2.

        if not self.usegeom:
            pg['tcov0'] = pg[['tcov0a','tcov0b','tcov0c']].mean(axis=1)
            # pg['tcov0'] = (2*pg['tcov0a']+pg['tcov0b']+pg['tcov0c'])/4. # weighted        
        else:
            pg['tcov0'] = N.power(pg['tcov0a']*pg['tcov0b']*pg['tcov0c'], 1/3.) # geometric mean
        pg.loc[pg['tcov0']<0,'tcov0'] = 0 # shouldn't really happen
        keys = [tuple(x) for x in p[idx][['tst','ted']].values]
        for f in ['tcov0','tcov0a','tcov0b','tcov0c']:
            p.loc[idx, f] = pg.set_index(['tst','ted']).ix[keys][f].values
        return pg[['chr','tst','ted','strand','tcov0']]
                
    def tcov_by_branchp(self, tst, ted, strand, tcov0):
        p = self.paths
        idx = (p['strand'].isin(A3.STRS[strand]))&(p['tst']==tst)&(p['ted']==ted)
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
        idxme = paths['name'].str.contains('\|')
        mepaths = paths[idxme].copy()
        sepaths = paths[~idxme].copy()

        self.paths = mepaths
        for s in ['+','-']:
            ps = mepaths[mepaths['strand'].isin(A3.STRS[s])]
            if len(ps)==0:
                continue
            # for chrom,st,ed in UT.union_contiguous(ps[['chr','st','ed']],returndf=False):
            poscols = ['chr','tst','ted']
            for chrom,st,ed in UT.union_contiguous(ps[poscols],pos_cols=poscols,returndf=False):
                pg = self.tcov_by_nnls(st,ed,s)
                if pg is not None:
                    for chrom,tst,ted,strand,tcov0 in pg.values:
                        self.tcov_by_branchp(tst,ted,strand,tcov0)

        e2c = UT.df2dict(self.exdf, 'name', 'ecov')
        sepaths['tcov'] = [e2c[x] for x in sepaths['name']]
        for f in ['tcov0','tcov0b']:
            sepaths[f] = sepaths['tcov']
        sepaths['tcov0a'] = 0.
        sepaths['tcov0c'] = 0.

        paths = PD.concat([mepaths, sepaths], ignore_index=True)
        paths.sort_values(['chr','st','ed'],inplace=True)
        self.paths = paths

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

def bundle_estimator(modelpre, bwpre, chrom, st, ed, dstpre, tcovth, usegeom):
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
    la = LocalEstimator(modelpre, bwpre, chrom, st, ed, dstpre, tcovth, usegeom)
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


def estimatecovs(modelpre, bwpre, dstpre, genome, tcovth=1, usegeom=True, np=6):
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
            args.append([modelpre, bwpre, chrom, st, ed, dstpre, tcovth, usegeom])
            bundles.append((chrom,st,ed))

    rslts = UT.process_mp(bundle_estimator, args, np=np, doreduce=False)

    concatenate_bundles(bundles, dstpre)


class CovEstimator(object):
    
    def __init__(self, modelpre, bwpre, dstpre, genome, tcovth=1, usegeom=False, np=6):
        self.modelpre = modelpre
        self.bwpre = bwpre
        self.dstpre = dstpre
        self.genome = genome
        self.tcovth = tcovth
        self.usegeom = usegeom
        self.np = np
        
    def run(self):
        self.server = server = TQ.Server(np=self.np)
        print('reading paths.withse.bed.gz')
        bed = GGB.read_bed(self.modelpre+'.paths.withse.bed.gz')
        chroms = bed['chr'].unique()
        csizedic = UT.df2dict(UT.chromdf(self.genome), 'chr', 'size')
        self.bundlestatus = bundlestatus = {}
        self.bundles = bundles = []

        with server:
            print('starting task server')
            subid = 0
            for chrom in chroms:
                print('chrom {0}'.format(chrom))
                sub = bed[(bed['chr']==chrom)]
                uc = UT.union_contiguous(sub[['chr','st','ed']], returndf=True)
                # total about 30K=> make batch of ~1000
                n = len(uc)
                nb = int(N.ceil(n/1000.))
                print(chrom,nb)
                for i in range(nb):
                    print('putting in bundle_estimator {0}.{1}'.format(chrom,subid))
                    sti = 1000*i
                    edi = min(1000*(i+1), len(uc)-1)
                    st = max(uc.iloc[sti]['st'] - 100, 0)
                    ed = min(uc.iloc[edi]['ed'] + 100, csizedic[chrom])
                    args = [self.modelpre, self.bwpre, chrom, st, ed, self.dstpre, self.tcovth, self.usegeom]
                    tname = 'bundle_estimator.{0}'.format(subid)
                    subid += 1
                    task = TQ.Task(tname, bundle_estimator, args)
                    server.add_task(task)
                    bundles.append((chrom,st,ed))
            nb = len(bundles)
            while server.check_error():
                try:
                    name, rslt = server.get_result(timeout=5)
                except TQ.Empty:
                    name, rslt = None, None
                if name is not None:
                    if name.startswith('bundle_estimator.'):
                        subid = name.split('.')[-1]
                        bundlestatus[subid] = rslt
                        if len(bundlestatus)==nb:
                            print('$$$$$$$$ putting in concatenate_bundles $$$$$$$$$$$')
                            tname='concatenate_bundles'
                            args = (bundles, self.dstpre)
                            task = TQ.Task(tname, concatenate_bundles, args)
                            server.add_task(task)
                    if name=='concatenate_bundles':
                        print('$$$$$$$$ concatenate_bundles done $$$$$$$$$$$')
                        break
            print('Exit Loop')
        print('Done')


############# Cov Collector ######################################################

class CovCollector(object):
    
    def __init__(self, covpres, dstpre, np=7):
        self.covpres = covpres
        self.modelpre = covpres[0]
        self.dstpre = dstpre
        self.np = np
        
    def run(self):
        self.server = server = TQ.Server(np=self.np)
        self.exdf = ex = UT.read_pandas(self.modelpre+'.covs.exdf.txt.gz', names=A3.EXDFCOLS)
        self.chroms = chroms = ex['chr'].unique()
        self.exstatus = exstatus = {}
        self.sjstatus = sjstatus = {}
        self.pastatus = pastatus = {}
        exdone=False
        sjdone=False
        padone=False
        n = len(self.covpres)
        nb = int(N.ceil(n/50.))
        with server:
            for subid in range(nb):
                covpressub = self.covpres[50*subid:50*(subid+1)]
                # ex
                tname = 'collect_ecov_subset.{0}'.format(subid)
                args = (self.modelpre, covpressub, self.dstpre, subid)
                task = TQ.Task(tname, collect_ecov_subset, args)
                server.add_task(task)
                # sj
                tname = 'collect_tcnt_subset.{0}'.format(subid)
                args = (self.modelpre, covpressub, self.dstpre, subid)
                task = TQ.Task(tname, collect_tcnt_subset, args)
                server.add_task(task)
                # path
                tname = 'collect_tcovs_subset.{0}'.format(subid)
                args = (self.modelpre, covpressub, self.dstpre, subid)
                task = TQ.Task(tname, collect_tcovs_subset, args)
                server.add_task(task)
            while server.check_error():
                try:
                    name, rslt = server.get_result(timeout=5)
                except TQ.Empty:
                    name, rslt = None, None
                if name is not None:
                    if name.startswith('collect_ecov_subset.'):
                        subid = name.split('.')[-1]
                        exstatus[subid] = rslt
                        if len(exstatus)==nb:
                            print('$$$$$$$$ putting in concatenate_ecov_subsets $$$$$$$$$$$')
                            for chrom in chroms:
                                tname='concatenate_ecov_subsets'
                                args = (self.modelpre, self.dstpre, range(nb), chrom)
                                task = TQ.Task(tname, concatenate_ecov_subsets, args)
                                server.add_task(task)
                    if name.startswith('collect_tcnt_subset.'):
                        subid = name.split('.')[-1]
                        sjstatus[subid] = rslt
                        if len(sjstatus)==nb:
                            print('$$$$$$$$ putting in concatenate_tcnt_subsets $$$$$$$$$$$')
                            for chrom in chroms:
                                tname='concatenate_tcnt_subsets'
                                args = (self.modelpre, self.dstpre, range(nb), chrom)
                                task = TQ.Task(tname, concatenate_tcnt_subsets, args)
                                server.add_task(task)
                    if name.startswith('collect_tcovs_subset.'):
                        subid = name.split('.')[-1]
                        pastatus[subid] = rslt
                        if len(pastatus)==nb:
                            print('$$$$$$$$ putting in concatenate_tcovs_subsets $$$$$$$$$$$')
                            for chrom in chroms:
                                tname='concatenate_tcovs_subsets'
                                args = (self.modelpre, self.dstpre, range(nb), chrom)
                                task = TQ.Task(tname, concatenate_tcovs_subsets, args)
                                server.add_task(task)
                    if name=='concatenate_ecov_subsets':
                        print('$$$$$$$$ concatenate_ecov_subsets done $$$$$$$$$$$')
                        exdone=True
                    if name=='concatenate_tcnt_subsets':
                        print('$$$$$$$$ concatenate_tcnt_subsets done $$$$$$$$$$$')
                        sjdone=True
                    if name=='concatenate_tcovs_subsets':
                        print('$$$$$$$$ concatenate_tcovs_subsets done $$$$$$$$$$$')
                        padone=True
                    if exdone&sjdone&padone:
                        break
            print('Exit Loop')
        print('Done')
                    
                
        
def collect_ecov_subset(modelpre, covpressub, dstpre, subid):
    return _collect_subset(modelpre, covpressub, dstpre, subid, 'ex')

def concatenate_ecov_subsets(modelpre, dstpre, subids, chrom):
    return _concatenate_subsets(modelpre, dstpre, subids, 'ex', chrom)

def collect_tcnt_subset(modelpre, covpressub, dstpre, subid):
    return _collect_subset(modelpre, covpressub, dstpre, subid, 'sj')

def concatenate_tcnt_subsets(modelpre, dstpre, subids, chrom):
    return _concatenate_subsets(modelpre, dstpre, subids, 'sj', chrom)
    
def collect_tcovs_subset(modelpre, covpressub, dstpre, subid):
    return _collect_subset(modelpre, covpressub, dstpre, subid, 'pa')

def concatenate_tcovs_subsets(modelpre, dstpre, subids, chrom):
    return _concatenate_subsets(modelpre, dstpre, subids, 'pa', chrom)

def _collect_subset(modelpre, covpressub, dstpre, subid, which):
    if which == 'ex':
        suf = 'exdf'
        flds = ['ecov']
        fsuf = 'ecovs'
        cols = A3.EXDFCOLS
    elif which == 'sj':
        suf = 'sjdf'
        flds = ['tcnt']
        fsuf = 'tcnts'
        cols = A3.SJDFCOLS
    else:
        suf = 'paths'
        flds = ['tcov0','tcov']
        fsuf = 'tcovs'
        cols = A3.PATHCOLS
    ex0 = UT.read_pandas(modelpre+'.covs.{0}.txt.gz'.format(suf), names=cols)
    chroms = ex0['chr'].unique()
    # read in exdf sort, transpose and write(append) to dst
    if all([os.path.exists(dstpre+'.{1}.{0}.txt.gz'.format(c,fsuf)) for c in chroms]):
        return []
    if all([os.path.exists(dstpre+'.{2}.{0}.{1}.txt.gz'.format(c,subid,fsuf)) for c in chroms]):
        return []
    ex0.sort_values(['chr','st','ed','strand'], inplace=True)
    names = []
    for pre in covpressub:
        name = pre.split('/')[-1]
        ex1 = UT.read_pandas(pre+'.covs.{0}.txt.gz'.format(suf), names=cols)
        ex1.sort_values(['chr','st','ed','strand'], inplace=True)
        for f in flds:
            cname = '{0}.{1}'.format(name, f)
            ex0[cname] = ex1[f].values
            names.append(cname)
    ex0.reset_index(inplace=True)
    files = []
    for chrom in ex0['chr'].unique():
        ex0chr = ex0[ex0['chr']==chrom].sort_values(['st','ed','strand'])
        dst = dstpre+'.{2}.{0}.{1}.txt.gz'.format(chrom,subid,fsuf)
        UT.write_pandas(ex0chr[names].T, dst, 'i')
        files.append(dst)
    return files

def _concatenate_subsets(modelpre, dstpre, subids, which, chrom):
    if which == 'ex':
        suf = 'exdf'
        fsuf = 'ecovs'
        cols = A3.EXDFCOLS
    elif which == 'sj':
        suf = 'sjdf'
        fsuf = 'tcnts'
        cols = A3.SJDFCOLS
    else:
        suf = 'paths'
        fsuf = 'tcovs'
        cols = A3.PATHCOLS
    
    ex0 = UT.read_pandas(modelpre+'.covs.{0}.txt.gz'.format(suf), names=cols)
    chroms = ex0['chr'].unique()
    files = []
    dstpath0 = dstpre+'.{1}.{0}.tmp.txt.gz'.format(chrom,fsuf)
    dstpath1 = dstpre+'.{1}.{0}.txt.gz'.format(chrom,fsuf)
    if not os.path.exists(dstpath1):
        with open(dstpath0, 'wb') as dst:
            for subid in subids:
                srcpath = dstpre+'.{2}.{0}.{1}.txt.gz'.format(chrom,subid,fsuf)
                with open(srcpath, 'rb') as src:
                    shutil.copyfileobj(src,dst)
                files.append(srcpath)
        ex0chr = ex0[ex0['chr']==chrom].sort_values(['st','ed','strand'])
        ex1chr = UT.read_pandas(dstpath0,names=ex0chr.index,index_col=[0]).T
        df = PD.concat([ex0chr, ex1chr],axis=1)
        UT.write_pandas(df, dstpath1, 'h')
        #os.unlink(dstpath0)
    for f in files:
        if os.path.exists(f):
            os.unlink(f)
    return dstpath1
            


def fix_i53completematch(exdf):
    # extend edge of 5'3' exons if they completely match to internal exons
    idxp = exdf['strand'].isin(A3.STRS['+'])
    idx5 = exdf['kind']=='5'
    idx3 = exdf['kind']=='3'
    idxi = exdf['kind']=='i'
    ileft = (idxp&idx5)|(~idxp&idx3)
    iright = (idxp&idx3)|(~idxp&idx5)
    steds = set([(x,y) for x,y in exdf[idxi][['st','ed']].values])

    idxm = N.array([(x,y) in steds for x,y in exdf[['st','ed']].values], dtype=bool)
    imleft = ileft&idxm
    imright = iright&idxm
    while (N.sum(imleft)+N.sum(imright))>0:
        exdf.loc[imleft,'st'] = exdf[imleft]['st']-10
        exdf.loc[imright, 'ed'] = exdf[imright]['ed']+10
        idxm = N.array([(x,y) in steds for x,y in exdf[['st','ed']].values], dtype=bool)
        imleft = ileft&idxm
        imright = iright&idxm











