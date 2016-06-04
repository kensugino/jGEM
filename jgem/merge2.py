"""

.. module:: merge2
    :synopsis: merge assemblies from different cell types
    jGEM version 2 merger

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
from jgem import assembler2 as A2
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
        exdone=False
        sjdone=False
        padone=False
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
                    if name=='prep_exbw':
                        print('$$$$$$$$ prep_exbw done $$$$$$$$$$$')
                        exdone=True
                    if name=='prep_sjbw':
                        print('$$$$$$$$ prep_sjbw done $$$$$$$$$$$')
                        sjdone=True
                    if name=='prep_sjpath':
                        print('$$$$$$$$ prep_sjpath done $$$$$$$$$$$')
                        padone=True
                    if exdone&sjdone&padone:
                        break
            print('Exit Loope')
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
        exdf = UT.read_pandas(pre+'.exdf.txt.gz',names=A2.EXDFCOLS)
        exdf = exdf[exdf['chr']==chrom]
        for s in ss:
            exsub = exdf[exdf['strand'].isin(s2s[s])]
            for st,ed,ecov in exsub[['st','ed','ecov']].values:
                a[s][st:ed] += ecov*scale
        sedf = UT.read_pandas(pre+'.sedf.txt.gz',names=A2.EXDFCOLS)
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
        sjdf = UT.read_pandas(pre+'.sjdf.txt.gz',names=A2.SJDFCOLS)
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
        paths = UT.read_pandas(pre+'.paths.txt.gz', names=A2.PATHCOLS)
        paths = paths[paths['chr']==chrom]
        for st,ed,name,s,tst,ted,tcov in paths[cols].values:
            pc = ','.join(name.split(',')[1:-1]) # trim 53exons => intron chain
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
    cmax = 9+N.log2(N.mean(scales))
    bed = A2.path2bed12(df, cmax)
    # reset sc1 to tcov (from log2(tcov+2)*100)
    bed['sc1'] = bed['tcov']
    GGB.write_bed(bed, path, ncols=12)
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

