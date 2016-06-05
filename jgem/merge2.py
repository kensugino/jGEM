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
    df = df.groupby('pc').first() # get rid of unstranded duplicates
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



SJFILTERPARAMS = dict(
    th_detected=1,
    th_maxcnt=1,
    th_maxoverhang=18,
    th_minedgeexon=18,
    th_sjratio2=1e-3,
)
class SJFilter(object):

    def __init__(self, bwsjpre, statspath, genome, np=10, *kw):
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
        rslts = UT.process_mp(filter_sj, args, np=self.np, doreduce=False)

        dstpath = self.bwsjpre+'.sjpath.filtered.bed.gz'
        with open(dstpath,'wb') as dst:
            for c in chroms:
                srcpath =  self.bwsjpre+'.sjpath.{0}.filtered.bed.gz'.format(c)
                with open(srcpath, 'rb') as src:
                    shutil.copyfileobj(src, dst)

def locus2pc(l):
    chrom,sted,strand = l.split(':')
    st,ed = sted.split('-')
    st = str(int(st)-1)
    if strand in ['+','.']:
        return  '|'.join([st,ed])
    return '|'.join([ed,st])

def filter_sj(bwsjpre, statspath, chrom, csize, params):
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
    # calculate sjratio, sjratio2
    sjexbw = A2.SjExBigWigs(bwsjpre, mixunstranded=False)
    for s in ['+','-']:
        idx = sj['strand']==s
        with sjexbw:
            sa = sjexbw.bws['sj'][s].get(chrom,0,csize)
            ea = sjexbw.bws['ex'][s].get(chrom,0,csize)
        a = sa+ea
        sj.loc[idx,'sjratio2'] = [x/N.mean(a[int(s):int(e)]) for x,s,e in sj[idx][['sc1','tst','ted']].values]
    sj = sj[sj['sjratio2']>params['th_sjratio2']]
    GGB.write_bed(sj, dstpath, ncols=12)




class LocalEstimator(A2.LocalAssembler):

    def __init__(self, bed12path, bwpre, chrom, st, ed, dstpre, tcovth):
        self.bed12path = bed12path
        self.tcovth = tcovth
        A2.LocalAssembler.__init__(self, bwpre, chrom, st, ed, dstpre, refcode=None)
        bed12 = GGB.read_bed(bed12path)
        idx = (bed12['chr']==chrom)&(bed12['tst']>=st)&(bed12['ted']<=ed)
        self.paths = bed12[idx].copy()
        sj = GGB.read_bed(bwpre+'.sjpath.bed.gz')
        idx0 = (sj['chr']==chrom)&(sj['tst']>=st)&(sj['ted']<=ed)        
        self.sjpaths0 = sj[idx0].copy()        

    def process(self):
        self.make_sjexdf()
        self.calculate_ecovs()
        self.calculate_scovs()
        self.estimate_abundance()
        self.write()

    def make_sjexdf(self):
        ap = self.paths
        dfs = []
        dfe = []
        chrom = self.chrom
        # allpaths name = e5,e5d|a,d|a,...,d|e3a,e3
        def _sgen():
            sted = set()
            for p,strand in ap[['name','strand']].values:
                for x in p.split(',')[1:-1]:
                    st,ed = [int(y) for y in x.split('|')]
                    if st>ed:
                        st,ed = ed,st
                    if (st,ed,strand) not in sted:
                        yield (chrom,st,ed,strand,x,'j')
                        sted.add((st,ed,strand))
        def _egen():
            sted = set()
            for p,strand in ap[['name','strand']].values:
                tmp = p.split('|')
                # 53
                st,ed = [int(y) for y in tmp[0].split(',')]
                if st>ed:
                    st,ed = ed,st
                if (st,ed,strand) not in sted:
                    yield (chrom,st,ed,strand,tmp[0],'5')
                    sted.add((st,ed,strand))
                st,ed = [int(y) for y in tmp[-1].split(',')]
                if st>ed:
                    st,ed = ed,st
                if (st,ed,strand) not in sted:
                    yield (chrom,st,ed,strand,tmp[-1],'3')
                    sted.add((st,ed,strand))
                # internal
                for x in tmp[1:-1]:
                    st,ed = [int(y) for y in x.split(',')]
                    if st>ed:
                        st,ed = ed,st
                    if (st,ed,strand) not in sted:
                        yield (chrom,st,ed,strand,x,'i')
                        sted.add((st,ed,strand))
        cols = ['chr','st','ed','strand','name','kind']
        sjdf = PD.DataFrame([x for x in _sgen()], columns=cols)
        exdf = PD.DataFrame([x for x in _egen()], columns=cols)
        A2.set_ad_pos(sjdf, 'sj')
        A2.set_ad_pos(exdf, 'ex')
        self.sjdf = sjdf
        self.exdf = exdf
        
    def estimate_abundance(self):
        # 1) 5-3 group by NNLS
        # 2) UTR difference
        # 3) within 5-3 group by tree branch prob
        paths = self.paths
        for s in ['+','-']:
            ps = paths[paths['strand'].isin(A2.STRS[s])]
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

    def write(self):
        pre = self.dstpre+'.{0}_{1}_{2}'.format(self.chrom,self.st,self.ed)
        # 1) exon, junctions, allpaths => csv (no header <= to concatenate bundles)
        ecols = A2.EXDFCOLS #['chr','st','ed','strand','name','kind','ecov']
        UT.write_pandas(self.exdf[ecols], pre+'.covs.exdf.txt.gz', '')
        scols = A2.SJDFCOLS #['chr','st','ed','strand','name','kind','tcnt'  ]#,'donor','acceptor','dp','ap']
        UT.write_pandas(self.sjdf[scols], pre+'.covs.sjdf.txt.gz', '')
        pcols = A2.PATHCOLS #['chr','st','ed','name','strand','tst','ted','tcov0','tcov1','tcov']
        UT.write_pandas(self.paths[pcols], pre+'.covs.paths.txt.gz', '')
        # write colored bed12 for tcov > th
        tgt = self.paths[self.paths['tcov']>=self.tcovth].copy()
        self.bed12 = A2.path2bed12(tgt, cmax=9, covfld='tcov')
        GB.write_bed(self.bed12, pre+'.covs.paths.bed.gz',ncols=12)

def bundle_estimator(bed12path, bwpre, chrom, st, ed, dstpre, tcovth):
    bname = A2.bundle2bname((chrom,st,ed))
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
    la = LocalEstimator(bed12path, bwpre, chrom, st, ed, dstpre, tcovth)
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
                    bname = A2.bundle2bname((chrom,st,ed))
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


def estimatecovs(bed12path, bwpre, dstpre, genome, tcovth=1, np=6):
    bed = GGB.read_bed(bed12path)
    chroms = bed['chr'].unique()
    csizedic = UT.df2dict(UT.chromdf(genome), 'chr', 'size')
    bundles = []
    for strand in ['+','-','.']:
        for chrom in chroms:
            sub = bed[(bed['chr']==chrom)&(bed['strand']==strand)]
            uc = UT.union_contiguous(sub[['chr','st','ed']], returndf=True)
            # total about 30K=> make batch of ~1000
            n = len(uc)
            nb = int(N.ceil(n/1000.))
            args = []
            for i in range(nb):
                sti = 1000*i
                edi = min(1000*(i+1), len(uc)-1)
                st = max(uc.iloc[sti]['st'] - 100, 0)
                ed = min(uc.iloc[edi]['ed'] + 100, csizedic[chrom])
                args.append([bed12path, bwpre, chrom, st, ed, dstpre, tcovth])
                bundles.append((chrom,st,ed))

    rslts = UT.process_mp(bundle_estimator, args, np=np, doreduce=False)
    concatenate_bundles(bundles, dstpre)




