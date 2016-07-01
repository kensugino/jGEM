"""

.. module:: convert
    :synopsis: convert sj,ex to bed, gtf  or gtf to sj,ex

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
import os
import tempfile
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)
import uuid
try:
    from itertools import izip
except:
    izip = zip

import pandas as PD
import numpy as N

from jgem import utils as UT
from jgem import gtfgffbed as GGB
from jgem import graph as GP

from jgem import assembler2 as A2
from jgem import assembler3 as A3

# Path(Bed12) => GTF

def gtf_from_bed12(modelpre, dstpath=None, source='.'):
    # path['gname'] contains gene id
    paths = GGB.read_bed(modelpre+'.paths.withse.bed.gz')
    ex = UT.read_pandas(modelpre+'.ex.txt.gz')
    ex['id'] = ex['chr']+':'+ex['name']
    n2gn = UT.df2dict(ex, 'id', 'gname')
    # n2gn = UT.df2dict(ex, 'name', 'gname') # there may be same st,ed in different chromosome
    paths['id'] = paths['chr']+':'+paths['name']
    paths['id0'] = paths['chr']+':'+paths['name'].str.split('|').str[0]
    paths['gname'] = [n2gn[x] for x in paths['id0']]
    g2cnt = {}
    tnames = []
    for x in paths['gname']:
        i = g2cnt.get(x,1)
        tnames.append('{0}.{1}'.format(x,i))
        g2cnt[x] = i+1
    paths['tname'] = tnames    
    txt = 'gene_id "{0}"; transcript_id "{1}"; exon_number "{2}";'
    def _gen():
        cols = ['chr','st','ed','gname','tname','esizes','estarts','strand']
        for c,s,e,gn,tn,esi,est,strand in paths[cols].values:
            esizes = [int(x) for x in esi.split(',')[:-1]]
            estarts = [int(x) for x in est.split(',')[:-1]]
            for i,(x,y) in enumerate(zip(esizes,estarts)):
                est = s+y
                eed = est+x
                extra = txt.format(gn,tn,i+1)
                yield (c,source,'exon',est+1,eed,'.',strand,'.',extra)
    df = PD.DataFrame([x for x in _gen()], columns=GGB.GTFCOLS)
    if dstpath is None:
        dstpath = bedpath.replace('.bed','.gtf')
    GGB.write_gtf(df, dstpath)
    
    idf = paths[['id','chr','name','tname','gname']]
    UT.write_pandas(idf, mdstpre+'.idmap.txt.gz','h')
    return df


# GTF <=> EX,SJ    ######################################################################


def as2exsj(dstpre, np=7):
    ex = UT.read_pandas(dstpre+'.exdf.txt.gz', names=A2.EXDFCOLS)
    sj = UT.read_pandas(dstpre+'.sjdf.txt.gz', names=A2.SJDFCOLS)
    se = UT.read_pandas(dstpre+'.sedf.txt.gz', names=A2.EXDFCOLS)
    paths = UT.read_pandas(dstpre+'.paths.txt.gz', names=A2.PATHCOLS)
    #ex.loc[ex['strand'].isin(['.+','.-']),'strand'] = '.'
    #sj.loc[sj['strand'].isin(['.+','.-']),'strand'] = '.'
    sj['st'] = sj['st']+1 
    cols = A2.EXDFCOLS
    ex = PD.concat([ex[cols],se[cols]],ignore_index=True)
    UT.set_info(sj,ex)
    UT.set_exon_category(sj, ex)

    # find genes (connected components) set '_gidx'
    graphpre = dstpre+str(uuid.uuid4())+'_'
    prefix = os.path.abspath(graphpre) # need unique prefix for parallel processing
    genes = GP.find_genes4(sj,ex,
        filepre=prefix,
        np=np,
        override=False,
        separatese=True)
    ex.loc[ex['kind']=='3','cat'] = '3'
    ex.loc[ex['kind']=='5','cat'] = '5'
    UT.write_pandas(ex, dstpre+'.ex.txt.gz', 'h')
    UT.write_pandas(sj, dstpre+'.sj.txt.gz', 'h')
    ci = UT.chopintervals(ex, dstpre+'.ci.txt.gz')
    return sj, ex

def as3exsj(dstpre, minelen=150, np=7):
    ex = UT.read_pandas(dstpre+'.exdf.txt.gz', names=A3.EXDFCOLS)
    sj = UT.read_pandas(dstpre+'.sjdf.txt.gz', names=A3.SJDFCOLS)
    se = UT.read_pandas(dstpre+'.sedf.txt.gz', names=A3.EXDFCOLS)
    paths = UT.read_pandas(dstpre+'.paths.txt.gz', names=A3.PATHCOLS)
    #ex.loc[ex['strand'].isin(['.+','.-']),'strand'] = '.'
    # sj.loc[sj['strand']=='.+','strand'] = '+'
    # sj.loc[sj['strand']=='.-','strand'] = '-'
    sj['st'] = sj['st']+1 
    cols = A3.EXDFCOLS
    ex = PD.concat([ex[cols],se[cols]],ignore_index=True)
    UT.set_info(sj,ex)
    UT.set_exon_category(sj, ex)

    # find genes (connected components) set '_gidx'
    graphpre = dstpre+str(uuid.uuid4())+'_'
    prefix = os.path.abspath(graphpre) # need unique prefix for parallel processing
    # genes = GP.find_genes4(sj,ex,
    #     filepre=prefix,
    #     np=np,
    #     override=False,
    #     separatese=True)
    genes = GP.find_genes3(sj,ex, # don't use exon overlap as connection
        np=np,
        override=False)
    ex.loc[ex['kind']=='3','cat'] = '3'
    ex.loc[ex['kind']=='5','cat'] = '5'

    # remove these with elen smaller than minelen
    ex['len'] = ex['ed']-ex['st']
    exsiz = ex.groupby('_gidx')['len'].sum()
    rgidx = exsiz[exsiz<minelen].index.values
    LOG.info('minelen filter #ex {0}=>{1}'.format(len(ex), len(ex)-len(rgidx)))
    ex2 = ex[~ex['_gidx'].isin(rgidx)]
    sj2 = sj[~sj['_gidx'].isin(rgidx)]
    # write
    UT.write_pandas(ex2, dstpre+'.ex.txt.gz', 'h')
    UT.write_pandas(sj2, dstpre+'.sj.txt.gz', 'h')
    ci = UT.chopintervals(ex2, dstpre+'.ci.txt.gz')
    return sj2, ex2

def as3genes(dstpre):
    # assume converted by as3exsj
    sj = UT.read_pandas(dstpre+'.sj.txt.gz')
    ex = UT.read_pandas(dstpre+'.ex.txt.gz')
    gpre = dstpre+'.genegraph-'
    mg = GP.MEGraph4(sj,ex,gpre)
    bw = BEDWriter(mg)
    path = dstpre+'.genes.bed.gz'
    LOG.info(' writing bed12 genes... {0}'.format(path))
    bw.write(path)

def gtf2exonsj(gtf, np=12, graphpre=None):
    """Extract exons and sj from GTF
    exon, junction coordinates = zero based (same as BED)
    junction start-1 = exon end
    junction end = exon start

    Args:
        gtf: Pandas.DataFrame 

    Returns:
        sj, ex: Pandas.DataFrames for splice junctions and exons

    """
    if len(gtf)==0: # edge case
        cols = GGB.BEDCOLS[:6]+['locus','_id','cat']
        sj = UT.make_empty_df(cols)
        ex = UT.make_empty_df(cols)
        return sj,ex
    exons = gtf[gtf['typ']=='exon'].sort_values(['chr','st','ed'])
    exons['_id'] = N.arange(len(exons))
    exons.sort_values(['transcript_id','st','ed'],inplace=True)
    # 5',3'
    ex_s = exons.groupby('transcript_id').size()
    tid_s = ex_s[ex_s==1].index
    id_m = ex_s[ex_s>1].index
    ex_m = exons[exons['transcript_id'].isin(id_m)]
    ex_m.sort_values(['transcript_id','st','ed'], inplace=True)
    ex_f = ex_m.groupby('transcript_id').first()
    ex_l = ex_m.groupby('transcript_id').last()
    if5 = list(ex_f[ex_f['strand']=='+']['_id'].values)
    if3 = list(ex_f[ex_f['strand']=='-']['_id'].values)
    il5 = list(ex_l[ex_l['strand']=='-']['_id'].values)
    il3 = list(ex_l[ex_l['strand']=='+']['_id'].values)
    exons['kind'] = 'i'
    exons.loc[exons['transcript_id'].isin(tid_s),'kind'] = 's'
    exons.loc[exons['_id'].isin(if5+il5),'kind'] = '5'
    exons.loc[exons['_id'].isin(if3+il3),'kind'] = '3'

    # find junctions
    def _igen():
        for k, g in exons.groupby('transcript_id'):
            if len(g)<2:
                continue
            g = g.sort_values(['st','ed'])
            chrom,strand,gid=g.iloc[0][['chr','strand','gene_id']]
            ists = g['ed'].values[:-1] + 1
            ieds = g['st'].values[1:] - 1
            for st,ed in izip(ists,ieds):
                # chr,st,ed,name=tid,sc1,strand,gene_id
                yield (chrom,st,ed,gid,0,strand)
    sj = PD.DataFrame([x for x in _igen()], columns=GGB.BEDCOLS[:6])
    sj['locus'] = UT.calc_locus_strand(sj)
    sj = sj.groupby('locus').first().reset_index()

    #cols = ['chr','st','ed','gene_id','sc1','strand']
    ex = exons #[cols]
    ex['locus'] = UT.calc_locus_strand(ex)
    ex = ex.groupby(['locus','kind']).first().reset_index()
    ex['name'] = ex['gene_id']
    ex['st'] = ex['st'] - 1

    if len(sj)==0:
        ex['_gidx'] = N.arange(len(ex))
        return sj, ex

    UT.set_info(sj,ex)
    UT.set_exon_category(sj, ex)

    # find genes (connected components) set '_gidx'
    if graphpre is None:
        graphpre = './'+str(uuid.uuid4())+'_'
    prefix = os.path.abspath(graphpre) # need unique prefix for parallel processing
    genes = GP.find_genes4(sj,ex,
        filepre=prefix,
        np=np,
        override=False,
        separatese=True)

    return sj, ex

# def gtf2exonsj(gtf, np=1, graphpre=None):
#     """Extract exons and sj from GTF
#     exon, junction coordinates = zero based (same as BED)
#     junction start-1 = exon end
#     junction end = exon start

#     Args:
#         gtf: Pandas.DataFrame 

#     Returns:
#         sj, ex: Pandas.DataFrames for splice junctions and exons

#     """
#     if len(gtf)==0: # edge case
#         cols = GGB.BEDCOLS[:6]+['locus','_id','cat']
#         sj = UT.make_empty_df(cols)
#         ex = UT.make_empty_df(cols)
#         return sj,ex
#     exons = gtf[gtf['typ']=='exon'].copy()
#     # 5',3'
#     exons.sort_values(['transcript_id','exon_number'])
#     ex_f = exons.groupby('transcript_id').first()
#     ex_l = exons.groupby('transcript_id').last()
#     if5 = ex_f[ex_f['strand']=='+'].index
#     if3 = ex_f[ex_f['strand']=='-'].index
#     il5 = ex_l[ex_l['strand']=='-'].index
#     il3 = ex_l[ex_l['strand']=='+'].index
#     ex_s = exons.groupby('transcript_id').size()
#     tid_s = ex_s[ex_s==1].index
#     exons['kind'] = 'i'
#     exons.loc[exons['transcript_id'].isin(tid_s),'kind'] = 's'
#     exons.loc[exons['transcript_id'].isin(list(if5)+list(il5)),'kind'] = '5'
#     exons.loc[exons['transcript_id'].isin(list(if3)+list(il3)),'kind'] = '3'

#     # find junctions
#     def _igen():
#         for k, g in exons.groupby('transcript_id'):
#             if len(g)<2:
#                 continue
#             g = g.sort_values(['st','ed'])
#             chrom,strand,gid=g.iloc[0][['chr','strand','gene_id']]
#             ists = g['ed'].values[:-1] + 1
#             ieds = g['st'].values[1:] - 1
#             for st,ed in izip(ists,ieds):
#                 # chr,st,ed,name=tid,sc1,strand,gene_id
#                 yield (chrom,st,ed,gid,0,strand)
#     sj = PD.DataFrame([x for x in _igen()], columns=GGB.BEDCOLS[:6])
#     sj['locus'] = UT.calc_locus_strand(sj)
#     sj = sj.groupby('locus').first().reset_index()

#     #cols = ['chr','st','ed','gene_id','sc1','strand']
#     ex = exons #[cols]
#     ex['locus'] = UT.calc_locus_strand(ex)
#     ex = ex.groupby('locus').first().reset_index()
#     ex['name'] = ex['gene_id']
#     ex['st'] = ex['st'] - 1

#     if len(sj)==0:
#         ex['_gidx'] = N.arange(len(ex))
#         return sj, ex

#     UT.set_info(sj,ex)
#     UT.set_exon_category(sj, ex)

#     # find genes (connected components) set '_gidx'
#     if graphpre is None:
#         graphpre = './'+str(uuid.uuid4())+'_'
#     prefix = os.path.abspath(graphpre) # need unique prefix for parallel processing
#     genes = GP.find_genes4(sj,ex,
#         filepre=prefix,
#         np=np,
#         override=False,
#         separatese=True)

#     return sj, ex

# def bed2exonsj(bed12, np=4, graphpre=None):
#     """Extract exons and junctions from BED12

#     Args:
#         bed12: Pandas.DataFrame containing BED12 data

#     Returns:
#         sj, ex: Pandas.DataFrames containing junction and exons

#     """
#     esizes = bed12['esizes'].apply(lambda x: N.array([int(y) for y in x.split(',') if y]))
#     estarts0 = bed12['estarts'].apply(lambda x: N.array([int(y) for y in x.split(',') if y]))
#     bed12['_estarts'] = bed12['st'] + estarts0
#     bed12['_eends'] = bed12['_estarts']+esizes
#     #istarts = eends[:-1]
#     #iends = estarts[1:]
#     cols =['chr','st','ed','tname','strand']
#     def _egen():
#         for chrom,tname,strand,est,eed in UT.izipcols(bed12,['chr','name','strand','_estarts','_eends']):
#             #for st,ed in izip(est,eed):
#             for st,ed in izip(est,eed):
#                 yield (chrom,st,ed,tname,0,strand)
#     def _igen():
#         for chrom,tname,strand,est,eed in UT.izipcols(bed12,['chr','name','strand','_estarts','_eends']):
#             #for st,ed in izip(eed[:-1],est[1:]):
#             for st,ed in izip(eed[:-1],est[1:]):
#                 yield (chrom,st+1,ed,tname,0,strand)
#                 # add 1 to match STAR SJ.tab.out 
#     ex = PD.DataFrame([x for x in _egen()], columns=GGB.BEDCOLS[:6])
#     ex['locus'] = UT.calc_locus_strand(ex)
#     ex = ex.groupby('locus').first().reset_index()
#     sj = PD.DataFrame([x for x in _igen()], columns=GGB.BEDCOLS[:6])
#     sj['locus'] = UT.calc_locus_strand(sj)
#     sj = sj.groupby('locus').first().reset_index()

#     UT.set_info(sj,ex)
#     UT.set_exon_category(sj, ex)

#     # find genes (connected components) set '_gidx'
#     if graphpre is None:
#         graphpre = './'+str(uuid.uuid4())+'_'
#     prefix = os.path.abspath(graphpre) # need unique prefix for parallel processing
#     genes = GP.find_genes4(sj,ex,
#         filepre=prefix,
#         np=np,
#         override=False,
#         separatese=True)

#     return sj, ex
    
# def kg2exonsj(kg, np=4, graphpre=None):
#     """Extract exons and junctions from UCSC knownGene.txt.gz 

#     Args:
#         kg: Pandas.DataFrame containing UCSC knownGene data

#     Returns:
#         sj, ex: Pandas.DataFrames containing junction and exons

#     """
#     kg['_ests'] = kg['exstarts'].apply(lambda x: N.array([int(y) for y in x.split(',') if y]))
#     kg['_eeds'] = kg['exends'].apply(lambda x: N.array([int(y) for y in x.split(',') if y]))
    
#     cols =['chr','st','ed','tname','strand']
#     def _egen():
#         for chrom,tname,strand,est,eed in UT.izipcols(kg,['chr','name','strand','_ests','_eeds']):
#             for st,ed in izip(est,eed):
#                 yield (chrom,st,ed,tname,0,strand)
#     def _igen():
#         for chrom,tname,strand,est,eed in UT.izipcols(kg,['chr','name','strand','_ests','_eeds']):
#             for st,ed in izip(eed[:-1],est[1:]):
#                 yield (chrom,st+1,ed,tname,0,strand)
#                 # add 1 to match STAR SJ.tab.out 
#     ex = PD.DataFrame([x for x in _egen()], columns=GGB.BEDCOLS[:6])
#     ex['locus'] = UT.calc_locus_strand(ex)
#     ex = ex.groupby('locus').first().reset_index() # remove dup
#     sj = PD.DataFrame([x for x in _igen()], columns=GGB.BEDCOLS[:6])
#     sj['locus'] = UT.calc_locus_strand(sj)
#     sj = sj.groupby('locus').first().reset_index() # remove dup

#     UT.set_info(sj,ex)
#     UT.set_exon_category(sj, ex)

#     # find genes (connected components) set '_gidx'
#     if graphpre is None:
#         graphpre = './'+str(uuid.uuid4())+'_'
#     prefix = os.path.abspath(graphpre) # need unique prefix for parallel processing
#     genes = GP.find_genes4(sj,ex,
#         filepre=prefix,
#         np=np,
#         override=False,
#         separatese=True)

#     return sj, ex

def make_sjexci(path, np):
    if path[-3:]=='.gz':
        bpath = path[:-3]
    else:
        bpath = path
    ext = bpath[-4:]
    if ext not in ['.gtf', '.bed', '.txt']:
        raise ValueError('unknown filetype {0}, should be either .gtf,.bed (bed12),.txt (ucsc knownGene)'.format(ext))
    pathprefix = bpath[:-4]

    if not os.path.exists(path):
        raise ValueError('{0} file does not exists'.format(ext))

    if ext=='.gtf':
        df = GGB.read_gtf(path).sort_values(['chr',])
        sj, ex = gtf2exonsj(df, np=np)
    elif ext=='.bed': 
        df = GGB.read_bed(path)
        sj, ex = bed2exonsj(df, np=np)
    elif ext=='.txt': # UCSC download
        if 'knownGene' in path:
            df = GGB.read_ucsc_knownGene(path)
            sj, ex = kg2exonsj(df, np=np)
        elif 'refGene' in path:
            df = GGB.read_ucsc_refGene(path)
            sj, ex = kg2exonsj(df, np=np) # same as kg
    
    # save
    LOG.info('saving sj to {0}'.format(pathprefix+'.sj.txt.gz'))
    UT.write_pandas(sj, pathprefix+'.sj.txt.gz', 'h')
    LOG.info('saving ex to {0}'.format(pathprefix+'.ex.txt.gz'))
    UT.write_pandas(ex, pathprefix+'.ex.txt.gz', 'h')    

    # make ci
    ci = UT.chopintervals(ex, pathprefix+'.ci.txt.gz')
    return sj, ex


def make_sjex(gtfpath, dstpre, np=12):
    if UT.isstring(gtfpath):
        gtf = GGB.read_gtf(gtfpath)
    else:
        gtf = gtfpath
    sj,ex = gtf2exonsj(gtf, np=np)
    print(ex.groupby(['kind','cat']).size())
    ex.loc[ex['kind']=='5','cat'] = '5'
    ex.loc[ex['kind']=='3','cat'] = '3'
    UT.write_pandas(ex, dstpre+'.ex.txt.gz', 'h')
    UT.write_pandas(sj, dstpre+'.sj.txt.gz', 'h')
    # make ci
    ci = UT.chopintervals(ex, dstpre+'.ci.txt.gz')
    return {'sj':sj,'ex':ex}


# convienience class ###########################################################

class GTF2SJEX(object):
    """Convenience class for making sj (junction), ex (exon), ci (chopped interval) from GTF.

    """

    def __init__(self, gtfpath, outdir=None):
        self.gtfpath = gtfpath
        if outdir is None:
            self.outdir = os.path.dirname(gtfpath) # use same place as GTF
        else:
            self.outdir = outdir
        if gtfpath[-3:]=='.gz':
            assert gtfpath[-7:]=='.gtf.gz'
            base = os.path.basename(gtfpath[:-7])
        else:
            assert gtfpath[-4:]=='.gtf'
            base = os.path.basename(gtfpath[:-4])
        self.prefix = os.path.join(self.outdir, base)

    def exists(self):
        return os.path.exists(self.gtfpath)

    def fname(self,suffix):
        return '{0}.{1}'.format(self.prefix,suffix)

    def sjexpaths(self):
        sjpath = self.fname('sj.txt.gz')
        expath = self.fname('ex.txt.gz')
        return sjpath, expath

    def cipath(self):
        return self.fname('ci.txt.gz')

    def sjex(self):
        sjpath, expath = self.sjexpaths()
        # check cached
        if os.path.exists(sjpath) and os.path.exists(expath):
            LOG.info('reading sj({0}),ex({1}) from cache...'.format(sjpath,expath))
            sj = UT.read_pandas(sjpath)
            ex = UT.read_pandas(expath)
            return sj,ex
        return self.make_sjex()
        
    def make_sjex(self, np=4):
        sjpath, expath = self.sjexpaths()
        if not os.path.exists(self.gtfpath):
            raise RuntimeError('file {0} does not exist'.format(self.gtfpath))
        LOG.info('making sj,ex...')
        gtf = GGB.read_gtf(self.gtfpath) # ~ 1.5 min => 
        # if 'cov' in gtf.iloc[0]['extra']:
        #     gtf['cov'] = GGB.get_gtf_attr_col(gtf, 'cov')
        # convert gtf to sjex
        pre = self.fname('graphpre{0}_'.format(uuid.uuid4()))
        sj, ex = gtf2exonsj(gtf, np=np, graphpre=pre)
        # save
        UT.write_pandas(sj, sjpath, 'h')
        UT.write_pandas(ex, expath, 'h')
        return sj,ex

    def sj(self):
        sjpath, expath = self.sjexpaths()
        if UT.notstale(sjpath):
            return UT.read_pandas(sjpath)
        sj,ex = self.sjex()
        return sj

    def ex(self):
        sjpath, expath = self.sjexpaths()
        if UT.notstale(expath):
            return UT.read_pandas(expath)
        sj,ex = self.sjex()
        return ex

    def ci(self):
        cicols = ['chr','st','ed','name','id']
        cipath = self.cipath()
        if os.path.exists(cipath):
            LOG.info('reading ci({0}) from cache...'.format(cipath))
            ci = UT.read_pandas(cipath, names=cicols)
            return ci
        if not os.path.exists(self.gtfpath):
            raise RuntimeError('file {0} does not exist'.format(self.gtfpath))
        LOG.info('making ci..')
        sj,ex = self.sjex()
        ci = UT.chopintervals(ex, cipath)
        return ci




# SJ,EX => GTF, BED ###########################################################


def sjex2bed12(sj,ex,dstpath):
    mg = GP.MEGraph4(sj,ex,dstpath+'genegraph-')
    w = BEDWriter(mg)
    LOG.info(' writing bed12 genes to {0}...'.format(dstpath))
    w.write(dstpath)

class BEDWriter(object):

    def __init__(self, mg, se=None):
        self.mg = mg
        self.ex = ex = mg.exons.set_index('_id')
        #self.genes = genes = ex.groupby('_gidx')['_id'].groups # _gidx => [_id] dict
        #WARNING above does not give the desired dict it maps to index instead of _id
        self.genes = ex.groupby('_gidx').groups
        self.se = se
        self.i2g = dict(UT.izipcols(mg.exons, ['_id','gname']))

    def gen_subpath(self, mg, starts, maxisonum):
        e = self.ex #mg.exons
        cnt = 0
        for x in starts: # starts is a list of list
            if cnt>maxisonum:
                break
            recurse = []
            for y in mg.ex_d_ex(x[-1]): # get downstream
                n = x+[y]
                #if not e.ix[y]['name'].endswith(']'):
                if e.ix[y]['d_id']==-1:
                    # reached end
                    cnt += 1
                    yield n
                    if cnt>maxisonum:
                        break
                else:
                    recurse.append(n)
            for z in self.gen_subpath(mg, recurse, maxisonum-cnt):
                cnt += 1
                yield z
                if cnt>maxisonum:
                    break
        
    def path2rec(self, exons, name):
        exons = exons.sort_values(['st','ed'])
        chrom,strand = exons.iloc[0][['chr','strand']]
        st = exons.iloc[0]['st']
        ed = exons.iloc[-1]['ed']
        nexons = len(exons)
        esizes = ','.join(map(str, (exons['ed']-exons['st']).values))+','
        estarts = ','.join(map(str, (exons['st']-st).values))+','
        yield (chrom,st,ed,name,0,strand,st,ed,0,nexons,esizes,estarts)

    def gen_iso_gene(self, mg, ids, gid, maxisonum=10):
        # kind = gtf or bed12
        # create bed12 records of isoforms
        exons = self.ex.ix[ids] #mg.exons.ix[ids]
        # starts
        #sidx = ~exons['name'].str.startswith('[')
        sidx = exons['a_id']==-1 # not acceptor
        starts =[ [x] for x in exons[sidx]['_id']]
        paths = [x for x in self.gen_subpath(mg, starts, maxisonum)]
        for i, path in enumerate(paths):
            name = '{0}.{1}'.format(gid, i+1)
            for x in self.path2rec(self.ex.ix[path], name):
                yield x
                                    
    def gen_one_gene(self, mg, ids, gid):
        # kind = gtf or bed12
        # create bed12 records of isoforms
        for x in self.path2rec(self.ex.ix[ids], gid):
            yield x
        
    def gen_iso_all(self, maxisonum=10):
        mg, genes = self.mg, self.genes
        ex = self.ex #mg.exons
        i2g = self.i2g 
        tot = len(genes)
        cnt = 0
        for i, exx in genes.items():
            cnt += 1
            if cnt % 1000==0:
                LOG.debug('   writing genes: {0}/{1}...'.format(cnt,tot))
            gid = i2g[exx[0]]
            for x in self.gen_iso_gene(mg, exx, gid, maxisonum):
                yield x
    
    def gen_all(self):
        mg, genes = self.mg, self.genes
        ex = self.ex #mg.exons
        i2g = self.i2g 
        for i, exx in genes.items():
            gid = i2g[exx[0]]
            for x in self.gen_one_gene(mg, exx, gid):
                yield x

    def write(self, fname):
        bedex = [x for x in self.gen_all()]
        self.df = beddf = PD.DataFrame(bedex, columns=GGB.BEDCOLS)
        GGB.write_bed(beddf, fname, ncols=12)

    def write_iso(self, fname, maxisonum):
        if fname[-3:]=='.gz':
            fname = fname[:-3]
        UT.makedirs(os.path.dirname(fname))
        with open(fname,'w') as fobj:
            for x in self.gen_iso_all(maxisonum=maxisonum):
                fobj.write('\t'.join(map(str,x))+'\n')
        UT.compress(fname)

class GTFWriter(BEDWriter):

    def path2rec(self, exons, name):
        exons = exons.sort_values(['st','ed'])
        chrom,strand = exons.iloc[0][['chr','strand']]
        # name = gname.tidx
        gname = name.split('.')[0]
        for i,ex in exons.iterrows():
            extra = 'gene_id:"{0}"; transcript_id:"{1}"; etype:"{2}"; eid:"{3}"'.format(
                     gname, name, ex['ptyp'], ex['name'])
            yield (chrom,'.', 'exon',str(ex['st']+1),str(ex['ed']),0,ex['strand'],0,extra)
                    
    def write(self, fname):
        # with open(fname,'w') as fobj:
        #     for x in self.gen_all():
        #         fobj.write('\t'.join(map(str,x))+'\n')
        # UT.compress(fname)
        gtfex = [x for x in self.gen_all()]
        self.df = gtfdf = PD.DataFrame(gtfex, columns=GGB.GTFCOLS)
        GGB.write_gtf(gtfdf, fname)

#######

# gene (exons & junctions) => tree traversal  (see write_gtf)

# from ponder import graph as PG
# from ponder import utils as UT

# import miscdefs as MD
# import eg as EG

class Writer2(object):

    def __init__(self, sjexpre, th=0.1):
        self.sjexpre = sjexpre
        self.th = th
        self.sj = sj = UT.read_pandas(sjexpre+'.sj.txt.gz')
        self.ex = ex = UT.read_pandas(sjexpre+'.ex.txt.gz')
        self.mg = mg = GP.MEGraph3(sj,ex) # only consider splice junction connections
        self.exg = ex.set_index('_gidx')
        self.exi = ex.set_index('_id')
        self.nullidx = UT.find_nullidx(self.ex)
        self.e2c = dict(UT.izipcols(ex, ['_id', 'cat']))
        self.precalc_branch_p()
        
    def precalc_branch_p(self):
        # junction branch: P(sj,d_id=>a_id)
        sj = self.sj[['d_id','_id','jcnt']].copy().set_index('_id')
        sj['tot'] = sj.groupby('d_id').sum().ix[sj['d_id'].values].values
        sj['p'] = sj['jcnt']/sj['tot']
        sj.fillna(0,inplace=True) # where tot is zero ???
        # exons which share an acceptor: P(ex,a_id)
        ex = self.ex[['a_id','_id','ecov']].copy().set_index('_id')
        ex['tot'] = ex.groupby('a_id').sum().ix[ex['a_id'].values].values
        ex['p'] = ex['ecov']/ex['tot']
        ex.fillna(0,inplace=True)
        j2nd = self.mg.j2nd
        j2nd['psj'] = sj.ix[j2nd['_id'].values]['p'].values
        j2nd['pex'] = ex.ix[j2nd['e_id_a'].values]['p'].values
        j2nd['p'] = j2nd['psj']*j2nd['pex']
        self.d2ep = j2nd.groupby('e_id_d')[['e_id_a','p']]
        
    def ex_d_e2p(self,eid):
        try:
            return dict(UT.izipcols(self.d2ep.get_group(eid), ['e_id_a','p']))
        except:
            return {}
        
    def _dfs(self, exs, p, top):
        # exs: list of exon ids
        x = exs[-1]
        th = self.th
        if (self.e2c[x] == '3'): # leaf
            yield exs,p
        else:
            e2p = self.ex_d_e2p(x)
            for y in e2p: # connected exons
                pn = p*e2p[y] # new probability
                if (y!=top) and pn>th:
                    for z in self._dfs(exs+[y], pn, top):
                        yield z
                    
    def depth_first_visit(self, gid):
        # gid: _gidx
        exons = self.exg.ix[[gid]] # exons in the gene
        if len(exons)==1: # single exon
            yield [exons.iloc[0]['_id']], 1.
        else:
            # starts
            for s in exons[exons['cat']=='5']['_id']:
                for t in self._dfs([s],1.,s):
                    yield t
        
    def path2gtfrec(self, eids, gid, cnt, p):
        exons = self.exi.ix[eids].reset_index()
        exons = exons.sort_values(['st','ed'])
        # name = gname.tidx
        tid = '{0}.{1}'.format(gid,cnt)
        tmpl = 'gene_id "{0}"; transcript_id "{1}"; exon_number "{6}"; etype "{2}"; eid "{3}"; p "{4:.4f}"; ecov "{5:.4f}";'
        for i, (cat,eid,chrom,st,ed,strand,ecov) in enumerate(exons[['cat','_id','chr','st','ed','strand','ecov']].values):
            extra = tmpl.format(gid, tid, cat, eid, p, ecov, i+1)
            if '_' not in chrom:
                #yield '\t'.join([chrom,'.', 'exon',str(st),str(ed),'0',strand,'0',extra])
                # 2016-02-25 : bug! one base extra at start
                yield '\t'.join([chrom,'.', 'exon',str(st+1),str(ed),'0',strand,'0',extra])
    
    def path2bedrec(self, eids, gid, cnt, p):
        exons = self.exi.ix[eids].reset_index()
        exons = exons.sort_values(['st','ed'])
        chrom,strand = exons.iloc[0][['chr','strand']]
        st = exons.iloc[0]['st']
        ed = exons.iloc[-1]['ed']
        gcov = exons.iloc[0]['gcov']
        name = '{0}.{1}'.format(gid,cnt)
        sc1 = '{:.4f}'.format(p)
        sc2 = '{:.1f}'.format(gcov)
        nexons = str(len(exons))
        esizes = ','.join(map(str, (exons['ed']-exons['st']).values))+','
        estarts = ','.join(map(str, (exons['st']-st).values))+','
        yield '\t'.join([chrom,str(st),str(ed),name,sc1,strand,str(st),str(ed),sc2,nexons,esizes,estarts])

    def gen_gtf(self, gid):
        for i, (subpath, p) in enumerate(self.depth_first_visit(gid)):
            for rec in self.path2gtfrec(subpath, gid, i+1, p):
                yield rec
                        
    def gen_bed(self, gid):
        for i, (subpath, p) in enumerate(self.depth_first_visit(gid)):
            for rec in self.path2bedrec(subpath, gid, i+1, p):
                yield rec


    def write(self, fname, gids=None, th=None):
        if fname[-4:]=='.gtf':
            formatter = self.gen_gtf
        else:
            formatter = self.gen_bed
        if th:
            self.th = th
        if gids is None:
            gids = sorted(set(self.ex['_gidx'].values))
        with open(fname,'w') as fobj:
            for gid in gids:
                recs = [x for x in formatter(gid)]
                if len(recs)>0:
                    fobj.write('\n'.join(recs)+'\n')
        
