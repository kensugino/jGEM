"""

.. module:: repeats
    :synopsis:  Repeats (transposon) related stuffs

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""

import csv
import subprocess
import os
import gzip
import glob
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)
import uuid

import pandas as PD
import numpy as N

from jgem import utils as UT
from jgem import fasta as FA
from jgem import filenames as FN
from jgem import bedtools as BT
from jgem import gtfgffbed as GGB




RMSKPARAMS = dict(
    np = 4,
    th_uexon=4,
    th_bp_ovl=50,
    th_ex_ovl=50,
    datacode='',
    gname='gname',
)

class RmskFilter(object):
    """Filter genes with by overlap to repeat masker.

    Args:
        sjexpre: path prefix to assembled ex.txt.gz, sj.txt.gz files (optionally unionex.txt.gz )
        code: identifier
        chromdir: direcotry which contains chromosomes sequences in FASTA format
        rmskviz: RepeatMasker viz track (UCSC) converted in BED7 (using jgem.repeats.rmskviz2bed7)
        outdir: output directory

    """

    def __init__(self, sjexpre, code, chromdir, rmskviz, outdir, **kw):
        self.sjexpre = sjexpre
        prefix = os.path.join(outdir, code)
        self.fnobj = FN.FileNamesBase(prefix)
        self.chromdir = chromdir
        self.rmskviz = rmskviz
        self.gfc = FA.GenomeFASTAChroms(chromdir)

        self.params = RMSKPARAMS.copy()
        self.params.update(kw)


        self.ex = UT.read_pandas(sjexpre+'.ex.txt.gz')
        self.sj = UT.read_pandas(sjexpre+'.sj.txt.gz')
        uexpath = sjexpre+'.unionex.txt.gz'        
        if os.path.exists(uexpath):
            self.uex = UT.read_pandas(uexpath)
        else:
            LOG.info('making union exons...saving to {0}'.format(uexpath))
            self.uex = UT.make_unionex(self.ex, '_gidx')
            UT.write_pandas(self.uex, uexpath, 'h')


    def calculate(self):
        """ Calculate base pair overlap to repeat using UCSC genome mask of repeats to lower case, 
        and exon level overlap to repeat using UCSC RepeatMaskerViz track. 
        ALso make a dataframe containing summary. 
        """
        pr = self.params
        uex = count_repeats_mp(self.uex, self.gfc, np=pr['np'], col='#repbp')
        uex = count_repeats_viz_mp(uex, self.rmskviz, np=pr['np'], idcol='_id', expand=0, col='repnames')
        self.ugb = self._make_gbed(self.ex, self.sj, uex, datacode=pr['datacode'], gname=pr['gname'])

    def _make_gbed(self, ex, sj, ugb, datacode='', gname='gname'):
        # rep%
        gr = ugb.groupby('_gidx')
        gb2 = gr[['chr',gname,'tlen','glen']].first()
        gb2['#repbp'] = gr['#repbp'].sum()
        gb2['rep%'] = 100.*gb2['#repbp']/gb2['tlen']
        # rmskviz, exon%
        gb2['#uexons'] = gr.size()
        gbsub = ugb[ugb['repnames']!='.(-1)'] # .(-1) == overlap
        gb2['#uexons_rmsk'] = gbsub.groupby('_gidx').size() # num exons overlapping rmskviz
        gb2.loc[gb2['#uexons_rmsk'].isnull(),'#uexons_rmsk'] = 0 
        gb2['rviz%'] = 100.*gb2['#uexons_rmsk']/gb2['#uexons']
        gb2['repnames'] = gbsub.groupby('_gidx')['repnames'].apply(lambda x: ';'.join(list(x)))
        # locus
        gb2['st'] = gr['st'].min()
        gb2['ed'] = gr['ed'].max()
        gb2['glocus'] = UT.calc_locus(gb2,'chr','st','ed')
        # rmskviz, class=[Simple_repeats, LINE, SINE, LTR, DNA]
        rcols = ['Simple_repeats','LINE','SINE','LTR','DNA']
        for k in rcols:
            gb2[k] = gb2['repnames'].str.contains('#'+k)

        dc = '_'+datacode if datacode else ''
        egr = ex.groupby('_gidx')
        gb2['#exons'] = egr.size()
        gb2['avgecov'] = egr['ecov'+dc].mean()
        gb2['gcov'] = egr['gcov'+dc].first() 

        sgr = sj.groupby('_gidx')
        gb2['ucnt'] = sgr['ucnt'+dc].sum()
        gb2['mcnt'] = sgr['mcnt'+dc].sum()
        gb2['minjcnt'] = sgr['mcnt'+dc].min()
        gb2['#junc'] = sgr.size()
        # gb2['lscore'] = N.log10(gb2['tlen']) - N.log10(gb2['glen']) + 2
        # gb2['jscore'] = N.log10(gb2['ucnt']) - N.log10(gb2['mcnt']) - 1.5    
        return gb2

    def filter(self, **kw):
        """ Filter genes.  
        base pair repeat overlap % >= th_bp_ovl (default 50)
        exon_repeat_overlap % >= th_ex_ovl (default 50)
        #union exon < th_uexon (default 4)

        That is, by default, it filters out 2,3 exon genes with both base pair and exon level
        overlap to repeats are greater or equal to 50%. Does not apply to single exons. 

        """
        d = self.ugb
        pr = self.params
        fn = self.fnobj

        idx1 = (d['rep%']>=pr['th_bp_ovl'])&(d['rviz%']>pr['th_ex_ovl'])
        idx2 = (d['#junc'].notnull())&(d['#uexons']<pr['th_uexon'])
        idx = ~(idx1&idx2)
        self.ugb2 = ugb2 = d[idx] # filtered
        gids = ugb2.index.values
        ex0 = self.ex
        sj0 = self.sj
        uex = self.uex
        # filter ex,sj,uex
        self.ex2 = ex2 = ex0[ex0['_gidx'].isin(gids)].sort_values(['chr','st','ed'])
        self.sj2 = sj2 = sj0[sj0['_gidx'].isin(gids)].sort_values(['chr','st','ed'])
        self.uex2 = uex2 = uex[uex['_gidx'].isin(gids)].sort_values(['chr','st','ed'])
        self.gbed2 = gbed2 = UT.unionex2bed12(uex2,name='gname',sc1='gcov',sc2='tlen')
        # write out filtered ex,sj,ci,unionex,gb2
        UT.write_pandas(ex2, fn.txtname('ex', category='output'), 'h')
        UT.write_pandas(sj2, fn.txtname('sj', category='output'), 'h')
        UT.write_pandas(uex2, fn.txtname('unionex', category='output'), 'h')
        UT.write_pandas(ugb2, fn.txtname('genes.stats', category='output'), 'h')
        UT.write_pandas(gbed2, fn.bedname('genes', category='output'), '') # BED12

    def save_params(self):
        UT.save_json(self.params, self.fnobj.fname('params.json', category='output'))

    def __call__(self):
        self.calculate()
        self.filter()
        self.save_params()


def plot_tlen_vs_glen_panels(gbed, fld='rep%', alpha=0.1, ms=0.8):
    fig,axr = P.subplots(2,5,figsize=(15,6), sharex=True, sharey=True)
    for i,t in enumerate(range(0,100,10)):
        ax = axr[int(i/5)][i % 5]
        _tvl(gbed, t, t+10, ax=ax, fld=fld, alpha=alpha, ms=ms)

def plot_tlen_vs_glen(gbed, title='', ax=None, mk='b.', ms=0.5, alpha=0.1):
    x = N.log10(gbed['glen'])
    y = N.log10(gbed['tlen'])
    if ax is None:
        fig, ax = P.subplots(1,1,figsize=(3,3))
    ax.plot(x.values, y.values, mk, ms=ms, alpha=alpha)
    ax.set_title('{0} (#{1})'.format(title,len(gbed)))
    ax.set_xlabel('log10(glen)')
    ax.set_ylabel('log10(tlen)')
    ax.set_xlim(1,7)
    ax.set_ylim(1.5,5.5)

def _tvl(gbed, t0, t1, alpha=0.1, ax=None, fld='rep%',ms=0.1):
    idx10 = (gbed[fld]>=t0)&(gbed[fld]<=t1)
    x = N.log10(gbed['glen'])
    y = N.log10(gbed['tlen'])
    if ax is None:
        fig, ax = P.subplots(1,1,figsize=(3,3))
    ax.plot(x[idx10].values, y[idx10].values, 'b.', ms=ms, alpha=alpha)
    ax.set_title('{0}<={3}<={1} ({2})'.format(t0,t1,N.sum(idx10),fld))
    #ax.axhline(3.2)
    ax.set_xlabel('log10(glen)')
    ax.set_ylabel('log10(tlen)')
    ax.set_xlim(1,7)
    ax.set_ylim(1.5,5.5)



def count_repeats(beddf, genomefastaobj, col='#repbp', returnseq=False, seqcol='seq'):
    """Looks up genome sequence and counts the number of lower characters.
    (RepeatMaker masked sequence are set to lower characters in UCSC genome)

    Args:
        beddf: Pandas DataFrame with chr,st,ed columns, when calculating repeats bp
         for genes, unioned bed should be used (use utils.make_unionex)
        genomefastaobj: an object with get(chr,st,ed) method that returns sequence
         (use fasta.GenomeFASTAChroms).
        col: column names where counts will be put in
        returnseq (bool): whether to return sequence or not (default False)
        seqcol: column where sequences are put in (default seq)

    Outputs:
        are put into beddf columns with colname col(default #repbp)

    """
    def _cnt(chrom,st,ed):
        seq = genomefastaobj.get(chrom,st,ed)
        return N.sum([x.islower() for x in seq])

    if returnseq:
        beddf[seqcol] = [genomefastaobj.get(*x) for x in beddf[['chr','st','ed']].values]
        beddf[col] = beddf[seqcol].apply(lambda x: N.sum([y.islower() for y in x]))
    else:
        beddf[col] = [_cnt(*x) for x in beddf[['chr','st','ed']].values]
    return beddf

def count_repeats_mp(beddf, genomefastaobj, col='#repbp', returnseq=False, seqcol='seq', idfld='_id', np=4):
    """ MultiCPU version of counts_repeats """
    # only send relevant part i.e. chr,st,ed,id
    if not idfld in beddf:
        beddf[idfld] = N.arange(len(beddf))
    # number per CPU
    n = int(N.ceil(len(beddf)/float(np))) # per CPU
    args = [(beddf.iloc[i*n:(i+1)*n],genomefastaobj,col,returnseq,seqcol) for i in range(np)]
    rslts = UT.process_mp(count_repeats, args, np=np, doreduce=False)
    df = PD.concat(rslts, ignore_index=True)
    i2c = UT.df2dict(df, idfld, col)
    beddf[col] = [i2c[x] for x in beddf[idfld]]
    if returnseq:
        i2s = UT.df2dict(df, idfld, seqcol)
        beddf[seqcol] = [i2s[x] for x in beddf[idfld]]
    return beddf


def count_repeats_viz_mp(beddf, rmskvizpath, idcol='_id', np=3, prefix=None, expand=0, col='repnames'):
    """Use rmsk-viz track and check each (unioned) exon overlaps with repeats and report repeat name(s).
    Uses Bedtools and calculates chromosome-wise.  

    Args:
        beddf: Pandas DataFrame with chr,st,ed cols, when calculating repeats bp
         for genes, unioned bed should be used (use utils.make_unionex)
        idcol: colname for unique row id (default _id)
        rmskvizpath: path to repeat masker viz BED7 file (created using rmskviz2bed7)
        np: number of CPU to use
        prefix: path prefix for temp file, if not None temp files are kept. (default None)
        expand: how many bases to expand exon region in each side (default 0)
        col: column name to put in overlapping repeat names (if multiple comma separated)

    Outputs:
        are put into beddf columns with colname col(default repnames)

    """
    cleanup = False
    if prefix is None:
        cleanup = True
        prefix = os.path.join(os.path.dirname(rmskvizpath), str(uuid.uuid4())+'_')

    # chrom-wise
    chroms = sorted(beddf['chr'].unique())
    # check whether rmskviz is already split
    splitrmsk=False
    for chrom in chroms:
        rpath = rmskvizpath+'.{0}.bed.gz'.format(chrom) # reuse
        if not os.path.exists(rpath):
            splitrmsk=True
            break
    if splitrmsk:
        rmsk = GGB.read_bed(rmskvizpath)

    args = []
    bfiles = []
    ofiles = []
    for chrom in chroms:
        bpath = prefix+'tgt.{0}.bed'.format(chrom) # don't compress
        rpath = rmskvizpath+'.{0}.bed.gz'.format(chrom) # reuse     
        if  expand>0:
            bchr = beddf[beddf['chr']==chrom].copy()
            bchr['st'] = bchr['st'] - expand
            bchr['ed'] = bchr['ed'] + expand
            bchr.loc[bchr['st']<0,'st'] = 0
        else:
            bchr = beddf[beddf['chr']==chrom]
        UT.write_pandas(bchr[['chr','st','ed',idcol]], bpath, '')
        bfiles.append(bpath)
        if splitrmsk:
            rchr = rmsk[rmsk['chr']==chrom]
            UT.write_pandas(rchr[['chr','st','ed','name','strand']], rpath, '')
        opath = prefix+'out.{0}.bed'.format(chrom)
        ofiles.append(opath)
        args.append([bpath, rpath, opath])

    rslts = UT.process_mp(count_repeats_viz_chr, args, np=np, doreduce=False)

    # gather outputs
    cols = ['name','repnames']
    outs = [UT.read_pandas(f, names=cols) for f in ofiles]
    df = PD.concat(outs, ignore_index=True)
    df['name'] = df['name'].astype(str)
    i2rn = UT.df2dict(df, 'name', 'repnames')
    beddf[col] = [i2rn[str(x)] for x in beddf[idcol]]

    # cleanup
    if cleanup:
        for f in bfiles:
            os.unlink(f)
        for f in ofiles:
            os.unlink(f)

    return beddf

def count_repeats_viz_chr(bedpath, rmskpath, outpath):
    c = BT.bedtoolintersect(bedpath, rmskpath, outpath, wao=True)
    cols = ['chr','st','ed','name','b_chr','b_st','b_ed','b_name','strand','ovl']
    df = UT.read_pandas(c, names=cols)
    df['rn'] = df['b_name']+'('+df['strand']+')'
    # group and concat repname
    dg = df.groupby('name')['rn'].apply(lambda x: ','.join(list(x))).reset_index()
    UT.write_pandas(dg, outpath, 'h')


def rmskviz2bed7(df):
    """Convert UCSC repeatmasker viz track download to BED7 format"""
    cols = ['chrom','chromStart','chromEnd','name','score','strand',
            'alignStart','alignEnd','blockSizes','blockRelStarts','id']
    # => chr[0],st(*),ed(*),name[3],sc1[4],strand[5],tst(id[-1])
    # st,ed: calculate from blocks, only use non -1 starts
    # st0[1],ed0[2],bsizes[-3],bstarts[-2]
    cols1 = ['chr','st','ed','name','sc1','strand','tst']
    def _gen():
        for x in UT.izipcols(df, cols):
            rec = [x[0],0,0,x[3],x[4],x[5],x[-1]]
            bsizes = [int(y) for y in x[-3].split(',')]
            bstarts = [int(y) for y in x[-2].split(',')]
            for y,z in zip(bstarts,bsizes):
                if y>=0:
                    rec[1] = x[1]+y
                    rec[2] = x[1]+y+z
                    yield rec.copy()
    rows = [x for x in _gen()]
    df = PD.DataFrame(rows, columns=cols1)
    return df
                    
                    