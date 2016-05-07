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

def count_repeats(beddf, genomefastaobj, col='#repbp', returnseq=False, seqcol='seq'):
    """Looks up genome sequence and counts the number of lower characters.
    (RepeatMaker masked sequence are set to lower characters in UCSC genome)

    Args:
        beddf: Pandas DataFrame with chr,st,ed columns, when calculating repeats bp
         for genes, unioned bed should be used (use utils.make_union_gene_bed)
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
         for genes, unioned bed should be used (use utils.make_union_gene_bed)
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
                    
                    