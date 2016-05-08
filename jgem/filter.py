"""

.. module:: filter
    :synopsis: module for filtering assembled genes

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""

import subprocess
import os
import gzip
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)
import shutil
import json

import pandas as PD
import numpy as N

from jgem import utils as UT
from jgem import gtfgffbed as GGB
from jgem import filenames as FN
from jgem import assembler as AS
from jgem import bedtools as BT
from jgem import bigwig as BW
from jgem import calccov as CC
from jgem import convert as CV
from jgem import fasta as FA
from jgem import repeats as RP


RMSKPARAMS = dict(
    np = 4,
    th_uexon=4,
    th_bp_ovl=50,
    th_ex_ovl=50,
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
            self.uex = UT.make_union_gene_bed(self.ex, '_gidx')
            UT.write_pandas(self.uex, uexpath, 'h')


    def calculate(self):
        """ Calculate base pair overlap to repeat using UCSC genome mask of repeats to lower case, 
        and exon level overlap to repeat using UCSC RepeatMaskerViz track. 
        ALso make a dataframe containing summary. 
        """
        pr = self.params
        uex = RP.count_repeats_mp(self.uex, self.gfc, np=pr['np'], col='#repbp')
        uex = RP.count_repeats_viz_mp(uex, self.rmskviz, np=pr['np'], idcol='_id', expand=0, col='repnames')
        self.ugb = self._make_gbed(self.ex, self.sj, uex )

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
        # write out filtered ex,sj,ci,unionex,gb2
        UT.write_pandas(ex2, fn.txtname('ex', category='output'), 'h')
        UT.write_pandas(sj2, fn.txtname('sj', category='output'), 'h')
        UT.write_pandas(uex2, fn.txtname('unionex', category='output'), 'h')
        UT.write_pandas(ugb2, fn.txtname('genes', category='output'), 'h')

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
        