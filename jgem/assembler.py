"""

.. module:: assembler
    :synopsis: assemble genes from RNASeq data (normalized genome coverage (bigwig) and junctions)

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
from jgem import filenames as FN
from jgem import graph as GP
from jgem import calccov as CC
from jgem import convert as CV


# Assemblers ###############################################################

PARAMS = dict(
    merging=False,
    genome='mm10', # genome id
    selectsj_ratio=1e-3, # ratio for selecting overlapping junctions
    checksjsupport=False, # whether check splice junctions without read support (necessary when merging or binth>0)
    binth=0, # bigwig to bed threshold

    jie_binth=320, #16,
    jie_sjth=100,    
    jie_ratio=0.01, # junctions in a exon ratio to surrounding junctions which define the exon

    ureadth=0, # SJ uniq read (can be -1 to use uread=0 but mread>0)
    mreadth=0, # SJ non-uniq read threshold
    # SJ is selected as (uread>ureadth | mread>mreadth), so ureadth should really be >=0
    maxexonsize=35000,# max exon size (Ttn has ~20kbp exon)
    edgesize=100, # temporary edge exon size

    mpth=0.95, # mapped% th for detecting gene boundary
    cutlen=350, # check candidate bounded exons larger than this for cutting
    gap=50,# fill this gap to distinguish gene boundary vs. exon (fill exon but not intergenic)

    gap5=50, #150,#300, # for 5' UTR extension
    gap3=50, #500, #1000,# for 3' UTR extension
    covfactor=0.1, # for gap filling: if coverage of the next interval is < covfactor*current cov
    # then don't fill

    binstrand='.',
    
    iret_mpth=0.98,# mapped% th for detecting intron retension
    iret_covratio=0.01, # min cov ratio between an iret and average of surrounding exons 
    iret_covth=0.1, #0.005, # if intron cov smaller than this, then ignore


    findsecovth=True, # whether to use adaptive secovth (pndr2)
    minsecovth=5, #0.1,# minimum single exon coverage (normalized to million alignments)
    secovth=10, #0.5, # default SE cov threshold if not using adaptive version (pndr1)
    se_gap=170, #50,# single exon gap fill
    se_sizeth=50,# single exon size th
    se_sizeth2=200, # for SELECTSEME
    se_binth=0.2, #0.01,
    # adaptive threshold is calculated and bigger of this and calculated value is used
    # se_th99, FINDSE_secovth in stats is the calculated value
    findsecovth_useref=True, # whether to use ref for finding secov, if not use ME
    savepndr1=True, # whether to save pndr1 (no adaptive secovth) when doing pndr2
    
    find53ir_covratio=0.2, # cov ratio threshold for FIND53IR
    find53ir_covth=0.6, #0.03, # cov threhold for FIND53IR

    remove_overlappingse=True,# whether to remove SE overlapping to ME
    remove_bad2exon=True,# whether to perform bad 2 exon gene removal
    me2exon_sjth=2, # 2exon genes with splice junction support less than this will be removed
    me2exon_sizeth=200,# 2exon genes terminal size th ~ 2 x read length
    me2exon_covth=10, #0.5,# 2exon genes cov th
    
    ed_window=15, # edge detector smooth window (in bp)
    ed_minth=0.5, # edge detector minth for smoothed derivative
    ed_sigma=3, # edge detector sigma for threshold (larger of minth, sigma*std is used)
    ed_covratio=0.001, # edge detector if cov ratio to surrounds is < covratio, discard
    # for abundant one covratio needs to be small to retrieve 
    ed_covth=2, #0.1,# edge detector abs cov threshold
    ed_smwinsize=151, # smooth window for trimming
    ed_minintsize=10, # window for merging peaks
    ed_aggratio=0.1, # 
    ed_mimath=0.15, # min-max ratio threshold for cut decision
    ed_mimath2=0.30,
    ed_triggerth=2, # for rise detection (detect rise when >= max*mimath*triggerth after <max*mima)

    #printerr=True, # mostly for FIXEDGES part
    override=False,# override previously saved calculation
    np=1, # number of process to spawn for connected component calculation
    writegtf=False,
    writeiso=False,
    maxisonum=10,
    useallconnected=True,

    )
# use different parameters for merging
MPARAMDIFF = dict(
    merging = True,
    checksjsupport=True,
    binth=0, #0.001
    mpth=0.9999,
    gap=0,
    gap5=1,
    gap3=1,

    # findsecovth=False, # 2016-04-30 turn it on
    minsecovth=5,
    # secovth=0.05,
    # se_binth=0.01,
    se_gap=0,
    # se_sizeth=10,

    ed_covratio=0.05,
    ed_minth=0.05,
    ed_mimath=0.20,
    ed_mimath2=0.75,
    ed_sigma=5,
    # ed_covth=0.001,

    iret_mpth=1, #0.9999,
    iret_covratio=0.1,
    # iret_covth=1, #0.01,

    find53ir_covratio=0.15,
    # find53ir_covth=0.03,
)
MPARAMS = PARAMS.copy()
MPARAMS.update(MPARAMDIFF)

# for documentation purpose
PARAMSDOC = PARAMS.copy()
for k,v in MPARAMS.items():
    PARAMSDOC[k+'_m'] = v

#LOG.debug(PARAMSDOC)

class Assembler(object):

    def __init__(self, fnobj, merging=False, saveintermediates=False, **kwargs):
        """
        Args:
            fnobj: FileNames object
            merging (bool): whether merging assembled models or not (default False)
            saveintermediates (bool): whether to save intermediates (default True)
            kwargs: to change any of the parameter values

        """
        self.fnobj=fnobj
        if merging:
            pr = MPARAMS.copy()
        else:
            pr = PARAMS.copy()

        self.saveintermediates = saveintermediates
        pr.update(kwargs)
        self.params = pr
        self.stats = {}

    # def delete_intermediates(self):
    #     "Delete intermediate files."
    #     fn = self.fnobj
    #     categories = [x for x in fn._fnames.keys() if x !='output']
    #     fn.delete(delete=categories, protect=['output'])

    def check_params(self):
        "Check parameter change (if run previously) and save current parameter set."
        fn = self.fnobj
        pr = self.params.copy()
        del pr['override']
        del pr['np']
        # parameter check
        if (pr['binth']>0) or (pr['merging']): 
            # always do consistency check if binth>0 (since support may be reduced)
            pr['checksjsupport'] = self.params['checksjsupport']=True
        prdf = PD.DataFrame(pr,index=['value']).T
        fname = fn.fname('assemble.params.txt', category='output')
        if os.path.exists(fname):
            prdf0 = UT.read_pandas(fname,index_col=[0])
            if len(prdf)==len(prdf0):
                if not all(prdf['value'].astype(str) == prdf0['value'].astype(str)):
                    self.params['override']=True
                    LOG.warning('parameter different overriding...')
                    for x in prdf.index:
                        p1 = str(prdf.ix[x]['value'])
                        p2 = str(prdf0.ix[x]['value'])
                        if p1 != p2:
                            LOG.info('  {0}: {1} <=> {2}'.format(x,p1,p2))
            else:
                self.params['override']=True
                LOG.warning('parameter set changed, overriding...')
                p1 = set(prdf.index.values)
                p2 = set(prdf0.index.values)
                a = p1.difference(p2)
                b = p2.difference(p1)
                LOG.info('  old({0}) new({1})'.format(a,b))
        # save parameters
        UT.write_pandas(prdf, fname, 'ih')
            
    def save_stats(self):
        """ Saves assemble related stats into (samplename).assemble.stats.txt. """
        df = PD.DataFrame(self.stats, index=['value']).T
        fname = self.fnobj.fname('assemble.stats.txt',category='stats')
        UT.write_pandas(df, fname, 'ih')

    def assemble(self):
        """ Perform the assembly """
        fn = self.fnobj
        pr = self.params
        st = self.stats
        self.check_params()

        self.sj = GGB.read_bed(fn.sjfile)

        if not pr['merging']:
            SELECTSJ(self)()
        if pr['checksjsupport']:
            CHECKSJSUPPORT(self)()

        REMOVEJIE(self)()
        SJ2EX(self)()
        # ADDJIE(self)() # no need to do this
        MERGEEXONS(self)()
        if pr['merging']:
            FINDEDGES2(self)()
            FIXSTRAND(self)()
            # FIND53IR only deals with exons with length > SE sizeth = 50bp)
            # gap smaller than this will be missed if we skip FINDIRETS
            FINDIRETS(self)()
            FINDSECOVTH(self)()
            FIND53IR(self)()
        else:
            FINDEDGES(self)()
            FIXSTRAND(self)()
            EDGEFIXER(self)()
            FINDIRETS(self)()
            FINDSECOVTH(self)()
            FINDSE(self)()

        CALCCOV(self)()
        SETINFO(self)() # SET ACCEPTOR/DONOR/EXON CATEGORY
        FINDGENES(self)()
        SELECTSEME(self)() # SELECT ME and SE, saves exname2, exname3, sjname2
        if not pr['merging']:
            FIXEDGES2(self)() # TRIM 3',5' edges

        CONSISTENTSJ(self)() # remove sj without ex support
        WRITESJEX(self)()
        WRITEGENES(self)()

        self.save_stats()
        if not self.saveintermediates:
            # self.delete_intermediates()
            fn.delete(delete=['temp'], protect=['output','stats'])


# assembler modules #######################################################

# def mp_worker(args):
#     func, arg = args
#     return func(*arg)

class SUBASE(object):
    """Base class of assembler modules."""

    def __init__(self, asm):
        self.asm = asm # assembler
        self.fnobj = asm.fnobj # FileNames object
        self.params = asm.params # dict 
        self.stats = asm.stats # dict
        self.info = ''

    def __call__(self, *args, **kwargs):
        LOG.info('{0} '.format(self.__class__.__name__)+'='*20)
        _sttime = time.time()
        rslt = self.call(*args, **kwargs)
        if self.info:
            LOG.info(' '+self.info)
        LOG.info(' time: {0:.3f}s'.format(time.time()-_sttime))

    def call(self, *args, **kwargs):
        raise NotImplementedError

    def bw2bed(self, binth):
        fn = self.fnobj
        pr = self.params
        binfile = BW.bw2bed(
            bwfile=fn.bwfile,
            bedfile=fn.bedname2('bw',binth),
            chroms=UT.chroms(pr['genome']),
            th=binth
        )
        return binfile

    def fillgap(self, gap, binfile):
        fn = self.fnobj
        if gap==0:
            gapfile=binfile
        else:
            gapfile = BT.fillgap(binfile,fn.bedname2('gap', gap), gap=gap)
        return gapfile

    def chroms(self, df):
        pr = self.params
        chroms0 = set(df['chr'].unique())
        return [x for x in UT.chroms(pr['genome']) if x in chroms0]

    # def _process_mp(self, func, args):
    #     np = self.params['np']
    #     rslts = []
    #     if np==1:
    #         for i, arg in enumerate(args):
    #             rslts += func(*arg)
    #             LOG.debug(' processing: {0}/{1}...'.format(i+1,len(args)))
    #     else:
    #         try:
    #             p = multiprocessing.Pool(np)
    #             a = zip(repeat(func), args)
    #             tmp = p.map(mp_worker, a)
    #         finally:
    #             LOG.debug('closing pool')
    #             p.close()
    #         rslts = reduce(iadd, tmp)
    #     return rslts

    def sjfile(self):
        fn = self.fnobj
        pr = self.params
        fname0 = fn.bedname('checksjsupport.sj')
        fname1 = fn.bedname('fixstrand.sj')
        if os.path.exists(fname1):
            return fname1
        if pr['checksjsupport'] and os.path.exists(fname0):
            return fname0
        return fn.sjfile


class SELECTSJ(SUBASE):
    """Select splice junction based on ratio of junction counts within overlapping junctions.
    Use ratio of unique counts and ratio of total counts (to incorporate cases where there are
    only non-unique counts). 

    Args:
        sj: junction DataFrame

    Returns:
        :sj: selected junctions dataframe

    Related Parameters:
        * selectsj_ratio: thredhold for the ratio, default:{selectsj_ratio}

    Files:
        * selectsj.bed.gz
        * selectsj.inte.txt.gz

    """

    def call(self):
        sj = self.asm.sj
        fn = self.fnobj
        pr = self.params

        a = b = GGB.write_bed(sj, fn.bedname('selectsj'), ncols=7)
        c = fn.txtname('selectsj.inte')
        c = BT.bedtoolintersect(a,b,c,wao=True) # somehow -s (force same strand) doesn't work
        o = BT.read_ovl(c, GGB.BEDCOLS[:7])
        idxstrand = o['strand']==o['b_strand']
        LOG.debug(a)
        LOG.debug(c)
        LOG.debug(o.head())
        sjgr = o.groupby(['chr','st','ed','strand'])
        # BEDCOLUMNS sc1(5th), tst(7th) contains unique count (ucnt) and non-unique count (mcnt)
        sj2 = sjgr[['name','sc1','tst']].first()
        sj2['ucnt_sum'] = sjgr['b_sc1'].sum()
        sj2['mcnt_sum'] = sjgr['b_tst'].sum()
        sj2['sum'] = sj2['ucnt_sum']+sj2['mcnt_sum']
        sj2['cnt'] = sj2['sc1']+sj2['tst']
        self.sj2 = sj2 = sj2.reset_index()
        sj2['locus'] = UT.calc_locus_strand(sj2)
        sj2['ratio'] = sj2['sc1']/sj2['ucnt_sum']
        sj2['ratio_m'] = sj2['tst']/sj2['mcnt_sum']
        sj2['ratio_a'] = sj2['cnt']/sj2['sum']        
        self.sj2 = sj2
        UT.write_pandas(sj2, fn.txtname('selectsj.sj2',category='output'), 'h') # TODO change cat to 'temp'

        # select 
        th_ratio = pr['selectsj_ratio']
        idx1 = (sj2['ratio']>=th_ratio)|(sj2['ratio_a']>=th_ratio)
        self.sj4 = sj4 = sj2[idx1]
        self.info = '#sj:{0}=>{1}'.format(len(sj), len(sj4))
        self.stats['SELECTSJ.#sj0'] = len(sj)
        self.stats['SELECTSJ.#sj'] = len(sj4)
        #return sj4
        self.asm.sj = sj4

class CHECKSJSUPPORT(SUBASE):
    """Check junction edges have >0 coverages. Remove junctions without support. 

    Args:
        sj: junction DataFrame

    Returns:
        :sj: compatible junctions

    Related Parameters:
        * binth: coverage threshold, default:{binth}, (for merging: {binth_m})
        * genome: genome version, default:{genome}

    TempFiles:
        * bw*.bed.gz
        * sjst.bed.gz
        * sjed.bed.gz
        * sjst.ovl.txt.gz
        * sjed.ovl.txt.gz

    """

    def call(self):
        sj = self.asm.sj
        fn = self.fnobj
        pr = self.params
        # get bindf
        binfile = self.bw2bed(pr['binth'])
        # write st,ed from sjs
        sj['st-1']=sj['st']-1
        sj['st-2']=sj['st']-2
        sj['ed+1']=sj['ed']+1
        sj['_id']=N.arange(len(sj))
        sjst = fn.bedname('checksjsupport.sjst')
        sjed = fn.bedname('checksjsupport.sjed')
        # BEDCOLS: chr,st,ed,name,sc1,strand,tst
        # UT.write_pandas(sj[['chr','st-2','st-1','_id','sc1','strand','tst']],sjst,'')
        # UT.write_pandas(sj[['chr','ed','ed+1','_id','sc1','strand','tst']],sjed,'')
        UT.write_pandas(sj[['chr','st-2','st-1','_id']],sjst,'')
        UT.write_pandas(sj[['chr','ed','ed+1','_id']],sjed,'')
        stovl = fn.txtname('checksjsupport.sjst.ovl')
        edovl = fn.txtname('checksjsupport.sjed.ovl')
        # bedtools intersect to bindf
        # ost = BT.calc_ovlratio(sjst, binfile, stovl, nacol=7, nbcol=3, idcol='name')
        # oed = BT.calc_ovlratio(sjed, binfile, edovl, nacol=7, nbcol=3, idcol='name')
        ost = BT.calc_ovlratio(sjst, binfile, stovl, nacol=4, nbcol=3, idcol='name')
        oed = BT.calc_ovlratio(sjed, binfile, edovl, nacol=4, nbcol=3, idcol='name')
        ost = ost.set_index('name')
        oed = oed.set_index('name')
        sjsupp = sj.set_index('_id')[(ost['ovlratio']>0)&(oed['ovlratio']>0)].copy()
        self.info ='#sj: {0}=>{1}'.format(len(sj), len(sjsupp))
        self.stats['CHECKSJSUPPORT.#sj'] = len(sj)
        self.stats['CHECKSJSUPPORT.#sjsupp'] = len(sjsupp)
        #return sjsupp
        fn.write_bed(sjsupp, 'checksjsupport.sj', ncols=7)
        self.asm.sj = sjsupp

class REMOVEJIE(SUBASE):
    """Remove Junctions In Exons. Often times exons with high coverage contain 
    noise junctions. 

    Args:
        sj: junction DataFrame

    Returns:
        :sj: junction dataframe without JIE
        :jie: junctions in exons

    Related Parameters:
        * jie_binth: coverage threshold,  default:{jie_binth}
        * jie_sjth: threshold for normalized read counts,  default:{genome}

    TempFiles:
        * bw*.bed.gz
        * sjsupp.bed.gz
        * jie.bw.ovl.txt.gz

    """

    def call(self):
        sj = self.asm.sj
        sjfile = self.sjfile()
        fn = self.fnobj
        pr = self.params

        # covarage file
        binfile = self.bw2bed(pr['jie_binth'])
        # if nothing in binfile then skip
        try:
            jiebw = GGB.read_bed(binfile)
        except:
            self.asm.jie = None
            self.info = 'nothing above jie_binth {0}'.format(pr['jie_binth'])
            self.stats['REMOVEJIE.#sj'] = len(sj)
            self.stats['REMOVEJIE.#jie'] = 0
            return            

        if len(jiebw)==0:
            self.asm.jie = None
            self.info = 'nothing above jie_binth {0}'.format(pr['jie_binth'])
            self.stats['REMOVEJIE.#sj'] = len(sj)
            self.stats['REMOVEJIE.#jie'] = 0
            return            

        sjmp = BT.calc_ovlratio(
            aname=sjfile, 
            bname=binfile, 
            tname=fn.txtname('removejie.bw.ovl'), 
            nacol=7, nbcol=3, 
            idcol=['chr','st','ed','strand']
        )
        # match records between sjmp and mg.sj
        sjmp['str_id'] = UT.calc_locus(sjmp)
        sj['str_id'] = UT.calc_locus(sj)
        sid2ovl = UT.df2dict(sjmp, 'str_id','ovlratio')
        sj['ovlratio'] = [sid2ovl.get(x,N.nan) for x in sj['str_id']]
        
        # should use count ratios instead of actual reads as threshold ?
        th = pr['jie_sjth']
        idx = (sj['ovlratio']==1)&(sj['sc1']<th)&(sj['tst']<th)
        

        sj1 = sj[~idx].copy() # use these for "nearest donor/acceptor" exon extraction
        jie = sj[idx].copy() # junctions in exon, add later
        self.info = '#sj:{0}=>{1}, jie {2}'.format(len(sj), len(sj1), len(jie))
        self.stats['REMOVEJIE.#sj'] = len(sj1)
        self.stats['REMOVEJIE.#jie'] = len(jie)
        #return sj1, jie
        self.asm.sj = sj1
        self.asm.jie = jie

class SJ2EX(SUBASE):
    """Find candidate exons from junctions.

    Args:
        sj: junction DataFrame

    Returns:
        :me: exons dataframe
        :sj: compatible junctions

    Related Parameters:
        * ureadth: threshold for junction unique counts, default:{ureadth} ({ureadth_m} for merging)
        * mreadth: threshold for non-unique junction counts, default:{mreadth} ({mreadth_m} for merging)
        * maxexonsize: maximum exon size, default:{maxexonsize}
        * edgesize: temporary edge exon size, default:{edgesize}


    """

    def call(self):
        sj = self.asm.sj
        fn = self.fnobj
        pr = self.params
        ureadth=pr['ureadth']
        mreadth=pr['mreadth']

        sj['st-1'] = sj['st']-1
        sj['ureads'] = sj['sc1']
        sj['mreads'] = sj['tst']
        sj1 = sj[(sj['ureads']>ureadth)|(sj['mreads']>mreadth)].copy()
        LOG.info('#sj:{0}=>{1} after ureadth, mreadth'.format(len(sj), len(sj1)))
        self.stats['SJ2EX.#sj_before_uth_mth'] = len(sj)
        self.stats['SJ2EX.#sj_after_uth_mth'] = len(sj1)

        ex = PD.DataFrame([x for x in self._sj2ex(sj1)],columns=GGB.BEDCOLS[:6])
        # there are small cases of duplicates
        ex = ex.sort_values(['chr','st','ed','strand'])
        ex['_ord'] = -ex['name'].str.len() # bring bounded both [...] to front (others lack)
        exg = ex.groupby(['chr','st','ed','_ord'])
        ex = exg.first().reset_index()
        self.info = '#ex:{0}, #sj:{1}'.format(len(ex), len(sj1))
        self.stats['SJ2EX.#exon_candidates'] = len(ex)
        #return ex, sj1
        self.asm.me = ex
        self.asm.sj = sj1

    def _sj2ex(self, sj):
        pr = self.params
        maxsize=pr['maxexonsize'] 
        if pr['merging']:
            edgesize = maxsize
        else:
            edgesize=pr['edgesize']

        counter = [0]
        # workaround for nonlocal closure variable for Python2
        # in Python3 use nonlocal keyword
                
        def sj2expn(strand):
            # strand +,-
            df = sj[(sj['strand']==strand)]
            # chrom wise
            for chrom, g in df.groupby('chr'):
                sts = sorted(set(g['st-1'].values)) # right end of exon (ed)
                eds = sorted(set(g['ed'].values)) # left end of exon (st)
                eds = [max(0,sts[0]-edgesize)]+list(eds)
                nsts = len(sts)
                neds = len(eds)
                j = 0 # current left end pos
                usededs = set()
                for i in range(nsts): # go through right ends
                    while (eds[j+1]<sts[i]): # find nearest left end
                        j += 1
                    if j==0:
                        usededs.add(j)
                        counter[0] +=1
                        name = '{0}{1}s]'.format(strand,counter[0])
                        yield (chrom, eds[j], sts[i], name, 0, strand)
                    elif (sts[i]-eds[j])>maxsize: # too far
                        counter[0] +=1
                        name = '{0}{1}f]'.format(strand,counter[0])
                        yield (chrom, sts[i]-edgesize, sts[i], name, 0, strand)
                    else: # bounded on both
                        usededs.add(j)
                        counter[0] +=1
                        name = '[{0}{1}a]'.format(strand,counter[0])
                        yield (chrom, eds[j], sts[i], name, 0, strand)
                unusededs = sorted(set(range(neds)).difference(usededs))
                i = 0
                sts = list(sts)+[eds[-1]+edgesize]
                for j in unusededs:
                    if j==0: # dummy ed
                        continue 
                    while (sts[i+1]<=eds[j]):
                        i += 1
                    if i+1==nsts:# len(sts)=nsts+1
                        counter[0] +=1
                        name = '[{0}{1}e'.format(strand,counter[0])
                        yield (chrom, eds[j], sts[i+1], name, 0, strand)
                    elif (sts[i+1]-eds[j])>maxsize: # too far
                        counter[0] +=1
                        name = '[{0}{1}f'.format(strand,counter[0])
                        yield (chrom, eds[j], eds[j]+edgesize, name, 0, strand)
                    else: # bounded on both sides
                        counter[0] +=1
                        name = '[{0}{1}b]'.format(strand,counter[0])
                        yield (chrom, eds[j], sts[i+1], name, 0, strand)

        def sj2exns():
            sj0 = sj
            df = sj[(sj['strand']=='.')]
            # chrom wise
            for chrom, g in df.groupby('chr'):
                sts = sorted(set(g['st-1'].values)) # right end of exon (ed) for '.'
                sjchr = sj0[sj0['chr']==chrom]#.copy()
                #sjchr['st'] = sjchr['st']-1
                tmp = [tuple(x) for x in sjchr[['st-1','strand']].values.tolist()]
                sts0 = sorted(set(tmp)) # all right end of exon (ed)+strand
                tmp = [tuple(x) for x in sjchr[['ed','strand']].values.tolist()]
                eds0 = sorted(set(tmp)) # all left end of exon (st)+strand
                sts0 = sts0+[(eds0[-1][0]+edgesize,'.')]
                eds0 = [(max(0,sts0[0][0]-edgesize),'.')]+eds0
                nsts = len(sts)
                j = 0 # current left end pos
                usededs = set()
                for i in range(nsts): # go through right ends
                    while (eds0[j+1][0]<sts[i]): # find nearest left end
                        j += 1
                    if j==0:
                        if eds0[j][1]=='.':
                            usededs.add(j)
                        counter[0] +=1
                        name = '.{0}s]'.format(counter[0])
                        yield (chrom, eds0[j][0], sts[i], name, 0, eds0[j][1])
                    elif (sts[i]-eds0[j][0])>maxsize: # too far
                        counter[0] +=1
                        name = '.{0}f]'.format(counter[0])
                        yield (chrom, sts[i]-edgesize, sts[i], name, 0,  eds0[j][1])
                    else: # bounded on both
                        if eds0[j][1]=='.':
                            usededs.add(j)
                        counter[0] +=1
                        name = '[.{0}a]'.format(counter[0])
                        yield (chrom, eds0[j][0], sts[i], name, 0, eds0[j][1])
                alleds = set([i for i,x in enumerate(eds0) if (x[1]=='.')&(i!=0)])
                # eds0[0] is a dummy record => don't include
                unusededs = sorted(alleds.difference(usededs))
                i = 0
                nsts0 = len(sts0)
                for j in unusededs:
                    while (sts0[i+1][0]<=eds0[j][0]):
                        i += 1
                    if i==nsts0:
                        counter[0] +=1
                        name = '[.{0}e'.format(counter[0])
                        yield (chrom, eds0[j][0], sts0[i+1][0], name, 0, sts0[i+1][1])
                    elif (sts0[i+1][0]-eds0[j][0])>maxsize: # too far
                        counter[0] +=1
                        name = '[.{0}f'.format(counter[0])
                        yield (chrom, eds0[j][0], eds0[j][0]+edgesize, name, 0, eds0[j][1])
                    else: # bounded on both sides
                        counter[0] +=1
                        name = '[.{0}b]'.format(counter[0])
                        yield (chrom, eds0[j][0], sts0[i+1][0], name, 0, sts0[i+1][1])

        for x in sj2expn('+'):
            yield x
        for x in sj2expn('-'):
            yield x
        for x in sj2exns():
            yield x
        LOG.debug(' total {0} exon candidates'.format(counter[0]))

class MERGEEXONS(SUBASE):
    """Merge overlapping exons.

    Args:
        me: exon DataFrame

    Returns:
        :me: exon dataframe with merged exons

    TempFiles:
        * sjex.bed.gz
        * sjex.inte.txt.gz

    """
    # TODO:
    #     Currently only two overlapping exons are merged. Generalized to n.


    def call(self):
        ex = self.asm.me
        fn = self.fnobj
        pr = self.params
        # ex vs. ex overlap
        a = b = GGB.write_bed(ex, fn.bedname('mergeexons.sjex'), ncols=6)
        c = fn.txtname('mergeexons.sjex.inte')
        c = BT.bedtoolintersect(a,b,c,wao=True) # somehow -s (force same strand) doesn't work
        cols0 = GGB.BEDCOLS[:6]
        o = BT.read_ovl(c, cols0)
        idxstrand = o['strand']==o['b_strand']
        # select ordered overlaps (to count overlap only once)
        idx1 = o['st']<o['b_st'] 
        idx2 = (o['st']==o['b_st'])&(o['ed']<o['b_ed'])
        o = o[idxstrand&(idx1|idx2)] # unique overlap pairs
        def _gen():
            cols = ['chr','st','ed','strand','name','b_st','b_ed','b_name']
            for c,s,e,t,n,bs,be,bn in UT.izipcols(o,cols):
                if e!=be:
                    yield (c,s,be,n+'+'+bn,0,t)
                if s!=bs:
                    yield (c,bs,e,bn+'+'+n,0,t)
        mrgd = PD.DataFrame([x for x in _gen()], columns=cols0)
        me = PD.concat([ex[cols0], mrgd[cols0]], ignore_index=True)
        me = me.groupby(['chr','st','ed','strand']).first().reset_index()
        self.info = '#ex:{0}'.format(len(me))
        #return me
        self.asm.me = me

class FINDEDGES2(SUBASE):
    """Find edges 

    Args:
        sj: junction DataFrame
        me: exon DataFrame

    Returns:
        :sj: compatible junctions
        :me: exon dataframe with edge exons added

    Related Parameters:
        * mpth: mapped% th for detecting gene boundary, default {mpth}
        * binth: bigwig to bed threshold, default {binth} (merging: {binth_m})
        * gap: fill this gap to distinguish gene boundary vs. exon (fill exon but not intergenic), default {gap} (merging: {gap_m})
        * edgesize: temporary edge exon size, default {edgesize}
        * maxexonsize: max exon size (Ttn has ~20kbp exon), default {maxexonsize}
        * cutlen: check candidate bounded exons larger than this for cutting, default {cutlen}
        * ed_sigma: edge detector sigma for threshold (larger of minth, sigma*std is used), 
          default {ed_sigma}  (merging: {ed_sigma_m})
        * ed_minth: edge detector minth for smoothed derivative, default {ed_minth} (merging: {ed_minth_m})
        * ed_covratio: edge detector if cov ratio to surrounds is < covratio, discard, 
          default {ed_covratio} (merging: {ed_covratio_m})
        * ed_window: edge detector smooth window (in bp), default {ed_window}
        * ed_covth: edge detector abs cov threshold, default {ed_covth}  (merging: {ed_covth_m})
        * ed_smwinsize: smooth window for trimming default {ed_smwinsize}
        * ed_minintsize: window for merging peaks, default {ed_minintsize}
        * ed_aggratio: default {ed_aggratio} 
        * ed_mimath: min-max ratio threshold for cut decision, default {ed_mimath} (merging: {ed_mimath_m})
        * ed_mimath2: default {ed_mimath2}  (merging: {ed_mimath2_m})
        * ed_triggerth: for rise detection (detect rise when >= max*mimath*triggerth after <max*mima), default {ed_triggerth}

    TempFiles:
        * fe2.sjbb.bed.gz
        * fe2.sjbb-ovl.txt.gz
        * fe2.me1.ci.txt.gz
        * fe2.me2p.ci.txt.gz
        * fe2.me2n.ci.txt.gz
        * fe2.sj.txt.gz
        * fe2.me.txt.gz

    """

    def call(self):
        sj = self.asm.sj
        me = self.asm.me
        fn = self.fnobj
        pr = self.params
        st = self.stats
        self.me = me
        self.sj = sj

        # bounded (kind: a,b)=> cut
        me['_id2'] = N.arange(len(me))
        idx = me['name'].str.contains('s|e|f') # unbounded
        me1 = me[~idx] # bounded cut these
        me2 = me[idx] # not (a,b) kind:(s,e,f) => fix (~20K)

        # do it at the level of ci 
        ci = UT.chopintervals(me1, fn.txtname('findedges2.me1.ci'), idcol='_id2')
        binfile = self.bw2bed(pr['binth'])
        sjbbname = UT.write_pandas(ci[['chr','st','ed','name','id']], fn.bedname('findedges2.sjbb'), '')
        bbg = BT.calc_ovlratio(
                aname=sjbbname, 
                bname=binfile, 
                tname=fn.txtname('findedges2.sjbb-ovl'), 
                nacol=5, nbcol=3
        )
        # calculate mp and len
        bbg['len'] = bbg['ed'] - bbg['st']
        bbg['name1'] = bbg['name'].astype(str).apply(lambda x:[int(y) for y in x.split(',')])
        # name1 contains list of ids
        ci1 = bbg[(bbg['ovlratio']<pr['mpth'])|(bbg['len']>=pr['cutlen'])] # ~ 16K
        # cut candidates
        eids1 = reduce(iadd, ci1['name1'].values, []) # exons being cut
        eids0 = sorted(set(me1['_id2'].values).difference(set(eids1)))
        me1i = me1.set_index('_id2') # me1 (bounded) indexed
        me1s = me1i.ix[eids1].reset_index() # cut targets
        me0 = me1i.ix[eids0].reset_index() # remaining (non-cut) exons
        LOG.debug('cutting {0} ci chrom-wise'.format(len(ci1)))
        me1r = self.process_mp(ci1,me1s,cutedges_m)

        # fix me2
        me2p = me2[me2['name'].str.startswith('[')]
        ci2p = UT.chopintervals(me2p, fn.txtname('findedges2.me2p.ci'), idcol='_id2')
        ci2p = ci2p.rename(columns={'id':'sc1'})
        ci2p['direction'] = '+'
        me2n = me2[me2['name'].str.endswith(']')]
        ci2n = UT.chopintervals(me2n, fn.txtname('findedges2.me2n.ci'), idcol='_id2')
        ci2n = ci2n.rename(columns={'id':'sc1'})
        ci2n['direction'] = '-'
        LOG.debug('fixing {0} cip chrom-wise'.format(len(ci2p)))
        me2pr = self.process_mp(ci2p,me2p,fixedges_m)
        LOG.debug('fixing {0} cin chrom-wise'.format(len(ci2n)))
        me2nr = self.process_mp(ci2n,me2n,fixedges_m)
        # concatenate
        cols = ['chr','st','ed','name','sc1','strand','_id2']
        me3 = PD.concat([me0[cols],me1r[cols],me2pr[cols],me2nr[cols]],ignore_index=True).sort_values(['chr','st','ed'])
        me3 = me3.groupby(['chr','st','ed','strand']).first().reset_index()

        # find consistent sj
        UT.set_info(sj,me3)
        # acceptor
        aids = set(me3['a_id'].values).intersection(set(sj['a_id'].values))
        dids = set(me3['d_id'].values).intersection(set(sj['d_id'].values))
        sja = sj.set_index('a_id').ix[aids].reset_index()
        dids = dids.intersection(set(sja['d_id'].values))
        sjd = sja.set_index('d_id').ix[dids].reset_index()
        #return sjd, me3
        fn.write_txt(sjd, 'findedges2.sj')
        fn.write_txt(me3, 'findedges2.me')
        self.asm.sj = sjd
        self.asm.me = me3

    def process_mp(self, ci, me, func):
        pr = self.params
        bwname = self.fnobj.bwfile
        args = []
        for chrom in self.chroms(me):
            mechr = me[me['chr']==chrom][['chr','st','ed','name','strand','_id2']].copy()
            cichr = ci[ci['chr']==chrom].copy()
            args.append((cichr,mechr,bwname,pr,chrom))
        rslts = UT.process_mp(func, args, pr['np'])
        # rslts = []
        # if pr['np']==1:
        #     for i,arg in enumerate(args):
        #         tmp = func(*arg) #cutedges_m(arg)
        #         LOG.debug('  processing {3}: {0}/{1} {2}...'.format(i+1,len(args),len(tmp),arg[-1]))
        #         rslts.append(tmp)
        #     return PD.concat(rslts,ignore_index=True)
        # try:
        #     p = multiprocessing.Pool(pr['np'])
        #     tmp = p.map(func, args)
        # finally:
        #     LOG.debug('  closing pool')
        #     p.close()
        # return PD.concat(tmp, ignore_index=True)
        return PD.DataFrame(rslts, columns=['chr','st','ed','name','sc1','strand','_id2'])

def cutedges_m(ci,me,bwname,pr,chrom):
    # input ci: chopped interval, me: exons, pr: params
    # return cut exons
    # convert ci => cicut  [+side col]
    edgedet = EdgeDetector(
                    bwname=bwname,
                    sigmath=pr['ed_sigma'],
                    minth=pr['ed_minth'],
                    covratio=pr['ed_covratio'],
                    winsize=pr['ed_window'],
                    covth=pr['ed_covth'],
                    gapth=pr['binth'],
                    gap=pr['gap'],
                    smwinsize=pr['ed_smwinsize'],
                    minintsize=pr['ed_minintsize'],
                    aggregateratio=pr['ed_aggratio'],
                    mimath=pr['ed_mimath'],
                    mimath2=pr['ed_mimath2'],
                    triggerth=pr['ed_triggerth'])
    cicols = ['chr','st','ed','name','sc1']
    def _gen():
        with edgedet:
            for chrom0,st0,ed0,name,cid in UT.izipcols(ci,cicols):
                cuts = edgedet.cutboth(chrom0,st0,ed0)
                # left
                for chrom1,st1,ed1 in cuts[0]:
                    yield (chrom,st1,ed1,name,cid, 0) # 0: left
                # right 
                for chrom1,st1,ed1 in cuts[1]:
                    yield (chrom,st1,ed1,name,cid, 1) # 1: right
                # both
                for chrom1,st1,ed1 in cuts[2]:
                    yield (chrom,st1,ed1,name,cid, 2) # 2: both
    tmp = [c for c in _gen()]
    #LOG.debug('cut ci {2}: {0}=>{1}'.format(len(ci),len(tmp),chrom))
    df = PD.DataFrame(tmp, columns=['chr','st','ed','name','cid','side'])
    # now put them together as exons
    # eid (_id2 in me) <=> cid (in df): encoded in df[name]
    df['eid'] = df['name'].astype(str).apply(lambda x:[int(y) for y in x.split(',')])
    dff = UT.flattendf(df, 'eid')
    # eid=>strand, eid=>name
    e2strand = dict(UT.izipcols(me, ['_id2','strand']))
    e2name = dict(UT.izipcols(me, ['_id2','name']))
    e2st = dict(UT.izipcols(me, ['_id2','st'])) 
    e2ed = dict(UT.izipcols(me, ['_id2','ed']))
    dff['strand'] = [e2strand[x] for x in dff['eid']]
    dff['ename'] = [e2name[x] for x in dff['eid']]
    dff['est'] = [e2st[x] for x in dff['eid']]
    dff['eed'] = [e2ed[x] for x in dff['eid']]
    dff = dff.sort_values(['eid','cid'])
    
    def _egen():
        for eid, gr in dff.groupby('eid',sort=False):
            if len(gr)==1 and gr.iloc[0]['side']==2:# no cut
                c,s,e,n,r,cid = gr.iloc[0][['chr','est','eed','ename','strand','cid']]
                yield (c,s,e,n,cid,r,eid)
                continue
            cids = sorted(gr['cid'].unique())
            gri = gr.set_index('cid')
            cnt = 1
            # left 
            for cid in cids:
                gci = gri.ix[[cid]].set_index('side')
                # side 0
                try:
                    test = gci.ix[0]
                    s0 = gci.ix[[0]]
                    for c,s,e,n,r in s0[['chr','est','ed','ename','strand']].values:
                        n1 = '{0}*(l{1})'.format(n,cnt)
                        cnt += 1
                        yield (c,s,e,n1,cid,r,eid)
                except:
                    pass
                # if not side 2 break
                try:
                    s2 = gci.ix[2]
                    if cid==cids[-1]: # last one all connected
                        c,s,e,n,r = s2[['chr','est','eed','ename','strand']].values
                        yield (c,s,e,n,-1,r,eid)
                except:
                    break
            # right
            for cid in cids[::-1]:
                gci = gri.ix[[cid]].set_index('side')
                # side 1
                try:
                    test = gci.ix[1]
                    s1 = gci.ix[[1]]
                    for c,s,e,n,r in s1[['chr','st','eed','ename','strand']].values:
                        n1 = '(r{1})*{0}'.format(n,cnt)
                        cnt += 1
                        yield (c,s,e,n1,cid,r,eid)
                except:
                    pass
                # if not side 2 break
                try:
                    s2 = gci.ix[2]
                except:
                    break

    tmp2 = [e for e in _egen()]
    #LOG.debug('cut ex {2}: {0}=>{1}'.format(len(me),len(tmp2),chrom))
    return tmp2
    # edf = PD.DataFrame(tmp2, columns=['chr','st','ed','name','sc1','strand','_id2'])
    # return edf
    
def fixedges_m(ci,me,bwname,pr,chrom):
    # input me: exons only need to be fixed in one direction (kind: s,e,f)
    # startswith [ => '+', endswith ] => '-'
    edgedet = EdgeDetector(bwname,
                        sigmath=pr['ed_sigma'],
                        minth=pr['ed_minth'],
                        covratio=pr['ed_covratio'],
                        winsize=pr['ed_window'],
                        covth=pr['ed_covth'],
                        gapth=pr['binth'],
                        gap=pr['gap'],
                        smwinsize=pr['ed_smwinsize'],
                        minintsize=pr['ed_minintsize'],
                        aggregateratio=pr['ed_aggratio'],
                        mimath=pr['ed_mimath'],
                        mimath2=pr['ed_mimath2'],                        
                        triggerth=pr['ed_triggerth'])
    cicols = ['chr','st','ed','name','sc1','direction']
    def _gen():
        with edgedet:
            for chrom0,st0,ed0,name,cid,d in UT.izipcols(ci,cicols):
                cuts = edgedet.cutone(chrom0,st0,ed0,d)
                # left
                for chrom1,st1,ed1 in cuts[0]:
                    yield (chrom,st1,ed1,name,cid, 0, d) # 0: left
                # right 
                for chrom1,st1,ed1 in cuts[1]:
                    yield (chrom,st1,ed1,name,cid, 1, d) # 1: right
                # both
                for chrom1,st1,ed1 in cuts[2]:
                    yield (chrom,st1,ed1,name,cid, 2, d) # 2: both
    tmp = [c for c in _gen()]
    #LOG.debug('cut ci2 {2}: {0}=>{1}'.format(len(ci),len(tmp),chrom))
    df = PD.DataFrame(tmp, columns=['chr','st','ed','name','cid','side','direction'])
    # now put them together as exons
    # eid (_id2 in me) <=> cid (in df): encoded in df[name]
    df['eid'] = df['name'].astype(str).apply(lambda x:[int(y) for y in x.split(',')])
    dff = UT.flattendf(df, 'eid')
    # eid=>strand, eid=>name
    e2strand = dict(UT.izipcols(me, ['_id2','strand']))
    e2name = dict(UT.izipcols(me, ['_id2','name']))
    e2st = dict(UT.izipcols(me, ['_id2','st'])) 
    e2ed = dict(UT.izipcols(me, ['_id2','ed']))
    dff['strand'] = [e2strand[x] for x in dff['eid']]
    dff['ename'] = [e2name[x] for x in dff['eid']]
    dff['est'] = [e2st[x] for x in dff['eid']]
    dff['eed'] = [e2ed[x] for x in dff['eid']]
    dff = dff.sort_values(['eid','cid'])
    def _egen():
        for eid, gr in dff.groupby('eid',sort=False):
            if len(gr)==1 and gr.iloc[0]['side']==2:# no cut
                c,s,e,n,r,cid = gr.iloc[0][['chr','est','eed','ename','strand','cid']]
                yield (c,s,e,n,cid,r,eid)
                continue
            cids = sorted(gr['cid'].unique())
            gri = gr.set_index('cid')
            cnt = 1
            direction = gr.iloc[0]['direction']
            # only do left or right depending on direction
            if direction=='+':
                # left 
                for cid in cids:
                    gci = gri.ix[[cid]].set_index('side')
                    # side 0
                    try:
                        test = gci.ix[0]
                        s0 = gci.ix[[0]]
                        for c,s,e,n,r in s0[['chr','est','ed','ename','strand']].values:
                            n1 = '{0}*(l{1})'.format(n,cnt)
                            cnt += 1
                            yield (c,s,e,n1,cid,r,eid)
                    except:
                        pass
                    # if not side 2 break
                    try:
                        s2 = gci.ix[2]
                        if cid==cids[-1]: # last one all connected
                            c,s,e,n,r = s2[['chr','est','eed','ename','strand']].values
                            yield (c,s,e,n,-1,r,eid)
                    except:
                        break
            else:
                # right
                for cid in cids[::-1]:
                    gci = gri.ix[[cid]].set_index('side')
                    # side 1
                    try:
                        test = gci.ix[1]
                        s1 = gci.ix[[1]]
                        for c,s,e,n,r in s1[['chr','st','eed','ename','strand']].values:
                            n1 = '(r{1})*{0}'.format(n,cnt)
                            cnt += 1
                            yield (c,s,e,n1,cid,r,eid)
                    except:
                        pass
                    # if not side 2 break
                    try:
                        s2 = gci.ix[2]
                        if cid==cids[0]: # last one all connected
                            c,s,e,n,r = s2[['chr','est','eed','ename','strand']].values
                            yield (c,s,e,n,-1,r,eid)                        
                    except:
                        break
    tmp2 = [e for e in _egen()]
    #LOG.debug('cut ex {2}: {0}=>{1}'.format(len(me),len(tmp2),chrom))
    return tmp2
    # edf = PD.DataFrame(tmp2, columns=['chr','st','ed','name','sc1','strand','_id2'])
    # return edf

class FINDEDGES(SUBASE):
    """Find edge exons. (Add * to the name when cut.)

    Args:
        me: exon DataFrame

    Returns:
        :me: exon dataframe with edge exons added

    Related Parameters:
        * mpth: mapped% th for detecting gene boundary, default:{mpth} (merging: {mpth_m})
        * binth: default {binth} (merging: {binth_m})
        * gap: fill this gap to distinguish gene boundary vs. exon (fill exon but not intergenic)
          default {gap} (merging: {gap_m})
        * edgesize: temporary edge exon size, default {edgesize}

    TempFiles:
        * fe.sjbb.bed.gz
        * fe.sjbb-ovl.txt.gz
        * fe.exons.bed.gz

    """

    def call(self):
        sjexdf = self.asm.me
        fn = self.fnobj
        pr = self.params
        gap = pr['gap']
        edgesize = pr['edgesize']

        # covarage file
        binfile = self.bw2bed(pr['binth'])
        gapfile = self.fillgap(gap, binfile)

        # write sjex bounded both
        idx = (sjexdf['name'].str.startswith('['))&(sjexdf['name'].str.endswith(']'))
        sjbb = sjexdf[idx]
        sjbbname = fn.write_bed(sjbb, 'findedges.sjbb', ncols=6)
        
        # bedtools intersect
        ovlname = fn.txtname('findedges.sjbb-ovl')
        bbg = BT.calc_ovlratio(sjbbname, gapfile, ovlname, nacol=6, nbcol=3)
        
        # write internal exons (>pr['mpth'])
        inexs = bbg[bbg['ovlratio']>=pr['mpth']]

        # turn ones with ratio<pr['mpth'] into edges and write together with 
        # other edge exons
        edexs1 = sjexdf[~idx]
        edges = bbg[bbg['ovlratio']<pr['mpth']]
        cols = GGB.BEDCOLS[:6] # ['chr','st','ed','name','sc1','strand']
        cols0 = ['chr','st','ed','strand','ovlratio','name']
        if len(edges)>0:
            def _iter_edges():
                for chrom,st,ed,strand,ratio,name in edges[cols0].values:
                    # name = [(strand)(counter)(kind)]
                    yield (chrom,st,st+edgesize,name[:-1]+'*',ratio,strand) 
                    # name = [(strand)(counter)(kind)*
                    yield (chrom,ed-edgesize,ed,'*'+name[1:],ratio,strand) 
                    # name = *(strand)(counter)(kind)]
            edexs2 = PD.DataFrame([x for x in _iter_edges()], columns = cols)
            edexs = PD.concat([edexs1[cols],edexs2[cols]], ignore_index=True)
        else:
            edexs = edexs1[cols]
        edexs = edexs.sort_values(['chr','st','ed','strand','name'])
        edexs = edexs.groupby(['chr','st','ed','strand'],sort=False).first().reset_index()
        
        # return concat
        exons = PD.concat([inexs[cols], edexs[cols]], ignore_index=True)
        #return exons
        fn.write_bed(exons, 'findedges.exons', ncols=6)
        self.asm.me = exons

class FIXSTRAND(SUBASE):
    """Assign strand to unstranded elements using information from connected components.

    Args:
        sj: junction dataframe
        me: exon dataframe

    Returns:
        :sj: junction dataframe with strand fixed
        :me: exon dataframe with strand fixed

    Related Parameters:
        * useallconnected (bool): whether to use all connected components or just direct neighbors
          default {useallconnected}

    TempFiles:
        * sj0.bed.gz

    """

    def call(self):
        sj = self.asm.sj
        me = self.asm.me
        fn = self.fnobj
        pr = self.params
        st = self.stats
        
        useallconnected=pr['useallconnected']
        # ureadth = pr['ureadth']
        # mreadth = pr['mreadth']
        # sj = sj[((sj['ureads']>ureadth)|(sj['mreads']>mreadth))].copy()
        
        mg = GP.MEGraph2(sj,me) # connections
        # fix unstranded exons
        idx = me['strand']=='.'
        tmp = me[idx]
        if not useallconnected:
            cnt = [Counter(me.ix[mg.ex_ex(x)]['strand'].values) for i,x in enumerate(tmp['_id'].values)]
        else:
            cnt = [Counter(me.ix[mg.connected(x)]['strand'].values) for i,x in enumerate(tmp['_id'].values)]
        LOG.debug('#unstranded={0}'.format(N.sum(idx)))
        st['FIXSTRAND.#unstranded_exons_before'] = N.sum(idx)
        me.loc[idx,'strand'] = [x.most_common()[0][0] for x in cnt]
        idx = me['strand']=='.'
        LOG.debug('#unstranded_exons_after={0}'.format(N.sum(idx)))
        st['FIXSTRAND.#unstranded_exons_after'] = N.sum(idx)
        # make sure there's no missing strand col
        me.loc[me['strand'].isnull(), 'strand'] = '.'
        # fix unstranded junctions
        idx = sj['strand']=='.'
        tmp = sj[idx]
        cnt = [Counter(mg.sj_ex(x,'strand')) for i,x in tmp.iterrows()]
        LOG.debug('#unstranded junctions={0}'.format(N.sum(idx)))
        st['FIXSTRAND.#unstranded_junctions_before'] = N.sum(idx)
        sj.loc[idx,'strand'] = [x.most_common()[0][0] if len(x)>0 else '.' for x in cnt]
        idx = sj['strand']=='.'
        LOG.debug('#still unstranded={0}'.format(N.sum(idx)))
        st['FIXSTRAND.#unstranded_junctions_after'] = N.sum(idx)
        if N.sum(idx)>0:
            LOG.warning('discarding {0} still unstranded junctions...'.format(N.sum(idx)))
            sj = sj[~idx].copy()
        #return sj, me
        fn.write_bed(sj, 'fixstrand.sj', ncols=7)
        self.asm.sj = sj
        self.asm.me = me

def fixedge(posarr,exs,strand,gap,utr,ignorefirstdonors,covfactor):

    def _gen():
        NAME = {True:{'gap_too_large':'==',
                      'assertfail_2_donor':'=D',
                      'assertfail_2_bined':'=B',
                      'acceptor':'_a',
                      'donor':'=d]',
                      'opposite':'_o'},
                False:{'gap_too_large':'==',
                      'assertfail_2_donor':'D=',
                      'assertfail_2_bined':'B=',
                       'acceptor':'a_',
                       'donor':'[d=',
                       'opposite':'o_'}}
        flag = ((strand=='+')&(utr=='3pr')) or ((strand=='-')&(utr=='5pr'))
        delta = 1 if flag else -1
        def _nposkind(i):
            npos, kind, cov = posarr.ix[i][['pos','kind','cov']]
            return npos, kind[1:], cov # remove initial one char which was for sorting

        def _step1(i,dcp):# look for bined or donor
            # (None(close) or newpos, reason(kind), newidx, coverage)
            if (i+1==len(posarr)):
                return None,'outofrange_1',i+1,0.
            npos, kind, cov = _nposkind(i+1)
            while((kind=='acceptor')|(kind=='binst')|((kind=='donor')&(delta*(dcp-npos)>0))): 
                i += 1
                if (i+1==len(posarr)):
                    return None,'outofrange_1',i+1, 0.
                npos, kind, cov = _nposkind(i+1)
            #if kind=='binst': # this cannot be the case<== actually if using different gap then happens
            #    return None,'assertfail_1_binst',i+1, cov
            if kind=='opposite': # stop here
                return npos, kind, i+1, cov
            # the rests are bined and donor, prioritize donor
            while(i+1<len(posarr)):
                npos2, kind2, cov2 = _nposkind(i+1)
                if npos2 != npos:
                    break
                i += 1
                #if kind2=='donor':
                #    kind = 'donor'
            return npos, kind, i, cov # i+1 is next pos so return i

        def _step2(i):# look for next interval
            if (i+1==len(posarr)):
                return None,'outofrange_2',i+1,0.
            npos, kind, cov = _nposkind(i+1)
            if kind in ['bined','donor']: # this cannot be the case
                return None,'assertfail_2_'+kind,i+1,0
            # the rests are binst, acceptor, or opposite
            return npos, kind, i+1, cov

        def _next(cpos,i0,donorcheckpos,cov1=None):
            # cpos: current end of interval could be None
            # i0: current index pos into posarr
            # donorcheckpos: (for SE attachment) ignore donor within this 
            # cov1: current coverage, could be None

            # cpos only updated at bined (end of interval) or at donor 

            npos1, kind1, i1, cov1n = _step1(i0,donorcheckpos) # step1 find ends
            if cov1 is None:
                cov1 = cov1n
            if npos1 is None: # didn't find end=> close at current end
                # (close?, cur_end_pos, reason, cur_idx, coverage)
                return (True, cpos, kind1, i1, cov1) 
            if kind1=='bined': # found end of an interval
                # but it could be an interval that we ignored, so check coverage
                if cov1*covfactor <= cov1n: # update pos/ use '<= ' so that cov1=cov1n=0 is ok
                    cpos = npos1
                # then find the start of the next interval
                npos2, kind2, i2, cov2 = _step2(i1)
                if npos2 is None:
                    if kind2=='outofrange_2': # no further intervals
                        return (True, cpos,'gap_too_large',i1, None) # close at current pos
                    return (True, cpos, kind2, i2, None) # assert fail
                # kind2 must be binst or acceptor or opposite
                if kind2=='binst':
                    if abs(npos2-cpos)<=gap: # keep moving
                        return (False, cpos, 'binst', i2, cov1)
                    return (True, cpos, 'gap_too_large', i1, None) # close
                # kind2 must be acceptor or opposite, close
                return (True, cpos, kind2, i1, None) 
            # must be donor or opposite, close
            if kind1=='donor':
                if cov1*covfactor <= cov1n: # update pos/ use '<= ' so that cov1=cov1n=0 is ok
                    cpos = npos1                
                return (True, cpos, kind1, i1, None) 
            # must be opposite
            return (True, cpos, kind1, i1, None) 


        # ex: ['chr','st','ed','name',...]
        if flag:
            tgt,tgt2 = 1,2 # st, ed
        else:
            tgt,tgt2 = 2,1 # ed, st
        for ex in exs.values:
            ex = list(ex)
            pos = ex[tgt]
            if ignorefirstdonors:
                donorcheckpos = ex[tgt2]
            else:
                donorcheckpos = ex[tgt]
            tmp = posarr[posarr['pos']==pos]['idx']
            if len(tmp)==0:
                continue
            i0 = tmp.iloc[-1]
            #LOG.debug('======= initial pos={0}'.format(pos))
            close,npos,kind,i0,cov = _next(pos,i0,donorcheckpos,None)
            while(not close):
                close,npos,kind,i0,cov = _next(npos,i0,donorcheckpos,cov)
            if (npos is not None) and (npos!=pos): # yield if there's end
                if flag:
                    ex[3] = ex[3]+NAME[flag][kind]
                else:
                    ex[3] = NAME[flag][kind]+ex[3]
                ex[tgt2] = npos
                if kind.startswith('assert'):
                    LOG.debug('  error:{0}/{1}/{2}/{3}/{4}'.format(UT.exid(ex), kind, npos, i0,ex[3]))
                    LOG.debug(posarr[i0-3:i0+3])
                yield ex
            else:
                #if printerr:
                LOG.debug('  error:{0}/{1}/{2}/{3}/{4}'.format(UT.exid(ex), kind, npos, i0,ex[3]))
                if not kind.startswith('outofrange'):
                    LOG.debug(posarr[i0-3:i0+3])
    
    rslt = [x for x in _gen()]    
    return rslt
        
class EDGEFIXER(SUBASE):
    """Fix edge exons by extending.

    Args:
        me: exon dataframe

    Returns:
        :me: exon dataframe with edge exons fixed
        :edgefixer: this instance for reuse later

    Related Parameters:
        * gap3: for 3' UTR extension, default {gap3} (merging: {gap3_m})
        * gap5: for 5' UTR extension, default {gap5} (merging: {gap5_m})
        * covfactor: for gap filling: if coverage of the next interval is < covfactor*current cov
          default {covfactor}

    TempFiles:
        * fixed5pr.bed.gz
        * fixed3pr.bed.gz
        * edgefixer.me.bed.gz

    """

    # [TODO]
    # - how to distinguish gaps in real 3'UTRs vs. gaps between genes or within exons?
    #   => use coverage information?
    #   => EDGEFIXER2 does this
    
    def call(self):
        me = self.asm.me[GGB.BEDCOLS[:6]] # part of the code assumes standard BED order
        fn = self.fnobj
        pr = self.params
        st = self.stats
        
        gap3 = pr['gap3']
        gap5 = pr['gap5']
        override = pr['override']
        covfactor = pr['covfactor']

        # make interval bed with strand info
        self.bindf = bindf = self.make_bindf(me) #sname, exons, binth=binth, override=override)
        
        # 3prime
        fname = fn.bedname('edgefixer.fixed3pr')
        if override or (not os.path.exists(fname)):
            LOG.debug(' fixing 3pr exons...')
            fixed3pr = self.process_edges(me, bindf, utr='3pr',gap=gap3,covfactor=covfactor)
            GGB.write_bed(fixed3pr, fname, ncols=6)
        else:
            LOG.debug(' reading cached fixed3pr ...')
            fixed3pr = GGB.read_bed(fname)
        # find reconnected and remove the other side
        
        def _remove(fixed, me):
            tmp0 = fixed[fixed['name'].str.endswith('*=d]')]
            tmp1 = fixed[fixed['name'].str.startswith('[d=*')]
            remove0 = dict([('*'+x[1:-4]+']',p) for x,p in UT.izipcols(tmp0,['name','ed'])])
            remove1 = dict([('['+x[4:-1]+'*',p) for x,p in UT.izipcols(tmp1,['name','st'])])
            LOG.debug('  removing {0}/{1}: me len={2}'.format(len(remove0),len(remove1),len(me)))
            me0a = me[me['name'].isin(remove0.keys())]
            me1a = me[me['name'].isin(remove1.keys())]
            LOG.debug('  me0a({0}), me1a({1})'.format(len(me0a),len(me1a)))
            me0b = me0a[me0a['ed']==[remove0[x] for x in me0a['name']]]
            me1b = me1a[me1a['st']==[remove1[x] for x in me1a['name']]]
            LOG.debug('  me0b({0}), me1b({1})'.format(len(me0b),len(me1b)))
            idx0 = list(me0b.index.values)+list(me1b.index.values)
            idx1 = me.index.isin(idx0)
            LOG.debug('  idx0({0}), idx1({1}))'.format(len(idx0),N.sum(idx1)))
            me1 = me[~idx1]
            LOG.debug('  after removing: me len={0} diff={1}'.format(len(me1), len(me)-len(me1)))
            return me1
            
        me = _remove(fixed3pr, me)
        
        # 5prime
        fname = fn.bedname('edgefixer.fixed5pr')
        if override or (not os.path.exists(fname)): 
            LOG.debug(' fixing 5pr exons...')
            fixed5pr = self.process_edges(me, bindf,utr='5pr',gap=gap5,covfactor=covfactor)
            GGB.write_bed(fixed5pr, fname, ncols=6)
        else:
            LOG.debug(' reading cached fixed5pr ...')
            fixed5pr = GGB.read_bed(fname)
        me = _remove(fixed5pr, me)
        
        UT.set_ptyp(me)
        inexs = me[me['ptyp']=='i']
        UT.set_ptyp(fixed3pr)
        UT.set_ptyp(fixed5pr)
        

        cols = GGB.BEDCOLS[:6]+['ptyp']
        nexons = PD.concat([inexs[cols], fixed3pr[cols], fixed5pr[cols]], ignore_index=True)
        self.nexons = nexons = nexons.groupby(['chr','st','ed','strand','ptyp']).last().reset_index()
        #return nexons
        fn.write_bed(nexons, 'edgefixer.me', ncols=6)
        self.asm.me = nexons
        self.asm.edgefixer = self

    def fixSEedge(self, me, targets, utr):
        pr = self.params
        if utr=='3pr':
            gap = pr['gap3']
        else:
            gap = pr['gap5']
        return self.process_edges_subset(me, self.bindf, targets, utr, gap=gap,covfactor=pr['covfactor'])

    def make_bindf(self, me):
        #sname, exons, binth=0, override=False
        fn = self.fnobj
        pr = self.params
        override = pr['override']

        # assign strand using exons
        bfile = fn.txtname('edgefixer.bindf')
        if (not os.path.exists(bfile)) or override:
            binfile = self.bw2bed(pr['binth'])
            if pr['binstrand']=='.':
                efile = fn.write_bed(me, 'edgefixer.exons', ncols=6)
                ofile = fn.txtname('edgefixer.bindf-exons')
                ofile = BT.bedtoolintersect(binfile,efile,ofile,wao=True)
                cols = ['chr','st','ed','ex_chr','ex_st','ex_ed','name','sc1','strand','ovl']
                df = UT.read_pandas(ofile,names=cols)
                tg = df.groupby(['chr','st','ed'])
                t2 = tg.first() # choose first, then fix intervals with multiple elements
                t2s = tg.size()
                tidx = t2s[t2s>1].index # one with more than 1 overlapping element
                tmp = df.set_index(['chr','st','ed']).ix[tidx]
                cnt = tmp.groupby(tmp.index)['strand'].apply(lambda x: len(set(x)))
                sidx = cnt[cnt>1].index # more than one strand
                if len(sidx)>0:
                    t2.ix[sidx,'strand'] = '.' # multiple strand over this inverval => set to '.'
                t3 = t2.reset_index()[['chr','st','ed','name','ovl','strand']]
                t4 = CC.calc_cov_mp(t3, fn.bwfile, bfile, pr['np'])
                # clean up
                os.unlink(efile)
                os.unlink(ofile)
                return t4 
            else:
                df = GGB.read_bed(binfile) # chr,st,ed
                df['strand'] = pr['binstrand']
                UT.save_tsv_nidx_whead(df, bfile)
                return df
        else:
            return UT.read_pandas(bfile)
        # subtract internal exons ...=> wasn't a good idea
        #     inexfile = os.path.join(MD.SJEXDIR, sname+'.sj.inex.bed.gz')
        #     binfile2 = os.path.join(MD.SJEXDIR, sname+'.bindf.bed')
        #     if (not os.path.exists(binfile2+'.gz')) or override:
        #         BT.bedtoolsubtract(binfile, inexfile, binfile2)
        #     return GGB.read_bed(binfile2+'.gz')

    def _make_arr(self, chrom, strand, bindf, exons, utr='3pr'):
        binchr = bindf[(bindf['chr']==chrom)&(bindf['strand'].isin([strand,'.']))]
        exchr = exons[(exons['chr']==chrom)&(exons['strand']==strand)]
        #exchr2 = exons[(exons['chr']==chrom)&(exons['strand']!=strand)] # opposite strand
        lex = exchr[exchr['name'].str.startswith('[')]
        rex = exchr[exchr['name'].str.endswith(']')]
        def _todf(df,tgt,name):
            tmp = PD.DataFrame(df[tgt].values, columns=['pos'])
            tmp['kind'] = name
            if 'cov' in df.columns:
                tmp['cov'] = df['cov'].values
            else:
                tmp['cov'] = 0.
            return tmp
        if ((strand=='+')&(utr=='3pr')) or ((strand=='-')&(utr=='5pr')) :
            posarr = PD.concat([_todf(binchr,'st','2binst'), # going right direction
                                _todf(binchr,'ed','2bined'), 
                                # make sure acceptor/donor comes before binst/ed (initial 1&2)
                                _todf(lex,'st','1acceptor'), 
                                # 5pr donor regard as acceptor for the purpose of edge extension
                                _todf(rex,'ed','1donor'),
                                #_todf(exchr2,'st','opposite')
                               ],
                               ignore_index=True).sort_values(['pos','kind'])
            eidx = (exchr['name'].str.startswith('['))&(~exchr['name'].str.endswith(']'))
        else:# (-,3') or (+,5')
            posarr = PD.concat([_todf(binchr,'st','1bined'), # going left direction
                                _todf(binchr,'ed','1binst'),
                                # since going opposite direction acceptor/donor after binst/ed
                                _todf(lex,'st','2donor'),
                                _todf(rex,'ed','2acceptor'),
                                #_todf(exchr2,'ed','opposite')
                               ],
                               ignore_index=True).sort_values(['pos','kind'], ascending=False)
            eidx = (~exchr['name'].str.startswith('['))&(exchr['name'].str.endswith(']'))
        exs = exchr[eidx]
        posarr.index = N.arange(len(posarr))
        posarr['idx'] = posarr.index.values #N.arange(len(posarr))
        return posarr, exs


    def process_edges(self, exons, bindf, utr, gap=300, covfactor=0.2):
        pr = self.params
        LOG.debug('  preparing data...')
        args = []
        for chrom in self.chroms(exons): # exons['chr'].unique():
            for strand in ['+','-','.']:
                posarr, exs = self._make_arr(chrom,strand,bindf,exons,utr)
                if len(exs)==0:
                    continue
                args.append((posarr, exs, strand, gap, utr, False, covfactor))#,printerr))
        rslts = UT.process_mp(fixedge, args, pr['np'])
        # rslts = []
        # if np==1:
        #     for arg in args:
        #         rslts += fixedge(arg)
        # else:
        #     try:
        #         p = multiprocessing.Pool(np)
        #         tmp = p.map(fixedge, args)
        #         LOG.debug('done fixEdge calculation: np={0}'.format(np))
        #     finally:
        #         LOG.debug('closing pool')
        #         p.close()
        #     rslts = reduce(iadd, tmp)
        return PD.DataFrame(rslts, columns=exons.columns)

    def process_edges_subset(self, exons, bindf, targets, utr, gap=300, covfactor=0.2):
        LOG.debug('  preparing data...')
        args = []
        for chrom in self.chroms(exons): #exons['chr'].unique():
            for strand in ['+','-','.']:
                idx = (targets['chr']==chrom)&(targets['strand']==strand)
                if N.sum(idx)==0:
                    continue
                posarr, exs = self._make_arr(chrom,strand,bindf,exons,utr)
                exs = targets[idx]
                args.append((posarr, exs, strand, gap, utr, True, covfactor))#,printerr))
        rslts = UT.process_mp(fixedge, args, self.params['np'])
        # rslts = []
        # if np==1:
        #     for arg in args:
        #         rslts += fixedge(arg)
        # else:
        #     try:
        #         p = multiprocessing.Pool(np)
        #         tmp = p.map(fixedge, args)
        #         LOG.debug('done fixEdge calculation: np={0}'.format(np))
        #     finally:
        #         LOG.debug('closing pool')
        #         p.close()
        #     rslts = reduce(iadd, tmp)            
        return PD.DataFrame(rslts, columns=targets.columns)

class FINDIRETS(SUBASE):
    """Find intron retentions.

    Args:
        sj: junction DataFrame
        me: exon DataFrame

    Returns:
        :me: exons with irets

    Related Parameters:
        * iret_mpth: mapped% th for detecting intron retension,
          default {iret_mpth} (merging: {iret_mpth_m})
        * iret_covth: if intron cov smaller than this, then ignore
          default {iret_covth} (merging: {iret_covth_m})
        * iret_covratio: min cov ratio between an iret and average of surrounding exons 
          default {iret_covratio} (merging: {iret_covratio_m})
        * binth: default {binth} (merging: {binth_m})


    TempFiles:
        * irets.bed.gz
        * irets.exons.txt.gz
        * irets.me.cov.txt.gz
        * irets.me.bed.gz
        * irets.me.covci.txt.gz
        * irets.me.ci.txt.gz
        * sj.iret.sub.bed.gz
        * sj.iret.ci.txt.gz
        * sj.iret.covci.txt.gz

    """

    def call(self):
        self.sj = sj = self.asm.sj
        self.me = me = self.asm.me
        fn = self.fnobj
        pr = self.params

        override = pr['override']
        iretmpth = pr['iret_mpth']
        binth = pr['binth']
        
        dfname = fn.bedname('findirets.irets')
        if override or (not os.path.exists(dfname)):
            LOG.debug(' finding iret...')
            irets = self.find_irets(sj, me, mpth=iretmpth, binth=binth)
        else:
            LOG.debug(' reading cached iret...')
            irets = GGB.read_bed(dfname)

        irets['ptyp'] = 'r'
        cols = GGB.BEDCOLS[:6]+['ptyp']
        if 'ptyp' not in me.columns:
            UT.set_ptyp(me)
        nexons = PD.concat([me[cols],irets[cols]], ignore_index=True)
        self.nexons = nexons = nexons.groupby(['chr','st','ed','strand','ptyp']).first().reset_index()
        #return nexons
        fn.write_txt(nexons, 'findirets.exons', fm='h')
        self.asm.me = nexons

    def find_irets(self, sj, me, mpth=0.95, binth=0):
        fn = self.fnobj
        pr = self.params
        st = self.stats
        override = pr['override']
        # covarage file
        binfile = self.bw2bed(binth)
        sjfile = self.sjfile()
        cname = fn.txtname('findirets.sj.mp')
        sjmp = BT.calc_ovlratio(
            aname=sjfile, 
            bname=binfile, 
            tname=cname, 
            nacol=7, 
            nbcol=3
        )
        # match records between sjmp and mg.sj
        sjmp['str_id'] = UT.calc_locus(sjmp)
        sj['str_id'] = UT.calc_locus(sj)
        sid2ovl = UT.df2dict(sjmp, 'str_id','ovlratio')
        sj['ovlratio'] = [sid2ovl.get(x,N.nan) for x in sj['str_id']]
        if '_id' not in sj.columns:
            UT.set_ids(sj)
        if '_id' not in me.columns:
            UT.set_ids(me)
        if ('st_id' not in sj.columns) or ('st_id' not in me.columns):
            UT.set_pos_info(sj,me) 

        sj['st-1'] = sj['st']-1
        sj2 = sj[sj['ovlratio']>=mpth].copy()

        # calc sj cov 1. subtract me, 2. calc_cov_ovl
        mefile = fn.write_bed(me, 'findirets.me', ncols=3)
        sjname = fn.bedname('findirets.sj')
        sjname = UT.write_pandas(sj2[['chr','st-1','ed','_id']], sjname, '')
        cname = fn.bedname('findirets.sj.sub')
        if override or (not os.path.exists(cname)):
            cname = BT.bedtoolsubtract(sjname, mefile, cname)
        sub = GGB.read_bed(cname) 
        # bed: idcol=name, moreover 
        # sj is reduced (some sj completely overlap with exon)
        sj2 = sj2.set_index('_id').ix[sub['name'].values].reset_index()

        ciname = fn.txtname('findirets.sj.ci')
        if override or (not os.path.exists(ciname)):
            ci = UT.chopintervals(sub, ciname, idcol='name') 
        else:
            ci = UT.read_pandas(ciname, names=['chr','st','ed','name','id'])

        covciname = fn.txtname('findirets.sj.covci')
        if override or (not os.path.exists(covciname)):
            covci = CC.calc_cov_mp(ci, fn.bwfile, covciname, np=pr['np'])
        else:
            covci = UT.read_pandas(covciname)

        covci['len'] = covci['ed']-covci['st']
        covci['val'] = covci['cov']*covci['len']            
        covci['sjid'] = covci['name'].apply(lambda x: [int(y) for y in x.split(',')])
        cov = UT.flattendf(covci[['id','sjid','len','val']], 'sjid')
        covg = cov.groupby('sjid')[['len','val']].sum().reset_index()
        covg['cov'] = covg['val']/covg['len']
        sj2cov = UT.df2dict(covg, 'sjid','cov')
        sj2['cov'] = [sj2cov[x] for x in sj2['_id']]
        sj2 = sj2[sj2['cov']>pr['iret_covth']]

        # calc me cov
        if override or not os.path.exists(fn.txtname('findirets.me.cov')):
            LOG.debug('calculating ME cov...')
            self.me = me = CC.calc_cov_ovl_mp(
                srcname=me, 
                bwname=fn.bwfile, 
                dstname=fn.txtname('findirets.me.cov'), 
                override=override, 
                np=pr['np'],
                covciname=fn.txtname('findirets.me.covci'),
                ciname=fn.txtname('findirets.me.ci'),
            )
        else:
            self.me = me = fn.read_txt('findirets.me.cov')


        self.irets = irets = sj2
        LOG.info('#irets candidates:{0}'.format(len(irets)))
        st['FINDIRETS.#irets_sj'] = len(irets)

        return self.process_mp(sj2, me, irets)

    def process_mp(self, sj, me, irets):
        fn = self.fnobj
        pr = self.params
        covratio = pr['iret_covratio']
        covth = pr['iret_covth']
        LOG.debug(' preparing data...')
        args = []
        for chrom in self.chroms(me):
            mechr = me[me['chr']==chrom][['chr','st','ed','name','_id','st_id','ed_id','cov']].copy()
            sjchr = sj[sj['chr']==chrom][['chr','st','ed','name','_id','st_id','ed_id']].copy()
            irchr = irets[irets['chr']==chrom][['chr','ovlratio','strand','st_id','ed_id','cov']]
            args.append((sjchr,mechr,irchr,chrom,covratio,covth))
        rslts = UT.process_mp(findirets, args, pr['np'])
        # rslts = []
        # np = pr['np']
        # if np==1:
        #     for i,arg in enumerate(args):
        #         tmp = findirets(arg)
        #         LOG.debug('  processing {3}: {0}/{1} {2}...'.format(i+1,len(args),len(tmp),arg[3]))
        #         rslts += tmp
        # else:
        #     try:
        #         p = multiprocessing.Pool(np)
        #         tmp = p.map(findirets, args)
        #     finally:
        #         LOG.debug('  closing pool')
        #         p.close()
        #     rslts = reduce(iadd, tmp)
        cols = GGB.BEDCOLS[:6] # ['chr','st','ed','name','sc1','strand']
        self.stats['FINDIRETS.#irets_ex'] = len(rslts)
        if len(rslts)>0:
            df = PD.DataFrame(rslts, columns = cols)
            dfname = fn.write_bed(df, 'findirets.irets', ncols=6)
            LOG.info('{0} irets found'.format(len(df)))
            return df
        LOG.warning('******************** NO IRET FOUND!***********************')
        return PD.DataFrame(N.zeros((0,len(cols))),columns=cols) # return empty dataframe
               
def findirets(sj,me,irets,chrom,covratio,covth):
    mg = GP.MEGraph2(sj,me)
    # turn irets junction into iret exons
    def _iter_irets():
        for chrom,ratio,strand,st_id,ed_id,cov in UT.izipcols(irets,['chr','ovlratio','strand','st_id','ed_id','cov']):
            for st,namel,covl in mg.sj_leftex(st_id,flds=['st','name','cov']):
                for ed,namer,covr in mg.sj_rightex(ed_id,flds=['ed','name','cov']):
                    if ((2*cov > (covr+covl)*covratio)) & (cov>covth):
                    #if (cov > min(covr,covl)*covratio) & (cov>covth):
                        name = namel+'_iret_'+namer
                        yield (chrom,int(st),int(ed),name,ratio,strand)

    cols = GGB.BEDCOLS[:6] # ['chr','st','ed','name','sc1','strand']
    recs = [x for x in _iter_irets()]
    #LOG.debug('{0}: len:{1}'.format(chrom,len(recs)))
    return recs

class FINDSE(SUBASE):
    """Find single exons.

    Args:
        me: exon DataFrame
        edgefixer: EdgeFixer instance
        secovth: se coverage threshold

    Returns:
        :ae: all exons
        :se: single exons
        :me: multi-exons

    Related Parameters:
        * minsecovth: minimum single exon coverage (normalized to million alignments)
          default {minsecovth} (merging: {minsecovth_m})
        * secovth: default SE cov threshold if not using adaptive version
          default {secovth} (merging: {secovth_m})
        * se_gap: single exon gap fill, default {se_gap}
        * se_binth: coverage threshold for SE finding, default {se_binth} (merging: {se_binth_m})
        * se_sizeth: single exon size th, default {se_sizeth} (merging: {se_sizeth_m})

    TempFiles:
        * se.cov.tmp.txt.gz
        * se.cov.all.txt.gz
        * secov.txt.gz
        * se.bed.gz
        * exons.afterfindse.txt.gz
        * me.exons.bed.gz

    """

    
    def call(self):
        me = self.asm.me
        edgefixer = self.asm.edgefixer
        secovth = self.asm.secovth

        fn = self.fnobj
        pr = self.params
        st = self.stats
        
        secovth = max(pr['minsecovth'], secovth)
        st['FINDSE.secovth'] = secovth
        
        gap = pr['se_gap']
        binth = pr['se_binth']
        sesizeth = pr['se_sizeth']
        override = pr['override']

        fname = fn.bedname('findse.se')
        if override or (not os.path.exists(fname)):
            LOG.debug(' finding SE...')
            se0 = self.find_se(me, covth=secovth,gap=gap,binth=binth,sizeth=sesizeth)
        else:
            LOG.debug(' reading cached SE {0}...'.format(fname))
            se0 = GGB.read_bed(fname)

        UT.set_ptyp(se0)
        idx = se0['ptyp']=='s'
        seme = se0[~idx] # 3' or 5' exons

        # fix 3' exons
        LOG.debug(' fixing 3pr extension 2nd phase ...')
        seme2 = edgefixer.fixSEedge(me, seme, utr='3pr')
        seme2 = edgefixer.fixSEedge(me, seme2, utr='5pr')
        #seme2 = edgefixer.fixSEedge(me, seme2, utr='se')
        se = se0[idx]
        st['FINDSE.#se'] = len(se)
        cols = GGB.BEDCOLS[:6]+['ptyp']
        me2 = PD.concat([me[cols],seme2[cols]], ignore_index=True)
        ae = PD.concat([me2[cols],se[cols]], ignore_index=True)
        #return ae, se, me2
        fn.write_bed(ae, 'assemble.exons0', ncols=6)
        UT.set_ids(ae)
        LOG.info('write exons ...=> assemble.exons0.txt.gz')
        fn.write_txt(ae, 'assemble.exons0', fm='h')
        self.asm.ae = ae
        self.asm.se = se
        self.asm.me = me2

    def find_se(self, exons, covth=0.5,gap=50,binth=0,sizeth=200,override=False):
        # from bin0gap50 subtract ME choose >200bp and cov>1
        # covarage file
        LOG.debug(' finding SE candidates...')
        fn = self.fnobj
        pr = self.params
        st = self.stats

        fname = fn.txtname('se.cov.tmp')
        aname = fn.txtname('se.cov.all')
        override = pr['override']

        if (not override) and (os.path.exists(fname)):
            LOG.debug('  reading cached SECOV {0}...'.format(fname))
            secov = UT.read_pandas(fname)
            secov['len'] = secov['ed']-secov['st']
            secov = secov[secov['len']>sizeth]
        elif os.path.exists(aname):
            # use cov calculated at FINDSECOVTH 
            secov = UT.read_pandas(aname)
            secov['len'] = secov['ed']-secov['st']
            secov = secov[secov['len']>sizeth]
        else:
            # if not using FINDSECOV then just calculate len>sizeth
            LOG.debug('  calculating SECOV...')
            binfile = self.bw2bed(binth)
            gapfile = self.fillgap(gap, binfile)
            mefile = fn.write_bed(exons, 'findse.me.exons', ncols=3)
            # subtract me from gapfilled
            LOG.debug(' calculating coverage...')
            cname = fn.bedname('findse.gap-sub-me') #os.path.join(MD.SJEXDIR,sname+'.gap-sub-me.bed.gz')
            if override or (not os.path.exists(cname)):
                BT.bedtoolsubtract(gapfile, mefile, cname)
            df = GGB.read_bed(cname)
            df['len'] = df['ed'] - df['st']
            df = df[df['len']>sizeth]
            #secov = BW.calc_cov(df, bwfile, fname) # ~30-40 sec
            secov = CC.calc_cov_mp(df, fn.bwfile, fname, pr['np']) # here secovtmp is saved
            os.unlink(mefile)

        secov = secov[secov['cov']>pr['minsecovth']] 
        # use lower threshold here so that you don't lose non-SE (i.e. alternative 3',5' UTR)
        # but can't process all of SE candidate since it will take forever, need to restrict to 
        # reasonable amount => minsecovth
        # select SE with proper secovth later at SELECTSEME
        fn.write_txt(secov, 'findse.secov', fm='h')
        LOG.info('#candidate SE:{0}'.format(len(secov)))
        st['FINDSE.#candidate_se'] = len(secov)
        # see if additional filtering is necessary for intronic SE
        # ==> covth 0.5 seems to be OK, ~4K overall
        # long ones are 3'UTR ==> needs to be fixed
        LOG.debug('fixing 3pr extension...')
        cols = GGB.BEDCOLS[:6]
        recs = [x for x in self._fix_se_ext(exons, secov)]
        LOG.info('#candidate SE after _fix_se_ext:{0}'.format(len(recs)))
        se = PD.DataFrame(recs,columns=cols)
        st['FINDSE.#candidate_se_after_fix_se_ext'] = len(se)
        fn.write_bed(se, 'findse.se', ncols=6)
        return se
    
    def _fix_se_ext(self, exons, secov):
        # generator
        # go through se and see if it is attached to 5' or 3'of me exon
        # if so create extended exon
        # otherwise just spit out itself
        est = exons.set_index(['chr','st'])
        eed = exons.set_index(['chr','ed'])
        def _get(tgt,chrom,pos):
            try:
                test0 = tgt.ix[chrom]
                test1 = test0.ix[pos] # if no match => except
                return test0.ix[[pos]] # always return DataFrame
            except: # constructing DataFrame for Null cases seems to take a lot of time
                return None
        LOG.info('  processing {0} SE candidates'.format(len(secov)))
        for i,(chrom,st,ed,cov) in enumerate(UT.izipcols(secov, ['chr','st','ed','cov'])):
            le = _get(eed,chrom,st) # for 3utr + strand find exons on left eed==st
            re = _get(est,chrom,ed) # for 3utr - strand find exons on right est==ed
            if (le is not None) and (re is not None):
                pass
                # this should have been detected in iret
                # for nal,srl,e_st in le[['name','strand','st']].values:
                #     for nar,srr,e_ed in re[['name','strand','ed']].values:
                #         if srl==srr:
                #             name = nal+'|SE{0}|'.format(i)+nar
                #             yield (chrom,e_st,e_ed,name,cov,sr)            
            elif le is not None:# chr, st, ed, name, sc1, strand
                for nal,srl,e_st in le[['name','strand','st']].values:
                    name = nal+'|SE{0}'.format(i)
                    yield (chrom,e_st,ed,name,cov,srl)
            elif re is not None:
                for nar,srr,e_ed in re[['name','strand','ed']].values:
                    name = 'SE{0}|'.format(i)+nar
                    yield (chrom,st,e_ed,name,cov,srr)
            else:
                yield (chrom,st,ed,'SE{0}'.format(i),cov,'.')
                
    def _fix_se_3prext(self, exons, secov):
        # generator
        # go through se and see if it is attached to 3'of me exon
        # if so create extended 3pr end exon
        # otherwise just spit out itself
        est = exons.set_index(['chr','strand','st'])
        eed = exons.set_index(['chr','strand','ed'])
        def _get(tgt,chrom,strand,pos):
            try:
                return tgt.ix[chrom].ix[strand].ix[pos]
            except:
                return None
        for chrom,st,ed,cov in secov[['chr','st','ed','cov']].values:
            le = _get(eed,chrom,'+',st) # for 3utr + strand find exons on left eed==st
            re = _get(est,chrom,'-',ed) # for 3utr - strand find exons on right est==ed
            if (le is not None) or (re is not None):
                if le is not None:# chr, st, ed, name, sc1, strand
                    if len(le.shape)==1:
                        name = le['name']+'|SE'
                        yield (chrom,le['st'],ed,name,cov,'+')
                    else:
                        for j,e in le.iterrows():
                            name = e['name']+'|SE'
                            yield (chrom,e['st'],ed,name,cov,'+')
                if re is not None:
                    if len(re.shape)==1:
                        name = 'SE|'+re['name']
                        yield (chrom,st,re['ed'],name,cov,'-')
                    else:
                        for j,e in re.iterrows():
                            name = 'SE|'+e['name']
                            yield (chrom,st,e['ed'],name,cov,'-')
            else:
                yield (chrom,st,ed,'SE',cov,'.')

class FIND53IR(SUBASE):
    """Find 5',3' exon (cut), intron retentions and single exons.

    Args:
        me: exon DataFrame

    Returns:
        :ae: exons with 5',3' cuts, irets, and single exons.

    Related Parameters:
        * minsecovth: minimum single exon coverage (normalized to million alignments)
          default {minsecovth} (merging: {minsecovth_m})
        * secovth: default SE cov threshold if not using adaptive version
          default {secovth} (merging: {secovth_m})
        * se_gap: single exon gap fill, default {se_gap}
        * se_binth: coverage threshold for SE finding, default {se_binth} (merging: {se_binth_m})
        * se_sizeth: single exon size th, default {se_sizeth} (merging: {se_sizeth_m})
        * find53ir_covratio: cov ratio threshold for FIND53IR, 
          default {find53ir_covratio}, (merging: {find53ir_covratio_m})
        * find53ir_covth: cov threshold for FIND53IR
          default {find53ir_covth}, (mergin: {find53ir_covth_m})

    TempFiles:
        * se.cov.tmp.txt.gz
        * se.cov.all.txt.gz
        * me.exons.bed.gz
        * bw*.bed.gz
        * gap-sub-me.bed.gz
        * assemble.exons0.txt.gz
        * assemble.exons0.bed.gz
    """

    def call(self):
        me = self.asm.me
        fn = self.fnobj
        pr = self.params
        st = self.stats
        
        secovth = max(pr['minsecovth'], pr['secovth'])
        st['FINDSE.secovth'] = secovth
        
        override = pr['override']

        fname = fn.bedname('find53ir')
        if override or (not os.path.exists(fname)):
            LOG.debug(' finding 53IR...')
            ae = self.find(me, secovth)#,gap=gap,binth=binth,sizeth=sesizeth)
        else:
            LOG.debug(' reading cached 53IR {0}...'.format(fname))
            ae = GGB.read_bed(fname)
        #return ae
        self.asm.ae = ae
        fn.write_bed(ae, 'assemble.exons0', ncols=6)
        UT.set_ids(ae)
        LOG.info('write exons ...=> assemble.exons0.txt.gz')
        fn.write_txt(ae, 'assemble.exons0', fm='h')

    def find(self, me, secovth):
        # TODO fixedge dse as well?

        fn = self.fnobj
        pr = self.params
        st = self.stats
                
        gap = pr['se_gap']
        sizeth = pr['se_sizeth']
        override = pr['override']
        np = pr['np']

        # calc SE candidates (subtract ME) and calc cov
        fname = fn.txtname('se.cov.tmp')
        aname = fn.txtname('se.cov.all')
        if (not override) and os.path.exists(fname):
            LOG.info('  reading cached SECOV {0}...'.format(fname))
            secov = UT.read_pandas(fname)
            secov['len'] = secov['ed']-secov['st']
            secov = secov[secov['len']>sizeth]
        elif (not override) and os.path.exists(aname):
            LOG.info('  reading cached SECOV {0}...'.format(aname))
            # use cov calculated at FINDSECOVTH 
            secov = UT.read_pandas(aname)
            secov['len'] = secov['ed']-secov['st']
            secov = secov[secov['len']>sizeth]
        else:
            # if not using FINDSECOV then just calculate len>sizeth
            LOG.info('  calculating SECOV...')
            binfile = self.bw2bed(pr['se_binth'])
            gapfile = self.fillgap(gap, binfile)
            mefile = fn.write_bed(me, 'find53ir.me.exons', ncols=3)
            LOG.debug(' calculating coverage...')
            cname = fn.bedname('find53ir.gap-sub-me')
            if override or (not os.path.exists(cname)):
                BT.bedtoolsubtract(gapfile, mefile, cname)
            df = GGB.read_bed(cname)
            df['len'] = df['ed'] - df['st']
            df = df[df['len']>sizeth]
            secov = CC.calc_cov_mp(df, fn.bwfile, fname, np) # here secovtmp is saved
            os.unlink(mefile)

        # calc ME cov
        if override or not os.path.exists(fn.txtname('find53ir.cov')):
            LOG.debug('calculating ME cov...')
            UT.set_ids(me)
            me = CC.calc_cov_ovl_mp(
                    srcname=me, 
                    bwname=fn.bwfile, 
                    dstname=fn.txtname('find53ir.cov'), 
                    np=np,
                    override=override,
                    covciname=fn.txtname('find53ir.covci'),
                    ciname=fn.txtname('find53ir.ci'))
        else:
            me = fn.read_txt('find53ir.cov')

        # write BED files and find attaching ME/SE 
        secov['st-1'] = secov['st']-1
        secov['ed+1'] = secov['ed']+1
        secov['name'] = ['SE{0}'.format(x) for x in N.arange(len(secov))]
        secov['_id2'] = N.arange(len(secov))
        secov['strand'] = '.'
        secols = ['chr','st-1','ed+1','cov','_id2','name','st','ed']
        mecols = ['chr','st','ed','name','_id','cov','strand']
        # somehow ['chr','st','ed','name','_id','strand','cov'] this order segfaults
        # when intersected with others
        mecols2 = ['echr','est','eed','ename','eid','ecov','strand']
        a = UT.write_pandas(secov[secols],fn.bedname('find53ir.se'),'')
        b = UT.write_pandas(me[mecols],fn.bedname('find53ir.me'),'')
        c = fn.txtname('find53ir.ovl')
        c = BT.bedtoolintersect(a,b,c,wao=True)
        cols = secols+mecols2+['ovl']
        d = UT.read_pandas(c, names=cols)
        d['attachleft'] = d['ed+1']-1==d['est']
        d['attachright'] = d['st-1']+1==d['eed']
        d['bound'] = (d['ename'].str.startswith('['))&(d['ename'].str.endswith(']'))
        # SE
        dse = d[d['echr']=='.'].copy() # no overlap
        dse = dse[dse['cov']>secovth]
        dse['sc1'] = dse['cov'] # for writing BED
        # ME
        dme = d[(d['attachleft']&d['bound'])|(d['attachright']&d['bound'])].copy()
        dme['covratio'] = dme['cov'].astype(float)/dme['ecov'].astype(float)
        dme = dme[dme['covratio']>pr['find53ir_covratio']] # and substantial portion

        # attachleft sid, attachright sid, iret side, etc.
        def _lri(strand):
            al0p = set(dme[(dme['attachleft'])&(dme['strand']==strand)]['_id2'].values)
            ar0p = set(dme[(dme['attachright'])&(dme['strand']==strand)]['_id2'].values)
            irp = al0p.intersection(ar0p)
            alp = al0p.difference(ar0p) 
            arp = ar0p.difference(al0p) 
            return alp,arp,irp
        alp,arp,irp = _lri('+')
        aln,arn,irn = _lri('-')
        ir = irp.union(irn)
        al = alp.union(aln)
        ar = arp.union(arn)
        LOG.debug('#al({0}) #ar({1}) #iret({2})'.format(len(al),len(ar),len(ir)))

        # fix edges, no chop necessary but need copy
        sei = secov.set_index('_id2')
        # cicols: ['chr','st','ed','name','sc1', 'direction']
        # mecols: ['chr','st','ed','name','strand','_id2']
        # ci.name == str(me._id2)
        def _calc(id2s,direction,func):
            if len(id2s)==0:
                return None
            me = sei.ix[id2s].reset_index()[['chr','st','ed','name','strand','_id2']].copy()
            ci = me[['chr','st','ed']].copy()
            ci['name'] = me['_id2'].astype(str)
            ci['sc1'] = me['_id2']
            ci['direction'] = direction
            return self.process_mp(ci,me,func)
        alr = _calc(al,'-',fixedges_m)
        arr = _calc(ar,'+',fixedges_m)
        irr = _calc(ir,'.',cutedges_m)
        # returned ['chr','st','ed','name','sc1','strand','_id2']

        # attach exons
        dmi = dme.set_index('_id2')
        i2c = UT.df2dict(secov, '_id2','cov')
        def _makedics(id2s, boolcol):
            t = dmi.ix[id2s].reset_index()
            t = t[t[boolcol]]
            tg = t.groupby('_id2')['eid'].apply(lambda x: list(x)).reset_index()
            s2eids = UT.df2dict(tg, '_id2','eid')
            tg2 = t.groupby('eid').first()
            e2row = dict(zip(tg2.index.values, UT.izipcols(tg2,['est','eed','ename','strand'])))
            return s2eids, e2row
        def _algen():
            if alr is None:
                return            
            s2eids,e2row = _makedics(al,'attachleft')
            for c,s,e,n,i2 in UT.izipcols(alr,['chr','st','ed','name','_id2']):
                # SE|[mel] # (chr,sst,eed,sn|en,i2,estrand)
                cov = i2c[i2]
                for eid in s2eids[i2]: # erow: st,ed,name,strand
                    erow = e2row[eid]
                    name = n+'|'+erow[2]
                    yield (c,s,erow[1],name,cov,erow[3])
        def _argen():
            if arr is None:
                return            
            s2eids,e2row = _makedics(ar,'attachright')
            for c,s,e,n,i2 in UT.izipcols(arr,['chr','st','ed','name','_id2']):
                # [mer]|SE, (chr,est,sed,name,i2,estrand)
                cov = i2c[i2]
                for eid in s2eids[i2]: # erow: st,ed,name,strand
                    erow = e2row[eid]
                    name = erow[2]+'|'+n
                    yield (c,erow[0],e,name,cov,erow[3])
        def _irgen():
            if irr is None:
                return
            s2le,le2row = _makedics(ir,'attachleft')
            s2re,re2row = _makedics(ir,'attachright')
            for c,s,e,n,i2 in UT.izipcols(irr,['chr','st','ed','name','_id2']):
                # [mer]|SE|[mel]
                cov = i2c[i2]
                for eidl in s2le[i2]:
                    erowl = le2row[eidl]
                    for eidr in s2re[i2]:
                        erowr = re2row[eidr]
                        if erowl[3]==erowr[3]:
                            name = erowr[2]+'|'+n+'|'+erowl[2]
                            yield (c,erowr[0],erowl[1],name,cov,erowl[3])
        alrecs = [x for x in _algen()]
        arrecs = [x for x in _argen()]
        irrecs = [x for x in _irgen()]
        cols = ['chr','st','ed','name','sc1','strand'] # se cov in 'sc1'
        medf = PD.DataFrame(alrecs+arrecs+irrecs,columns=cols)
        ex = PD.concat([dse[cols],medf[cols]], ignore_index=True)
        fn.write_bed(ex,'find53ir',ncols=6)

        UT.set_ptyp(ex)
        UT.set_ptyp(me)
        idxse = ex['ptyp']=='s'
        st['FINDSE.#se'] = N.sum(idxse)
        cols = GGB.BEDCOLS[:6]+['ptyp']
        me['sc1'] = me['cov']
        ae = PD.concat([me[cols],ex[cols]], ignore_index=True)
        ae = ae.groupby(['chr','st','ed','strand']).first().reset_index()
        return ae

    def process_mp(self, ci, me, func):
        pr = self.params
        bwname = self.fnobj.bwfile
        args = []
        for chrom in self.chroms(me):
            mechr = me[me['chr']==chrom][['chr','st','ed','name','strand','_id2']].copy()
            cichr = ci[ci['chr']==chrom].copy()
            args.append((cichr,mechr,bwname,pr,chrom))
        rslts = UT.process_mp(func, args, pr['np'])
        # rslts = []
        # if pr['np']==1:
        #     for i,arg in enumerate(args):
        #         tmp = func(*arg)
        #         LOG.debug('  processing {3}: {0}/{1} {2}...'.format(i+1,len(args),len(tmp),arg[-1]))
        #         rslts += tmp
        #     return PD.concat(rslts,ignore_index=True)
        # try:
        #     p = multiprocessing.Pool(pr['np'])
        #     #tmp = p.map(cutedges_m, args)
        #     tmp = p.map(func, args)
        # finally:
        #     LOG.debug('  closing pool')
        #     p.close()
        # return PD.concat(tmp, ignore_index=True)
        return PD.DataFrame(rslts, columns=['chr','st','ed','name','sc1','strand','_id2'])

class FINDSECOVTH(SUBASE):
    """Find single exon coverage threshold (secovth).

    Args:
        me: exon DataFrame

    Returns:
        :secovth: threshold

    Related Parameters:
        * findsecovth_useref: whether to use ref for finding secov, if not use ME
          default {findsecovth_useref}
        * minsecovth: minimum single exon coverage (normalized to million alignments)
          default {minsecovth} (merging: {minsecovth_m})
        * se_gap: single exon gap fill, default {se_gap}

    Files:
        * findsecovth.params.txt.gz
        * findsecovth.pdf

    TempFiles:
        * findsecovth.refex.covci.txt.gz
        * findsecovth.refex.cov.txt.gz
        * se.cov.all.txt.gz
        * irets.me.cov.txt.gz
        * me.covci.txt.gz
        * me.ci.txt.gz
        * findsecovth.gap-sub-me.bed.gz


    """

    bins = N.arange(0.4,1.8,0.05)
    """ gamma values to search """

    fitrange = (6,9)
    """ x range to fit a line (fitmax,resmax)"""
    
    def call(self):
        ex = self.asm.me
        fn = self.fnobj
        pr = self.params
        st = self.stats
        self.ex = ex  
        if pr['findsecovth']:
            if pr['findsecovth_useref'] and fn.refgtf.exists():
                self.calc_refexcov()
            self.calc_allsecov(ex)
            secovth = self.find_secovth()
            self.save_params()
        else:
            secovth = pr['secovth']
        secovth = max(pr['minsecovth'], secovth)
        st['FINDSECOVTH.secovth'] = secovth
        self.info = 'SECOVTH={0}'.format(secovth)
        #return secovth
        self.asm.secovth = secovth
    
    def save_params(self):
        fn = self.fnobj
        st = {}
        # gamma, a1
        # residue, th99 for (ref, me, se)
        attrs = ['gamma','a1','ref_res','ref_th99bn','me_res','me_th99bn','se_res',
                'ref_nf','me_nf','se_nf','se_th99bn','se_th99']
        pr = dict([(x, getattr(self, x)) for x in attrs])
        for a in attrs:
            st['FINDSECOVTH.'+a] = pr[a]
        prdf = PD.DataFrame(pr,index=[fn.sname]).T
        fname = fn.fname('findsecovth.params.txt', category='stats')
        UT.write_pandas(prdf, fname, fm='ih')
        
    def calc_refexcov(self):
        fn = self.fnobj
        pr = self.params

        excovname = fn.txtname('findsecovth.refex.cov')
        if  os.path.exists(excovname):
            rsj,rex = fn.refsjex()
            self.rsj = rsj
            self.rex = UT.read_pandas(excovname)
        else:
            rsj,rex = fn.refsjex()
            self.rsj = rsj
            self.rex = CC.calc_cov_ovl_mp(
                srcname=rex, 
                bwname=fn.bwfile, 
                dstname=excovname, 
                np=pr['np'], 
                covciname=fn.txtname('findsecovth.refex.covci'), 
                ciname=fn.refname('ci'), 
                colname='cov', 
                override=pr['override']
            )
    
    def calc_allsecov(self, me):
        fn = self.fnobj
        pr = self.params

        aname = fn.txtname('se.cov.all')
        if not os.path.exists(aname):
            self._calc_allsecov(me)
        LOG.debug('calculating ME cov...')
        mname = fn.txtname('irets.me.cov')

        if os.path.exists(mname):
            # use previously calculated cov (w/o irets)
            # since only used for making figs/fitting w/o irets is OK
            self.ex = UT.read_pandas(mname)
        else:
            UT.set_ids(me)
            self.ex = CC.calc_cov_ovl_mp(
                srcname=me, 
                bwname=fn.bwfile, 
                dstname=mname,
                np=pr['np'], 
                covciname=fn.txtname('findsecovth.me.covci'), 
                ciname=fn.txtname('findsecovth.me.ci'), 
                colname='cov', 
                override=pr['override']
            )        

    def _calc_allsecov(self, me):
        LOG.debug(' finding SE candidates...')
        fn = self.fnobj
        pr = self.params

        binfile = self.bw2bed(pr['binth'])
        gapfile = self.fillgap(pr['se_gap'], binfile) 
        # pr['gap'] 50 vs. pr['se_gap'] 170 makes huge difference...
        mefile = fn.write_bed(me, 'findsecovth.me.exons', ncols=3)
        # subtract me from gapfilled
        LOG.debug(' calculating coverage...')
        cname = fn.bedname('findsecovth.gap-sub-me')
        cname = BT.bedtoolsubtract(gapfile, mefile, cname)
        df = GGB.read_bed(cname)
        fname = fn.txtname('se.cov.all')
        #secov = BW.calc_cov(df, bwfile, fname) # slow ~500K targets
        secov = CC.calc_cov_mp(bed=df, bwname=fn.bwfile, fname=fname, np=pr['np'])
        #os.unlink(cname)
    
    def _calc_th99(self, xf,yo,yf,gamma,i0=0,i1=40):
        idx = slice(i0,i1)
        cp = N.sum(2**yf[idx]-gamma)   # condition positive
        cn = N.sum(2**yo[idx]-gamma)-cp # condition negative
        if cn<0:
            cn = 0
        fn = N.cumsum(2**yf[idx]-gamma) # false negative (condition positive but < th)
        tp = cp - fn # true positive (condition positive and >= th)
        tn = N.cumsum(2**yo[idx]-gamma) - fn  # true negative
        tn[tn<0]=0
        fp = cn - tn
        fp[fp<0]=0
        tpr = tp/cp
        fpr = fp/cn
        fpr[N.isnan(fpr)]=0
        tf = tpr-fpr
        tmp = N.nonzero(fpr<=0.01)[0]
        th99x = xf[N.min(tmp)]
        th99 = 2**th99x - gamma
        return th99x,th99

    def _fit_hist(self, refcov, gamma, a1=None, ax=None, title='', nf=1.):
        # find best gamma for gen4
        fma,xma = self.fitrange
        fmi = N.log2(gamma)
        nbins = float(max(min(100, len(refcov)/50), 50))
        width = (xma-fmi)/nbins
        bins = N.arange(fmi,xma,width)
        tot = float(len(refcov))
        v = N.log2(refcov+gamma)
        h,b = N.histogram(v, bins=bins)
        x = b[1:]
        y = N.log2(h+1) # logscale
        # st,ed = 0,60 # range of fit
        st = 0 
        ed = int(nbins*(float(fma-fmi)/(xma-fmi)))
        xf = x[0:ed]
        yo = y[0:ed]
        if a1 is None:
            x3 = x[st:ed]
            y3 = y[st:ed]
            a1,a0 = N.polyfit(x3,y3,1)
        else:
            LOG.info('ed={0} x[ed]={1}'.format(ed, x[ed]))
            # st = 15 # throw away initial 14 points
            # x3 = x[st:ed]
            # y3 = y[st:ed]
            # a0 = N.mean(y3-a1*x3) # directly calculate intercept
            # change st and calculate ave_res take a0 for min ave_res
            use = int(nbins*0.4) if nbins < 50 else 20
            sts = range(st,ed-use)
            mres = []
            a0s = []
            for sttmp in sts:
                x3 = x[sttmp:ed]
                y3 = y[sttmp:ed]
                a0tmp = N.mean(y3-a1*x3)
                yftmp = a0tmp+a1*x3
                a0s.append(a0tmp)
                mres.append(N.mean((y3-yftmp)**2))
            min_res = N.min(mres)
            min_idx = mres.index(min_res)
            a0 = a0s[min_idx]
            min_st = sts[min_idx]
            LOG.info('FINDSECOV: SE fit a0={0}, st={1}, x[st]={2}'.format(a0,min_st,x[min_st]))

        yf = a0+a1*xf
        res = N.sum((yo-yf)**2)
        
        if ax is not None: # also calculate th99
            th99x,th99 = self._calc_th99(xf,yo,yf,gamma,0,40)

            ax.plot(x,y,'o',ms=3)
            ax.plot(x,a0+a1*x,'m-',lw=1)
            ax.plot([th99x,th99x],[0,16],'r--',lw=0.5)
            ax.plot([x[ed],x[ed]],[0,8],'k--',lw=0.5)
            #ax.plot([x[0],x[0]],[0,a0+a1*x[0]],'k--',lw=0.5)
            ax.set_title(title)
            ax.set_xlabel('log2(cov+gamma)')
            ax.set_ylabel('log2(count+1)')
            ax.text(4,14.5,'res={:.1f}'.format(res))
            ax.text(4,13.,'slope={:.2f}'.format(a1))
            ax.text(4,11.5,'ncovth={:.1f}'.format(th99))
            ax.text(4,10.,'covth={:.1f}'.format(th99*nf))
            ax.text(4,8.5,'nf={:.1f}'.format(nf))
            ax.set_ylim([0,17])
            
            return res, th99
        # [TODO] return types are different (confusing) => separate into two methods
        return gamma,a1,a0,res
    
    def find_gamma_slope(self):
        gammas = self.bins
        refcov, ref_nf = self._get_refcov()
        cols = ['gamma','a1','a0','res']
        self.refcovfits = fits = PD.DataFrame([self._fit_hist(refcov, x, nf=ref_nf) for x in gammas], columns=cols)
        fits['-gamma'] = -fits['gamma']
        gamma,a1 = fits.sort_values(['res','-gamma']).iloc[0][['gamma','a1']]
        self.gamma=gamma
        self.a1 = a1
        return gamma, a1

    def _get_refcov(self):
        pr = self.params
        fn = self.fnobj
        if pr['findsecovth_useref'] and fn.refgtf.exists():
            cdf = self.rex
        else: # use ME cov
            ex = self.ex
            if 'cat' not in ex.columns:
                sj = self.asm.sj
                UT.set_exon_category(sj,ex)
            if 'len' not in ex.columns:
                ex['len'] = ex['ed']-ex['st']
            cdf = ex[ex['cat']!='s']
        cdf = cdf[cdf>0] # only use expressed 
        # normalize 1e4 factor is to match to 100bp readlen normalized to 1M reads
        self.ref_nf = ref_nf = self.calc_normfactor(cdf) # N.sum(refcov)/1e4
        self.refcov = refcov = cdf['cov'].values/ref_nf
        return refcov, ref_nf

    def calc_normfactor(self, cdf):
        # cdf['cov', and 'len',or 'st','ed']
        # calculate total bp = sum(cov*len) ~ aligned*100
        # so normalizing to aligned/1M ~ sum(cov*len)/1e8
        if not 'len' in cdf.columns:
            cdf = cdf[['cov','st','ed']].copy()
            cdf['len'] = cdf['ed']-cdf['st']
        totbp = (cdf['cov']*cdf['len']).sum()
        return totbp/1e8

    def find_secovth(self):
        fn = self.fnobj
        pr = self.params

        if pr['findsecovth_useref'] and fn.refgtf.exists():
            useref = True
        else:
            useref = False
        
        # find gamma, slope from reference coverage histogram
        gamma, a1 = self.find_gamma_slope()
        self.gamma = gamma
        self.a1 = a1
        
        # use gamma, slope to find secovth
        fig, axr = P.subplots(2,2,figsize=(6,6))
        # fig, axr = P.subplots(2,3,figsize=(9,6))
        P.subplots_adjust(wspace=0.4,hspace=0.4,top=0.9)
        fig.suptitle(fn.sname)
        # panel1 refcov fits
        fits = self.refcovfits
        ax = axr[0][0]
        ax.plot(fits['gamma'].values, fits['res'].values, '.-')
        ylim = ax.get_ylim()
        ax.plot([gamma,gamma],ylim, 'r--', lw=0.5)
        ypos=ylim[0]+0.9*(ylim[1]-ylim[0])
        ax.text(1,ypos, 'min at {0:.2f}'.format(gamma))
        ax.set_xlabel('gamma')
        ylabel = 'REF cov residual' if useref else 'sample ME cov residual'
        ax.set_ylabel(ylabel)
        ax.set_title('finding gamma and slope')
        
        # panel2 refcov optimum
        ax = axr[0][1]
        refcov, ref_nf = self._get_refcov()
        title = 'ref ME & SE' if useref else 'sample ME'
        self.ref_res,self.ref_th99bn = self._fit_hist(refcov, gamma, a1, ax=ax, title=title, nf=ref_nf)
        
        # refcov SE
        # if pr['findsecovth_useref'] and fn.refgtf.exists():
        #     rsj,rex = self.rsj, self.rex
        #     if 'cat' not in rex.columns:
        #         UT.set_exon_category(rsj,rex)
        #     rse = rex[rex['cat']=='s']
        #     rse = rse[rse['cov']>0]['cov'].values
        #     rse_nf = N.sum(rse)/1e4
        #     self.rse = rse = rse/rse_nf
        #     ax = axr[0][2]
        #     tmp = self._fit_hist(rse, gamma, a1, ax=ax, title='ref SE', nf=rse_nf)

        # panel3 me fits
        ax = axr[1][0]
        ex = self.ex
        if 'cat' not in ex.columns:
            sj = self.asm.sj
            UT.set_exon_category(sj,ex)
        if 'len' not in ex.columns:
            ex['len'] = ex['ed']-ex['st']
        me = ex[ex['cat']!='s']
        me = me[me['cov']>0]
        # normalize
        self.me_nf = me_nf = self.calc_normfactor(me) # N.sum(me)/1e4
        self.mecov = mecov = me['cov'].values/me_nf
        self.me_res,self.me_th99bn = self._fit_hist(mecov, gamma, a1, ax=ax,title='ME', nf=me_nf)
        
        # panel4 seall
        ax = axr[1][1]
        se = fn.read_txt('se.cov.all')
        if 'len' not in se.columns:
            se['len'] = se['ed']-se['st']
        se = se[se['cov']>0]
        # normalize
        self.se_nf = se_nf = self.calc_normfactor(se) #N.sum(secov)/1e4
        self.secov = secov = se['cov'].values/se_nf
        self.se_res,self.se_th99bn = self._fit_hist(secov, gamma, a1, ax=ax,title='SE candidates', nf=se_nf)
        self.se_th99 = se_nf*self.se_th99bn
        fname = fn.fname('findsecovth.pdf', category='stats')
        try:
            fig.savefig(fname)
        except:
            LOG.warning('failed to save secovthfig')
        return self.se_th99
                
class SELECTSEME(SUBASE):
    """Remove noise.

    Args:
        sj: junction dataframe
        ae: exon (all exon = single+multi exons) dataframe

    Returns:
        :sj: junction dataframe without spurious junctions
        :ae: exon dataframe without spurious exons

    Related Parameters:
        * se_sizeth2: remove SE with length smaller than this, default {se_sizeth2}
        * np: number of CPU to use, default {np}
        * me2exon_sjth: 2exon genes with splice junction support less than this will be removed
          default {me2exon_sjth}
        * me2exon_sizeth: 2exon genes terminal size th ~ 2 x read length, default {me2exon_sizeth}
        * me2exon_covth: 2exon genes cov th, default {me2exon_covth} 


    Files:
        * sj.txt.gz
        * ex.txt.gz

    TempFiles:
        * selectseme.rm.se.bed.gz
        * selectseme.rm.me.bed.gz
        * selectseme.rm.ovl.txt.gz

    """
    
    def call(self):
        sj = self.asm.sj
        ae = self.asm.ae
        fn = self.fnobj
        pr = self.params
        
        # select SE cov > secovth
        secovth = self.stats.get('FINDSE.secovth', pr['secovth'])
        ae['len'] = ae['ed']-ae['st']
        me, se = UT.mese(ae)
        se2 = se[(se['cov']>secovth)&(se['len']>pr['se_sizeth2'])]
        # sesizeth 50 (shorter in case of concatenation), sesizeth2 200 (real size threshold)
        #secov thresholding is also done at "assemble" FINDSE part, but there may be ME origin SE (irets)
        sj1, me1 = self._remove_bad_2exon_genes(sj, me)
        sj2, ex2 = self._remove_overlapping_se(sj1, se2, me1)
        # ~400 SE (~15%) overlaps with ME
        # this will remove irets for 2 exon genes, which often can be noise
        # so removing corresponding iret will lose true model
        # only benefit of doing this is on the exon-match index, so not necessary
        #sj2 = sj1
        #ex2 = PD.concat([me1,se2], ignore_index=True)
        ex3 = ex2.groupby(['chr','st','ed','strand']).first().reset_index()
        
        # SAVE
        fn.write_txt(sj2, 'sj', category='output')
        fn.write_txt(ex3, 'ex', category='output')
        # if pr['findsecovth']: # save exname3
        #     fn.write_txt(ex3, 'ex2', category='output')
        # else: # save in exname2
        #     fn.write_txt(ex3, 'ex', category='output')
        #return sj2, ex3
        self.asm.sj = sj2
        self.asm.ae = ex3
    
    def _remove_bad_2exon_genes(self, sj, ex):
        # [TODO] also remove ME (>2 exons) with low SJ support
        # e.g. 3 exon genes but uniq SJ support <2 etc. 
        
        pr = self.params
        # find 2 exon genes
        exg = ex.groupby('_gidx')
        exgs = exg.size()
        exons = ex.set_index('_gidx').ix[exgs[exgs==2].index].reset_index()
        LOG.debug('  number of 2 exon exons:{0}'.format(len(exons)))
        self.stats['SELECTSEME_num_2exonexons'] = len(exons)
        # find junction support
        exons = exons.sort_values(['_gidx','a_id'])
        aids = exons.groupby('_gidx',sort=False)['a_id'].last() # -1,... 
        jsup = sj.set_index('a_id').ix[aids.values]['sc1']
        jsup = jsup.groupby(jsup.index).sum()
        aid2jsup = dict(zip(jsup.index, jsup.values))
        jsup = PD.Series([aid2jsup[x] for x in aids.values], name='jsup', index=aids.index)
        # find total length
        tlen = exons.groupby('_gidx')['len'].sum()
        # find coverages
        e = exons.set_index('_gidx')
        val = e['len']*e['cov']
        #tval = val.groupby(val.index).sum()
        #tcov = tval/tlen
        tval = val.groupby(val.index).min() # in case one overlaps with an exon of another gene
        tcov = 2*tval/tlen
        # which ones to remove?
        sjth = pr['me2exon_sjth'] # 1 or 2
        lenth = pr['me2exon_sizeth']
        covth = pr['me2exon_covth']
        idx0 = jsup<sjth # if less than two uniq junction for 2 exon gene remove regardless of len/cov
        idx1 = (jsup==sjth)&((tlen<lenth)|(tcov<covth))
        idx = idx0 | idx1
        r_aids = aids[idx].values
        r_gids = idx[idx].index.values
        ex1 = ex[~ex['_gidx'].isin(r_gids)]
        #sj1 = sj[~sj['a_id'].isin(r_aids)]
        sj1 = sj
        LOG.info('ex({0}) sj({1}) ==> ex({2}) sj({3})'.format(len(ex),len(sj),len(ex1),len(sj1)))
        self.stats['SELECTSEME.#bad2exons'] = len(ex)-len(ex1)
        return sj1, ex1
        
    def _remove_overlapping_se(self, sj, se, me):
        # remove SE overlapping ME exons
        # 1. write SE,ME bed
        # 2. bedtools intersect
        # 3. remove overlapping SE
        fn = self.fnobj
        st = self.stats
        
        cols = ['chr','st','ed','name','_id','strand']
        a = UT.write_pandas(se[cols], fn.bedname('selectseme.rm.se'), fm='')
        b = UT.write_pandas(me[cols], fn.bedname('selectseme.rm.me'), fm='')
        c = BT.bedtoolintersect(a,b,fn.txtname('selectseme.rm.ovl'),wao=True) # get original entry
        cdf = BT.read_ovl(c, GGB.BEDCOLS[:6])
        if len(cdf)>0:
            cdfg = cdf.groupby('sc1') # sc1 <= _id
            ovl = cdfg['ovl'].sum()
            siz = cdfg.size()
            nonovl_se = ovl[(ovl==0)|(siz>1)].index.values 
            # remove overlapping to only one ME exon ~((ovl>0)&(siz==1))
            se2 = se.set_index('_id').ix[nonovl_se].reset_index()
        else: # no overlap
            se2 = se
        ex2 = PD.concat([me,se2], ignore_index=True)
        LOG.info('se({0}) => se({1}) {2} removed'.format(len(se),len(se2),len(se)-len(se2)))
        st['SELECTSEME.#removed_se'] = len(se)-len(se2)
        return sj, ex2

def fixedge2(ex, find, params, chrom):
    ed = EdgeDetector(**params)
    # make sure first 3 cols are chr,st,ed
    cols = list(ex.columns)
    ist = cols.index('st')
    ied = cols.index('ed')
    icov = cols.index('cov')
    def _gen():
        with ed:
            if find=='drop': # left (st) => right (ed)
                for e in ex.values:
                    edges = ed.find(e[0],e[1],e[2],find,False,e[icov]) # [(interval, cov),...]
                    for (x0,x1), cov in edges:
                        e1 = e.copy()
                        e1[ied] = x1
                        e1[icov] = cov
                        yield e1
            else: # right => left
                for e in ex.values:
                    edges = ed.find(e[0],e[1],e[2],find,False,e[icov]) # [(interval, cov),...]
                    for (x0,x1), cov in edges:
                        e1 = e.copy()
                        e1[ist] = x0
                        e1[icov] = cov
                        yield e1
    tmp = [e for e in _gen()]
    #msg = '{2} found {0}/{1}'.format(len(tmp),len(ex),chrom)
    #LOG.debug(msg)
    #open('fixedge2-{0}.txt'.format(chrom),'w').write(msg)
    return tmp

class FIXEDGES2(SUBASE):
    """Find genes (connected components).

    Args:
        sj: junction dataframe
        ae: exon (all exon = single+multi exons) dataframe

    Related Parameters:
        * np: number of CPU to use, default {np}

    TempFiles:
        * genegraph-*
        * genes.cache.pic

    """

    # look at cov and detect sharp drop/rise
    
    def call(self):
        ex = self.asm.ae
        # clean up 3', 5' ends by detecting sharp drop and removing them if covratio is < threshold
        fn = self.fnobj
        pr = self.params

        cache = fn.txtname('fixedge2')
        if (not pr['override']) and (os.path.exists(cache)):
            return UT.read_pandas(cache)

        self.ex = ex
        idx3 = ex['cat']=='3'
        idx5 = ex['cat']=='5'
        idxp = ex['strand']=='+'
        idxn = ~idxp
        idxd = (idx3&idxp)|(idx5&idxn) # => find sharp drop
        idxr = (idx3&idxn)|(idx5&idxp) # <= find sharp rise
        self.exd = exd = ex[idxd]
        self.exr = exr = ex[idxr]
        self.exo = exo = ex[(~idx3)&(~idx5)] # 'i', 's'
        self.nexd = nexd = self.process_mp(exd,'drop')
        self.nexr = nexr = self.process_mp(exr,'rise')
        nex = PD.concat([exo,nexd,nexr],ignore_index=True)
        nexg = nex.groupby(['chr','st','ed','strand']).first().reset_index() # remove dup, in terms of pos&strand
        nexg['_id0'] = nexg['_id'] # old id used for gene identification
        UT.set_ids(nexg) # since there are dups
        self.nexg = nexg # for debug
        UT.write_pandas(nexg, cache, 'h')
        #return nexg
        self.info = '#ae:{0}'.format(len(nexg))
        self.stats['FIXEDGES2.#ae'] = len(nexg)
        self.asm.ae = nexg
        
    def process_mp(self, ex, find):
        pr = self.params
        fn = self.fnobj
        LOG.debug(' preparing data...')
        args = []
        params = dict(bwname=fn.bwfile,
                      sigmath=pr['ed_sigma'],
                      minth=pr['ed_minth'],
                      delta=pr['ed_window'],
                      covratio=pr['ed_covratio'],
                      winsize=pr['ed_window'])
        for chrom in self.chroms(ex):
            exc = ex[ex['chr']==chrom]
            args.append((exc, find, params,chrom))
        rslts = UT.process_mp(fixedge2, args, pr['np'])
        # rslts = []
        # np = p['np']
        # if np==1:
        #     for i,arg in enumerate(args):
        #         LOG.debug('  processing {0}/{1}...'.format(i+1,len(args)))
        #         rslts += fixedge2(arg)
        # else:
        #     try:
        #         p = multiprocessing.Pool(np)
        #         tmp = p.map(fixedge2, args)
        #     finally:
        #         LOG.debug('  closing pool')
        #         p.close()
        #     rslts = reduce(iadd, tmp)
        return PD.DataFrame(rslts, columns=ex.columns)
        
# class ADDJIE(SUBASE):
    # def call(self):
        # me = self.asm.me
        # sj = self.asm.sj
        # jie = self.asm.jie
        # fn = self.fnobj
        # pr = self.params        
        # # 1. find overlapping exons (discard jie if not found)
        # # 2. find surrounding junctions
        # # 3. select jie by ratio to surrounding junctions
        # a = GGB.write_bed(jie, fn.bedname('jie.sj'), ncols=6)
        # b = GGB.write_bed(me, fn.bedname('jie.ex'), ncols=6)
        # c = fn.txtname('jie.inte')
        # c = BT.bedtoolintersect(a,b,c,wao=True,f=1) # completely included
        # cols0 = GGB.BEDCOLS[:6]
        # cols = cols0+['b_'+x for x in cols0]+['ovl']
        # o = UT.read_pandas(c, names=cols)
        # o['len'] = o['ed'] - o['st'] # should be == ovl
        # UT._set_pos(sj, 'd_pos', 'st-1', 'ed')
        # UT._set_pos(sj, 'a_pos', 'ed', 'st-1')
        # UT._set_pos(o, 'd_pos', 'b_ed', 'b_st')
        # UT._set_pos(o, 'a_pos', 'b_st', 'b_ed')
        # # a_pos => max jcnt 
        # # TODO fix from sc1 to appropriate column
        # sjg = sj.groupby('a_pos')['sc1'].max().reset_index()
        # a2maxj = UT.df2dict(sjg, 'a_pos', 'sc1')
        # sjg = sj.groupby('d_pos')['sc1'].max().reset_index()
        # d2maxj = UT.df2dict(sjg, 'd_pos', 'sc1')
        # o['d_max'] = [d2maxj.get(x,0) for x in o['d_pos']]
        # o['a_max'] = [a2maxj.get(x,0) for x in o['a_pos']]
        # o['sjmax'] = o[['d_max','a_max']].max(axis=1)
        # o['locus'] = UT.calc_locus_strand(o)
        # jie['locus'] = UT.calc_locus_strand(jie)
        # l2j = UT.df2dict(jie, 'locus', 'sc1')
        # o['jcnt'] = [l2j[x] for x in o['locus']]
        # o['ratio'] = o['jcnt']/o['sjmax']
        # o = o[(o['strand']==o['b_strand'])&(o['len']==o['ovl'])&(o['ratio']>pr['jieratio'])]
        # usedjie = o.groupby(['chr','st','ed','strand']).first().reset_index()[cols0]
        # def _gen():
        #     cols = ['chr','st','ed','strand','name','b_st','b_ed','b_name']
        #     for c,s,e,t,n,bs,be,bn in UT.izipcols(o,cols):
        #         if bs<s:
        #             yield (c,bs,s,bn+'/'+n,0,t)
        #         if e<be:
        #             yield (c,e,be,n+'/'+bn,0,t)
        # ie = PD.DataFrame([x for x in _gen()], columns=cols0)
        # ex = PD.concat([me[cols0], ie[cols0]], ignore_index=True)
        # ex = ex.groupby(['chr','st','ed','strand']).first().reset_index()
        # sj1 = PD.concat([sj[cols0], usedjie[cols0]], ignore_index=True)
        # sj1 = sj1.groupby(['chr','st','ed','strand']).first().reset_index()
        # LOG.info('add junctions in exon...')        
        # _sttime = time.time()
        # n0 = len(self.sj)
        # n1 = len(self.me)
        # LOG.info(' sj {0}=>{1} ex {2}=>{3}'.format(n0, len(self.sj), n1, len(self.me)))
        # LOG.info(' time: {0:.3f}s'.format(time.time()-_sttime))        
        ## return ex, sj1
        # self.asm.me = ex
        # self.asm.sj = sj1

class CALCCOV(SUBASE):
    """Calculate coverages.

    Args:
        ae: exons dataframe

    SideEffects:
        add column (cov) to ae dataframe

    TempFiles:
        * cov.txt.gz
        * cov.bed.gz
        * assemble.cov.txt.gz
        * assemble.exons0.covci.txt.gz
        * assemble.exons0.ci.txt.gz
    """

    def call(self):
        ex = self.asm.ae
        fn = self.fnobj
        pr = self.params
        exons = CC.calc_cov_ovl_mp(
                    srcname=ex, #fn.txtname('assemble.exons0'), # cnt 2
                    bwname=fn.bwfile, 
                    dstname=fn.txtname('assemble.exons0.cov'), 
                    np=pr['np'],
                    covciname=fn.txtname('assemble.exons0.covci'), # cnt 1
                    ciname=fn.txtname('assemble.exons0.ci'), # cnt 1
                    colname='cov'
                    )
        self.asm.ae = exons
        fn.write_txt(exons, 'cov')
        # write BED format as well for checking purpose
        exons['sc1'] = exons['cov']
        fn.write_bed(exons, 'cov', ncols=6)

class SETINFO(SUBASE):
    """Set acceptor/donor id, exon categories.

    Args:
        sj: junction dataframe
        ae: exon dataframe

    SideEffects:
        add columns to sj, ae dataframe

    """

    def call(self):
        sj = self.asm.sj
        ae = self.asm.ae
        UT.set_ad_info(sj, ae)
        UT.set_exon_category(sj, ae)

class FINDGENES(SUBASE):
    """Find genes (connected components).

    Args:
        sj: junction dataframe
        ae: exon (all exon = single+multi exons) dataframe

    Related Parameters:
        * np: number of CPU to use, default {np}

    TempFiles:
        * genegraph-*
        * genes.cache.pic

    """

    def call(self):
        sj = self.asm.sj
        ae = self.asm.ae
        fn = self.fnobj
        pr = self.params
        st = self.stats
        genes = GP.find_genes4(
            sj=sj, 
            ae=ae, 
            filepre=fn.fname('genegraph-'), 
            cachename=fn.fname('genes.cache.pic'), 
            np=pr['np']
        )
        numme = len([v for v in genes if len(v)>1])
        numse = len(genes) - numme
        self.info = '#genes:ME{0}, SE(before selection){1}'.format(numme, numse)
        self.stats['FINDGENES.#me_genes'] = numme
        self.stats['FINDGENES.#se_genes'] = numse
        #return genes
        self.asm.genes = genes

class WRITESJEX(SUBASE):
    """Write sj, ex and ci files.

    Args:
        sj: junction dataframe
        ae: exon dataframe

    Related Parameters:
        * findsecovth (bool): whether secovth is found from data

    Files:
        * sj.txt.gz
        * ex.txt.gz
        * ci.txt.gz
        * sj.bed.gz
        * ex.bed.gz

    """

    def call(self):
        sj = self.asm.sj
        ae = self.asm.ae

        fn = self.fnobj
        pr = self.params

        fn.write_txt(sj, 'sj', category='output')
        fn.write_bed(sj, 'sj', category='output', ncols=6)
        fn.write_txt(ae, 'ex', category='output')
        fn.write_bed(ae, 'ex', category='output', ncols=6)
        cipath = fn.txtname('ci', category='output')
        ci = UT.chopintervals(ae, cipath)

        # if pr['findsecovth']: # save exname3
        #     fn.write_txt(ae, 'ex2', category='output')
        #     fn.write_bed(ae, 'ex2', category='output', ncols=6)
        # else: # save in exname2
        #     fn.write_txt(ae, 'ex', category='output')
        #     fn.write_bed(ae, 'ex', category='output', ncols=6)

class WRITEGENES(SUBASE):
    """Write genes (and isoform) BED (and GTF) files for viewing on browsers.

    Args:
        sj: junction dataframe
        ae: exon dataframe

    Related Parameters:
        * writeiso (bool): whether to write isoforms or not, default {writeiso}
        * maxisonum (int):; number of maximum isoforms to write, default {maxisonum}
        * writegtf (bool): whether to write GTF file, default {writegtf}

    Files:
        * genes.bed.gz
        * genes.iso*.bed.gz
        * genes.gtf.gz
        * genes.iso*.gtf.gz

    """

    def call(self):
        sj = self.asm.sj
        ae = self.asm.ae
        fn = self.fnobj
        pr = self.params

        mg = GP.MEGraph4(sj,ae,fn.fname('genegraph-'))
        self.bedwriter = bw = CV.BEDWriter(mg)
        LOG.info(' writing bed12 genes...')
        bw.write(fn.bedname('genes', category='output'))
        maxiso = pr['maxisonum']
        if pr['writeiso']:
            LOG.info(' writing bed12 isoforms...')
            fname = fn.fname('genes.iso{0}.bed.gz'.format(maxiso), category='output')
            bw.write_iso(fname, maxiso)
        if pr['writegtf']:
            self.gtfwriter = gw = CV.GTFWriter(mg)
            LOG.info(' writing gtf genes...')
            gw.write(fn.fname('genes.gtf.gz', category='output'))
            if pr['writeiso']:
                LOG.info(' writing gtf isoforms...')
                fname = fn.fname('genes.iso{0}.gtf.gz'.format(maxiso), category='output')
                gw.write_iso(fname, maxiso)

class CONSISTENTSJ(SUBASE):
    """Remove junctions without connections to exons

    Args:
        sj: junction dataframe
        ae: exon dataframe

    Returns:
        :sj: modify assembler sj

    """
    def call(self):
        sj = self.asm.sj
        ae = self.asm.ae
        dids = set(ae['d_id'].values)
        aids = set(ae['a_id'].values)
        idx = sj['a_id'].isin(aids) & sj['d_id'].isin(dids)
        sj2 = sj[idx]
        self.info = ' #sj {0} => {1}'.format(len(sj), len(sj2))
        self.stats['CONSISTENTSJ.#sj'] = len(sj2)
        self.asm.sj = sj2



# EdgeDetector ##########################################################

class EdgeDetector(object):
    """
    """
    # TODO
    # double mimath: lower always just cut, between lower and higher: also return entire interval

    def __init__(self, bwname, sigmath=3, minth=0.5, delta=15, covratio=0.1, 
                       winsize=15, covth=0.005, smwinsize=151, minintsize=10,
                       mimath=0.15, triggerth=2, mimath2=None, mincutlen=50,
                       aggregateratio=0.1, verbose=False, gap=0, gapth=0):
        """
        Args:
            bwname (str): path to coverage bigwig file 
            sigmath (float): sigma threshold (default 3)
            minth (float):  (default 0.5)
            delta (int): (default 15)
            covratio (float): (default 0.1)
            winsize (int): (default 15)
            covth (float): (default 0.005)
            smwinsize (int): (default 151)
            minintsize (int): (default 10)
            mimath (float): min, max threshold (default 0.15)
            triggerth (float): (default 2)
            mimath2 (float): (default None)
            mincutlen (int): (default 50)
            aggregateratio (float): (default 0.1)
            verbose (bool): (default False)
            gap (int): (default 0)
            gapth (float): (default 0.)

        """
        self.bwname = bwname
        self.sigmath = sigmath
        self.minth = minth
        self.delta = delta
        self.covratio=covratio
        self.covth=covth
        if winsize % 2 ==0:
            winsize += 1 # make sure to use odd number
        self.winsize = winsize
        if smwinsize % 2 ==0:
            smwinsize += 1
        self.smwinsize = smwinsize
        self.minintsize = minintsize # ignore interval smaller than this
        self.mimath = mimath # min/max ratio threshold
        self.triggerth = triggerth # once under mimath, trigger rise if cov > triggerthth*mimath 
        self.swin = N.ones(smwinsize)
        self.win = N.ones(winsize)
        if mimath2 is None:
            self.mimath2 = min(0.5, mimath*2.5)
        else:
            self.mimath2 = mimath2
        self.mincutlen = mincutlen
        self.aggregateratio = aggregateratio
        self.verbose = verbose
        self.gap = gap
        self.gapth = gapth


    def __enter__(self):
        self.fobj = fobj = open(self.bwname, 'rb')
        self.bw = BW.BigWigFile(fobj)
        return self
        
    def __exit__(self, exc_type, exc_value, traceback):
        self.fobj.close()
        
    def get(self, chrom, st, ed):
        # bw.get(chr,st,ed) => [ (st,ed,val), ... ]
        #d = self.delta
        #st0 = st-d
        st0 = st
        #ed0 = ed+d
        ed0 = ed
        vals = N.zeros(ed0-st0,dtype=float)
        for x0,x1,v in self.bw.get(chrom,st0,ed0):
            i0 = x0-st0
            i1 = x1-st0
            vals[i0:i1] = v
        return vals
        
    def findpeak(self, dtmp, th, sign):
        if sign ==-1: # negative peak
            idx = N.nonzero(dtmp<-th)[0]
        else:
            idx = N.nonzero(dtmp>th)[0]
        # if non found at this stage return []
        #LOG.debug 'len(idx)=',len(idx)
        if len(idx)==0:
            return []
        # group into contiguous indices
        groups = []
        x0 = idx[0]
        cgr = [x0]
        for x in idx[1:]:
            if x==(x0+1): # next index
                cgr.append(x)
                x0=x
            else:
                groups.append(N.array(cgr))
                x0 = x
                cgr = [x0]
        groups.append(cgr) # last group
        #LOG.debug 'groups=',groups
        # find min(sign -1) or max(sign +1) within each group
        if sign == -1:
            rslt = [g[N.argmin(dtmp[g])] for g in groups]
        else:
            rslt = [g[N.argmax(dtmp[g])] for g in groups]
        return rslt
        
    def cutQ(self, chrom, st, ed, v=None):
        if v is None:
            v = self.get(chrom, st, ed)
        mima= N.min(v)/N.max(v)
        if mima<=self.mimath:
            return True
        return False

    def plotcut(self,chrom,st,ed,direction='+',ylim=None,xlim=None):
        self.verbose = True
        with self:
            itvs = self.fix(chrom,st,ed,direction)
            v = self.get(chrom,st,ed)
        sm = self.getsm(v)
        dm = self.getdm(v)
        lv = N.log2(v+1)
        ma = N.max(v)
        th1 = self.covth # absolute min
        #th2 = 2**(N.max(N.log2(v+1))*self.covratio)-1 # relative conservative
        th2 = ma*self.covratio
        th3 = ma*self.mimath # last resort
        th3b = th3*self.triggerth
        thd = max(self.minth, self.sigmath*dm.std())
        #thd = 2**thd - 1
        LOG.debug('plotcut variables:ma({ma}),th1({th1}),th2({th2}),th3({th3}),th3b({th3b})'.format(locals()))
        fig,axr = P.subplots(1,1,figsize=(15,3))
        x = N.arange(st,ed)
        P.plot(x,v)
        #P.plot(x,2**N.abs(dm)-1)
        P.plot(x,dm)

        P.plot([x[0],x[-1]], [thd,thd], 'r')
        P.plot([x[0],x[-1]], [th1,th1],'r--')
        P.plot([x[0],x[-1]], [th2,th2],'g--')
        P.plot([x[0],x[-1]], [th3,th3],'m--')
        P.plot([x[0],x[-1]], [th3b,th3b],'m--')

        if ylim:
            axr.set_ylim(ylim)
        else:
            ylim = P.ylim()
        if xlim:
            axr.set_xlim(xlim)
        if len(itvs)>0:
            if direction=='+':
                itvs = sorted(itvs)
                chrom0,st0,ed0 = itvs[0]
                for chrom1, st1, ed1 in itvs:
                    ist = st0-st
                    ied = ed1-st
                    cov = N.mean(v[ist:ied])
                    st0 = ed1
                    #LOG.debug((chrom1, st1, ed1), cov, (ist,ied))
                    P.plot([ed1,ed1],ylim,'r--')
                    P.text(ed1+1,ylim[1]*0.9,'{0:.4f}'.format(cov))
            else:
                vi = v[::-1]
                itvsi = sorted(itvs, key=lambda x:x[2])
                chrom0,st0,ed0 = itvsi[0]
                for chrom1, st1, ed1 in itvsi:
                    ist = ed - ed0
                    ied = ed - st1
                    cov = N.mean(vi[ist:ied])
                    ed0 = st1
                    #LOG.debug((chrom1, st1, ed1), cov, (ist,ied))
                    P.plot([st1,st1],ylim,'r--')
                    P.text(st1+1,ylim[1]*0.9,'{0:.4f}'.format(cov))

        return fig

    def cutboth(self, chrom, st, ed):
        if (ed-st)<self.mincutlen:
            return [],[],[(chrom,st,ed)]
        v = self.get(chrom,st,ed)
        mima = N.min(v)/N.max(v)
        if mima>=self.mimath2:
            return [], [], [(chrom,st,ed)]
        sws = self.smwinsize # smooth window for abs cov th detection
        swin = self.swin
        v0 = N.concatenate([swin*v[0], v, swin*v[-1]])
        sm = N.convolve(v0, swin/float(sws), 'same')[sws:-sws]
        # left and right
        cleft = self.fix(chrom,st,ed,'+',v,sm)
        cright = self.fix(chrom,st,ed,'-',v,sm)
        if (mima>self.mimath) or ((len(cleft)==0) and (len(cright)==0)):
            return cleft, cright, [(chrom,st,ed)]
        return cleft, cright, []

    def cutone(self, chrom, st, ed, direction):
        if (ed-st)<self.mincutlen:
            return [],[],[(chrom,st,ed)]
        v = self.get(chrom,st,ed)
        mima = N.min(v)/N.max(v)
        if mima>=self.mimath2:
            return [], [], [(chrom,st,ed)]
        sws = self.smwinsize # smooth window for abs cov th detection
        swin = self.swin
        v0 = N.concatenate([swin*v[0], v, swin*v[-1]])
        sm = N.convolve(v0, swin/float(sws), 'same')[sws:-sws]
        # left and right
        if direction=='+':
            cleft = self.fix(chrom,st,ed,'+',v,sm)
            if (mima>self.mimath) or (len(cleft)==0):
                return cleft, [], [(chrom,st,ed)]
            return cleft,[],[]
        else:
            cright = self.fix(chrom,st,ed,'-',v,sm)
            if (mima>self.mimath) or (len(cright)==0):
                return [], cright, [(chrom,st,ed)]
            return [], cright, []

    def getdm(self, v):
        ws = self.winsize # smooth window for derivative
        hws = int(ws/2) # ws has to be odd number
        win = self.win 
        #lv = N.log2(v+1)
        lv = v
        lv1 = N.concatenate([win*lv[0], lv, win*lv[-1]])
        sm = N.convolve(lv1, win/float(ws), 'same')
        dm = (sm[ws:]-sm[:-ws])[(hws+1):-hws]
        return dm

    def getsm(self, v):
        sws = self.smwinsize # smooth window for abs cov th detection
        swin = self.swin
        v0 = N.concatenate([swin*v[0], v, swin*v[-1]])
        return N.convolve(v0, swin/float(sws), 'same')[sws:-sws]

    def _findmax(self,v):
        gap = self.gap
        minidx = N.nonzero(v<=self.gapth)[0]
        if len(minidx)==0:
            return len(v)
        cidx = minidx[0]
        size = 0
        for idx in minidx[1:]:
            if idx-cidx==1: # contiguous
                size +=1
                if size>gap:
                    return min(len(v),cidx+1)
            else:
                size = 0
            cidx = idx
        # all gaps were < gap size th
        return len(v)

    def fix(self, chrom, st, ed, direction, v=None, sm=None):
        if v is None:
            v = self.get(chrom,st,ed)
        if sm is None:
            sm = self.getsm(v)
        if direction=='-':
            v = v[::-1]
            sm = sm[::-1]
        # gap?
        olen = len(v)
        maidx = self._findmax(v)
        v = v[:maidx]
        sm = sm[:maidx]
        if self.verbose:
            LOG.debug('size={0}=>{1}'.format(olen,maidx))

        mima = N.min(v)/N.max(v)
        l0,l1,l2,l3,l4 = 0,0,0,0,0
        eds0 = self.detect_sharp_drops(v, olen)
        l0 = len(eds0)
        eds1 = self.detect_rise(v)
        l1 = len(eds1)
        if l1==0:
            eds1 = self.detect_low_level(v)
            l2 = len(eds1)
            if l2==0 and l0==0:
                eds1 = self.detect_min(v)
                l3 = len(eds1)
                if l3==0:# nothing detected, last resort
                    eds1 = self.detect_rise2(v)
                    l4 = len(eds1)
        if self.verbose:
            LOG.debug('drop({0}), rise({1}), low({2}), min({3}), rise2({4})'.format(l0,l1,l2,l3,l4))
        if len(eds1)>0:
            eds1 = self.trim(v,sm,eds1)
        eds = eds0+eds1  
        if len(eds)==0 and (olen>maidx):
            eds = [maidx]
        eds = self._aggregate(olen,eds)
        if len(eds)==0:
            return []
        if direction=='-':
            return [(chrom, ed-x, ed) for x in eds] # (st=ed-idx, ed)
        return [(chrom, st, st+x) for x in eds] # (st, st+idx)

    def _aggregate(self, maxlen, bs):
        if len(bs)==0:
            return []
        bs = sorted(bs)
        mis = self.minintsize
        while(len(bs)>0 and bs[0]<mis):
            bs = bs[1:]
        while(len(bs)>0 and maxlen-bs[-1]<mis):
            bs = bs[:-1]
        if len(bs)>1:
            bsi = bs[::-1]
            ar = self.aggregateratio # 10% of longer
            cur = bsi[0]
            mis2 = ar*cur
            bs2 = [cur]
            for b in bsi[1:]:
                if (cur-b)>mis2:#ar*cur:
                    bs2.append(b)
                    cur = b
            bs = sorted(bs2)
        # if len(bs)>1: # aggregate close by
        #     bs = [bs[0]] + [x1 for x0,x1 in zip(bs[:-1], bs[1:]) if x1-x0>mis]
        return bs

    def detect_rise(self, v):
        ma = N.max(v)
        mima = N.min(v)/ma
        if mima<self.mimath:
            th1 = ma*self.mimath
            th2 = th1*self.triggerth
        elif mima<self.mimath2:
            th1 = ma*self.mimath2
            th2 = th1+(ma*self.mimath)
        else:
            return []
        ist = N.nonzero(v<th1)[0][0]
        mis = self.minintsize
        if (ist<mis) or (len(v)-ist<mis):
            return []
        #LOG.debug('detect_rise startidx:{0}'.format(ist))
        cmin = v[ist]
        imin = ist
        for i in range(ist+1,len(v)):
            if cmin>v[i]:
                cmin = v[i]
                imin = i
            elif v[i]>th2: # rise detected
                #LOG.debug(' rise at:{0},trigger at:{1}'.format(imin, i))
                return [imin]
        return []

    def detect_rise2(self, v):
        ma = N.max(v)
        mi = N.min(v)
        mima = mi/ma
        if mima<self.mimath:
            th1 = ma*self.mimath
            th2 = th1*self.triggerth
        elif mima<self.mimath2:
            th1 = ma*self.mimath2
            th2 = th1+(ma*self.mimath)
        else:
            return []
        ist = N.nonzero(v>=th1)[0][0]
        mis = self.minintsize
        if (ist<mis) or (len(v)-ist<mis):
            return []
        #LOG.debug('detect_rise startidx:{0}'.format(ist))
        cmin = v[ist]
        imin = ist
        for i in range(ist+1,len(v)):
            if cmin>v[i]:
                cmin = v[i]
                imin = i
            elif v[i]>th2: # rise detected
                #LOG.debug(' rise at:{0},trigger at:{1}'.format(imin, i))
                return [imin]
        return []

    def detect_low_level(self, v):
        ma = N.max(v)
        l = len(v)
        mis = self.minintsize
        th1 = self.covth # absolute min
        th2 = ma*self.covratio
        th3 = 2*th2
        idx = N.nonzero(v<th1)[0]
        def _chk(idx):
           return (len(idx)==0) or (idx[0]<mis) or (l-idx[0]<mis) 
        if _chk(idx):
            idx = N.nonzero(v<th2)[0]
            if _chk(idx):
                idx = N.nonzero(v<th3)[0]
                if _chk(idx):
                    return []
        return [idx[0]]

    def detect_min(self, v):
        idx = N.argmin(v)
        l = len(v)
        mis = self.minintsize
        if (idx<mis) or (l-idx<mis) :
            return []
        return [idx]


    def trim(self, v, sm, eds):
        # just trim left most
        th1 = self.covth # absolute min
        #th2 = 2**(N.max(N.log2(v+1))*self.covratio)-1 # relative conservative
        eidx=eds[-1]
        if eidx==0:
            return eds
        ma = N.max(v[:eidx]) # max of the region
        th2 = ma*self.covratio
        #acovth = max(th1,th2)
        #acovth=th2
        #self._acovth = acovth
        while(eidx>0 and sm[eidx]<th2):
            eidx = eidx-1
        while(eidx>0 and v[eidx]<th2):
            eidx = eidx-1
        return [x for x in eds[:-1] if x<eidx]+[eidx]

    def detect_sharp_drops(self, v, maxlen):
        ws = self.winsize # smooth window for derivative
        hws = int(ws/2) # ws has to be odd number
        win = self.win 
        #lv = N.log2(v+1)
        lv = v
        lv1 = N.concatenate([win*lv[0], lv, win*lv[-1]])
        sm = N.convolve(lv1, win/float(ws), 'same')
        dm = (sm[ws:]-sm[:-ws])[(hws+1):-hws]
        sigma = dm.std()
        th = max(self.minth, self.sigmath*sigma)
        bs = self.findpeak(dm, th, -1) # detect drop
        mis = self.minintsize
        return self._aggregate(maxlen,bs)


    def calc_stats(self, chrom, st, ed, v=None):
        if v is None:
            v = self.get(chrom,st,ed) # raw cov between (st,ed)
        th = self.covth # abs th
        #d = self.delta # not use 
        sws = self.smwinsize # smooth window for abs cov th detection

        ws = self.winsize # smooth window for derivative
        hws = int(ws/2) # ws has to be odd number
        swin = self.swin
        win = self.win 

        lv = N.log2(v+1)
        
        v0 = N.concatenate([swin*v[0], v, swin*v[-1]])
        ssm = N.convolve(v0, swin/float(sws), 'same')[sws:-sws]

        lv1 = N.concatenate([win*lv[0], lv, win*lv[-1]])
        sm = N.convolve(lv1, win/float(ws), 'same')
        dm = (sm[ws:]-sm[:-ws])[(hws+1):-hws]
        #LOG.debug len(dm),len(lv),len(sm)
        mi = N.min(v)
        ma = N.max(v)
        mima = mi/ma
        #dmax = N.max(N.abs(dm))
        #rdmax = dmax/N.max(sm)
        sigma = dm.std()
        th = max(self.minth, self.sigmath*sigma)
        return v, dm, th, ssm, mima

        #ltmp = N.log2(tmp+1)
        #w = N.ones(sws)/float(sws)
        #sm = N.convolve(ltmp, w, 'same')
        #dtmp = N.array(ltmp[1:]-ltmp[:-1])
        #hws = int(ws/2)
        #dtmp = N.array([0]+[(ltmp[i:i+ws].mean()-ltmp[max(0,i-ws):i].mean()) for i in range(1,len(ltmp))])    
        #w2 = N.ones(ws)/float(ws)
        #sm2 = N.convolve(ltmp, w2, 'same')
        #dtmp1 = sm2[ws:]-sm2[:-ws]
        #dtmp = N.concatenate([N.zeros(hws+1),dtmp1, N.zeros(hws)])
        #sigma = dtmp.std()
        #th = max(self.minth, self.sigmath*sigma)
        #return tmp, dtmp, th, sm
        
    def find(self, chrom, st, ed, find='both',returnall=True, ocov=0., v=None):
        #tmp, dtmp, th, sm = self._calc_dtmp(chrom,st,ed)
        tmp, dtmp, th, sm, mima = self.calc_stats(chrom,st,ed,v)
        # trim edge
        #d = self.delta
        #s0 = st-d 
        s0 = st
        # acovth = self.covth
        # use adaptive abs cov
        acovth = max(self.covth, 2**(N.max(N.log2(tmp+1))*self.covratio)-1)
        self._acovth = acovth
        if find=='both' or find=='drop': # trim right side
            eidx=len(tmp)-1
            while(eidx>0 and sm[eidx]<acovth):
                eidx = eidx-1
            while(eidx>0 and tmp[eidx]<acovth):
                eidx = eidx-1
            e1 = min(ed,s0+eidx)
        else:
            e1 = ed
        if find=='both' or find!='drop': # trim left side
            sidx=0
            n=len(tmp)
            while(sidx<n-1 and sm[sidx]<acovth):
                sidx+=1
            while(sidx<n-1 and tmp[sidx]<acovth):
                sidx+=1
            s1 = max(st, s0+sidx)
        else:
            s1 = st
        self._st = st
        self._ed = ed
        self._s1 = s1
        self._e1 = e1

        if (e1<=st) or (s1>=ed): # entire region below covth return original
            if returnall:
                return [((st,ed),ocov,False)]
            return [((st,ed),ocov)]

        if find=='both':
            idx0 = self.findpeak(dtmp, th, -1) #N.nonzero(dtmp<-th)[0]+1
            idx1 = self.findpeak(dtmp, th, +1) #N.nonzero(dtmp>th)[0]
        elif find=='drop':
            idx0 = self.findpeak(dtmp, th, -1) #N.nonzero(dtmp<-th)[0]+1
            idx1 = []
        else:
            idx1 = self.findpeak(dtmp, th, +1) #N.nonzero(dtmp>th)[0]
            idx0 = []
        idx = sorted(set(idx0).union(set(idx1)))
        
        bs = [s0+x for x in idx if ((s0+x)>=s1) and ((s0+x)<=e1)] # boundaries
        mis = self.minintsize
        while(len(bs)>0 and abs(bs[0]-s1)<mis):
            bs = bs[1:]
        while(len(bs)>0 and abs(bs[-1]-e1)<mis):
            bs = bs[:-1]
        if len(bs)>1: # aggregate close by
            bs = [bs[0]] + [x1 for x0,x1 in zip(bs[:-1], bs[1:]) if x1-x0>mis]
        # for each ed1 in bs make (chr,st,ed1)
        # bs=[b1,b2,b3,...,bn]
        # intervals = (st,b1),(b1,b2),(b2,b3),...,(bn-1,bn),(bn,ed)
        if len(bs)>0:
            itvs =[(s1,bs[0])]+[x for x in zip(bs[:-1],bs[1:])]+[(bs[-1],e1)]
        else:
            itvs = [(s1,e1)]
        # for each interval, calculate cov (just sum of tmp/len)
        covs = N.array([N.mean(tmp[x-s0:y-s0]) for x,y in itvs])
        #covs = [N.sum(tmp[x-s0:y-s0])/(y-x) for x,y in itvs]
        #rth = self.covratio
        #th = self.covth
        acovth2 = 2**(N.max(N.log2(covs+1))*self.covratio)-1
        self._acovth2 = acovth2
        if len(covs)>1:
            #i1 = set([i for i in range(len(covs)-1) if ((covs[i]/covs[i+1])>rth)&(covs[i]>th)])
            #i2 = set([i for i in range(1,len(covs)) if ((covs[i]/covs[i-1])>rth)&(covs[i]>th)])
            #idx = sorted(i1.union(i2))
            idx = [i for i in range(len(covs)) if covs[i]>acovth2]
        else:
            idx = [0]
        if returnall:
            return [(itvs[x], covs[x], x in idx) for x in range(len(itvs))]
        if len(idx)<3:
            return [(itvs[x], covs[x]) for x in idx]
        return [(itvs[0], covs[0]), (itvs[idx[-1]], covs[idx[-1]])] # only return first & last
        
    def plotone(self,chrom,st,ed,strand='+',utr='3', find=None, ylim=None, xlim=None):
        if find is None:
            if utr=='3':
                find = 'drop' if strand=='+' else 'rise'
            else:
                find = 'drop' if strand=='-' else 'rise'
        with self:
            #tmp, dtmp, th, sm = self._calc_dtmp(chrom,st,ed)
            tmp, dtmp, th, sm, mima = self.calc_stats(chrom,st,ed)
            itcv = self.find(chrom,st,ed,find)
        ltmp = N.log2(tmp+1)
        #covs = N.log2(N.array([x[1] for x in itcv])+1)
        #lacovth = N.max(ltmp)*self.covratio
        #lacovth2 = N.max(lcovs)*self.covratio
        self.lacovth = lacovth = N.log2(self._acovth+1)
        self.lacovth2 = lacovth2 = N.log2(self._acovth2+1)

        #         cs0 = N.cumsum(tmp)
        #         tot = N.sum(tmp)
        #         cs1 = tot - cs0
        #         cs0 = cs0/tot
        #         cs1 = cs1/tot
        
        fig,axr = P.subplots(1,1,figsize=(15,3))
        #d = self.delta
        #x = N.arange(st-d,ed+d)
        x = N.arange(st,ed)
        P.plot(x,ltmp)
        P.plot(x,N.abs(dtmp))
        #         P.plot(x,cs0)
        #         P.plot(x,cs1)
        
        P.plot([x[0],x[-1]], [th,th], 'r')
        P.plot([x[0],x[-1]], [lacovth,lacovth],'r--')
        P.plot([x[0],x[-1]], [lacovth2,lacovth2],'g--')
        if ylim:
            axr.set_ylim(ylim)
        else:
            ylim = P.ylim()
        if xlim:
            axr.set_xlim(xlim)
        P.text(x[0], ylim[1]*0.8, '{0}:{1}-{2}:{3}'.format(chrom,st,ed,strand))
        for (x1,x2),cov, sel in itcv:
            LOG.debug( (x1,x2), cov, sel )
            if sel:
                P.plot([x1,x1],ylim,'r--')
                P.plot([x2,x2],ylim,'r--')
                P.text(x1+1,ylim[1]*0.9,'{0:.4f}'.format(cov))
            else:
                P.plot([x1,x1],[0,ylim[1]*0.5],'c--')
                P.plot([x2,x2],[0,ylim[1]*0.5],'c--')                
                P.text(x1+1,ylim[1]*0.45,'{0:.4f}'.format(cov))
        return fig

try: # only Python3
    for klass in SUBASE.__subclasses__():
        klass.__doc__ = klass.__doc__.format(**PARAMSDOC)
except:
    pass
