"""

.. module:: merge
    :synopsis: module for merging multiple assemblies

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""

import subprocess
import os
import gzip
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)
import shutil

import pandas as PD
import numpy as N

from jgem import utils as UT
from jgem import gtfgffbed as GGB
from jgem import filenames as FN
from jgem import assembler as AS
from jgem import bedtools as BT
from jgem import bigwig as BW
from jgem import calccov as CC

MERGECOVPARAM = dict(
    np = 1, # number of CPU to use
    genome = 'mm10', # UCSC genome name
    covth = 0, # exon read digitization threshold
    covdelta = 1, # exon read digitization unit
    uth=0, # unique count threshold
    mth=5, # non-unique count threshold
    th_ratio=1e-3, # discard junctions less than this portion within overlapping junctions
    th_detected=2, # at least observed in 3 samples
    th_maxcnt1=0.1, # max read count should be larger than this
    th_maxcnt2=1, # if max read count is larger than this then single sample observation is OK
)
MERGEASMPARAM = dict(
    se_maxth=0.5,   # SE maxcov threshold
    se_gidstart=50000, # SE gidx start
)


class MergeInputNames(FN.FileNamesBase):
    """Filelname manager for generating inputs for merging process.

    Attributes:
        sampleinfo: sample info dataframe (with columns: name, sjpath, expath, bwfile, sjfile)
        code: merge identifier (for merge input generation part)
        outdir: output directory

    All outputs and temporary files are prefixed by **outdir/code**

    SampleInfo Columns:
        name: sample name (unique)
        sjexpre: path prefix to SJ/EX files (sj.txt.gz, ex.txt.gz will be added)
        bw_path: original bigwig coverage
        sjbed_path: original junction bed file (converted from SJ.out.tab)

    """

    def __init__(self, sampleinfo, code, outdir):
        self.si = sampleinfo
        self.code = code
        self.outdir = outdir

        prefix = os.path.join(outdir, code)
        super(MergeInputNames, self).__init__(prefix)

    def expaths(self):
        return [(n, '{0}.ex.txt.gz'.format(x)) for n,x in self.si[['name','sjexpre']].values]

    def sjopaths(self):
        return [(n, '{0}.sj.txt.gz'.format(x)) for n,x in self.si[['name','sjexpre']].values]

    def sjpaths(self):
        return self.si[['name','sjbed_path']].values # BED

    def bwpaths(self):
        return self.si[['name','bw_path']].values # BIGWIG

    def sj0_bed(self):
        return self.bedname('sj0', category='output')

    def sj_bed(self, strand):
        """SJ output, strand = p,n """
        return self.bedname('sj.{0}'.format(strand), category='output')

    def sj2_txt(self):
        """intermediate output with selection parameters """
        return self.txtname('sj2', category='output')

    def allsj_txt(self):
        return self.txtname('allsj', category='output')

    def allsj_stats(self):
        return self.txtname('allsj.stats', category='output')

    def ex_bw(self, k):
        """BW output, strand = p,n """
        return self.fname('ex.{0}.bw'.format(k), category='output')

    def ex_bed(self, k):
        return self.bedname('ex.{0}'.format(k))

    def agg_bw(self):
        return self.fname('allsample.bw', category='output')

    def snames(self):
        return list(self.si['name'])

class MergeAssemblyNames(FN.FileNamesBase):
    """Filelname manager for the assembling part of the merging process.

    Attributes:
        code: merge identifier (for merge input generation part)
        outdir: output directory
        refgtf: reference gtf path if using it for finding SE cov threshold

    All outputs and temporary files are prefixed by **outdir/code**

    """

    def __init__(self, code, outdir, refgtf='.gtf'):
        self.code = code
        self.outdir = outdir
        self.refgtf = refgtf

        prefix = os.path.join(outdir, code)
        super(MergeAssemblyNames, self).__init__(prefix)

    def ex_out(self, ext='txt'):
        return self.fname('ex.{0}.gz'.format(ext), category='output')

    def sj_out(self, ext='txt'):
        return self.fname('sj.{0}.gz'.format(ext), category='output')

    def ci_out(self):
        return self.fname('ci.txt.gz', category='output')

    def genes_out(self, ext='txt'):
        return self.fname('genes.{0}.gz'.format(ext), category='output')


class MergeInputs(object):
    """Generates inputs to merge assembling. Creates multiple bigwig files and a junction file.
    Multiple bigwig file are for exons in each strand and for single exons. 

    1. pepare 3 merged bigwigs (ME+/ME-/SE) from all assemblies
    2. prepare an aggregated junction file
    3. select junctions
    4. calculate average bigwig from original bigwigs

    (Note: at 1. bigwigs are from exon model outputs from each assembly, at 4. bigwig is just an
    average of all original bigiwigs which were inputs to each assembly.)
    """
    def __init__(self, fnobj, genome='mm10', **kw):
        """
        Args:
            fnobj: MergeInputNames object

        Other keywords arguments:
            * genome: UCSC genome name
            * np: number of CPU to use
            * covth: only use exons with coverage > covth (default 0)
            * covdelta: exon coverage quantization delta
            * uth: unique reads threshold for junctions (default 0)
            * mth: non-unique reads threshold for junctions (default 5)
            * th_ratio: threshold for selecting junctions by overlapping ratio 
              (default 0.001, i.e. if a junction's read is less than 1/1000 of the sum of 
              all the reads of overlapping junctions then the junction is discarded)
            * th_detected: junctions need to be detected in more than this number of samples (default 2)
              unless max junction reads across samples is larger than th_maxcnt2
            * th_maxcnt1: max junction reads across sample has to be larger than this (default 0.1)
            * th_maxcnt2: if max junction reads across samples is larger than this, ignore th_detected
              (default 4)

        """
        self.fnobj = fnobj
        self.params = MERGECOVPARAM.copy()
        self.params.update(kw)
        self.params['genome'] = genome
        self.chromsizes = UT.chromsizes(genome)
        self.chroms = UT.chroms(genome)
        self.tgts =  ['mep','men','se'] #'sep','sen']

    def prepare(self):
        """ Prepare merged bigwig coverage files, merged junction files and aggregated bigwig file (average cov)."""
        self.make_ex_bigwigs()
        self.make_sj_bed()
        self.aggregate_bigwigs()
        fn.delete(delete=['temp'],protect=['output'])

    def make_ex_bigwigs(self):
        """Make 5 bigwig files from assembly exon outputs by treating each exon as reads weighted by coverage. """
        fn = self.fnobj
        # first make BED files
        LOG.debug('making BEDS...')
        self._make_ex_beds()
        # then convert them to BIGWIGs
        for k in self.tgts:
            bedpath = fn.ex_bed(k) #fn.fname('ex.{0}.bed.gz'.format(k))
            bwpath = fn.ex_bw(k) #fn.fname('ex.{0}.bw'.format(k), category='output')
            LOG.debug('converting {0} to BIGWIG...'.format(bedpath))
            #totbp,covbp = BT.get_total_bp_bedfile(bedpath, bed12=False)
            #scale = float(covbp)/totbp # normalize average coverage to 1.
            #LOG.info('{0}:totbp={1},covbp={2},scale={3}'.format(bedpath,totbp,covbp,scale))
            # scale = 1e8/totbp # old normalization = 1e6/totaligned when readlen=100bp
            # ^==== TODO: Is normalizing to totaligned good? 
            # If complexity (#genes) is bigger then per element cov is smaller. 
            # This will affect the average cov level and noise level. 
            # ====> make average coverage constant

            # [Q] normalize chrom-wise? Mouse chr11,19,X seems higher than average
            #     in addition to the obvious low expressing chrY

            # 2016-04-28: don't scale just reflect read depth all the way through

            #BT.bed2bw(bedpath, self.chromsizes, bwpath, scale=scale)
            BT.bed2bw(bedpath, self.chromsizes, bwpath, scale=None)

        # delete temp files
        fn.delete(delete=['temp'],protect=['output'])

    def _make_ex_beds(self):
        fn = self.fnobj
        pr = self.params
        np = pr['np']
        chroms = self.chroms
        expaths = fn.expaths() # [(name, expath), ...]
        dcode = 'ex.'
        dstpre = fn.fname(dcode)
        tgts = self.tgts # ['mep','men','se'] #sep','sen']
        th = pr['covth']
        delta = pr['covdelta']
        UT.makedirs(os.path.dirname(dstpre))
        args = [(expaths, dstpre, x, th, delta, tgts) for x in chroms]

        rslts = UT.process_mp(make_ex_bed_chr, args, np, doreduce=False)

        # concatenate
        LOG.debug('concatenating chroms...')
        UT.makedirs(os.path.dirname(fn.ex_bed(tgts[0])))
        for k in tgts:
            bf = fn.ex_bed(k) #fn.fname('ex.{0}.bed.gz'.format(k))
            with open(bf,'wb') as dst:
                for x in chroms:
                    sf = fn.fname('{0}{1}{2}.gz'.format(dcode,x,k))
                    with open(sf,'rb') as src:
                        shutil.copyfileobj(src, dst)
                        #dst.write(src.read())
        
    def make_sj_bed(self):
        """Make merged junction input file. """
        self.make_sj0_bed() # aggregate junctions
        self.collect_sj() # collect sample junction counts
        self.select_sj() # select junctions ==> do selection in the assembler? (keep it for now )
        self.write_sj() # write sj.p, sj.n 
        self.fnobj.delete(delete=['temp'],protect=['output'])

    def make_sj0_bed(self):
        """Aggregate all junctions in the samples. Ucnt, mcnt will be the sum over all samples. """
        pr = self.params
        fn = self.fnobj
        uth = pr['uth']
        mth = pr['mth']
        np = pr['np']
        chroms = self.chroms
        sjpaths = fn.sjpaths() # [(name, sjpath),..]
        scode = 'sjbed.gz'
        asjpath = fn.fname(scode) # aggregated sj file
        args = [(sjpaths, asjpath, x, uth, mth) for x in chroms]
        UT.makedirs(os.path.dirname(asjpath))
        rslts = UT.process_mp(make_sj_bed_chr, args, np, doreduce=False)

        # concatenate
        LOG.debug('merge sj: concatenating chroms...')
        UT.makedirs(os.path.dirname(asjpath))
        with open(asjpath,'wb') as dst:
            for x in chroms:
                sf = fn.fname('{0}{1}.gz'.format(scode,x))
                with open(sf,'rb') as src:
                    dst.write(src.read())

        dstpath = fn.sj0_bed() #fn.fname('sj0.bed.gz', category='output') # before selection

        # group same junctions
        msj = UT.read_pandas(asjpath, names=['chr','st','ed','strand','src','ucnt','mcnt'])
        # average <= 2016-04-28 don't average just aggregate
        # scale = 1/float(len(sjpaths)) 
        # msj['ucnt'] = scale*msj['ucnt']
        # msj['mcnt'] = scale*msj['mcnt']
        # unique junctions
        msjg = msj.groupby(['chr','st','ed','strand'])[['ucnt','mcnt']].sum().reset_index()
        # msjg['sc1'] = msjg['ucnt'] # for BED
        # msjg['tst'] = msjg['mcnt'] # for BED
        # u = msjg['ucnt'].map('{:.2}'.format)
        # m = msjg['mcnt'].map('{:.2}'.format)
        u = msjg['ucnt'].astype(str)
        m = msjg['mcnt'].astype(str)
        msjg['_id'] = N.arange(len(msjg))
        msjg['name'] = msjg['_id'].astype(str)+'_u:'+u+'_m:'+m
        cols = GGB.SJCOLS #GGB.BEDCOLS[:7] # chr,st,ed,name,sc1,strand,tst
        UT.write_pandas(msjg[cols], dstpath, '') # BED file
        self.sj0 = msjg

    def collect_sj(self):
        """From aggregated junctions, collect original reads for each samples.

        Inputs:
            Aggregated junction file ('sj0.bed.gz')
            Sample junction files (input to original assemblies ~ SJ.out.tab)

        Outputs:
            'allsj.txt.gz'
        """
        fn = self.fnobj
        if hasattr(self, 'sj0'):
            msj = self.sj0
        else:
            self.sj0 = msj = GGB.read_sj(fn.sj0_bed())
        # sc1: ucnt, tst: mcnt

        #sjpaths = fn.sjopaths() # [(name,sjpath),...], output of assembly (restricted)
        sjpaths = fn.sjpaths() # [(name,sjpath),...], input of assembly (all observed junctions)
        snames = [x[0] for x in sjpaths]

        msj['locus'] = UT.calc_locus_strand(msj)
        for i, (sname,spath) in enumerate(sjpaths):
            #scale = 1e6/float(aligned)
            sj = GGB.read_sj(spath)
            sj['locus'] = UT.calc_locus_strand(sj)
            sj['cnt'] = (sj['ucnt']+sj['mcnt']) #*scale <== sjbed is already normalized
            l2u = dict(UT.izipcols(sj, ['locus','cnt']))
            msj[sname] = [l2u.get(x,0) for x in msj['locus']]  
        msj['#detected'] = (msj[snames]>0).sum(axis=1) # number of samples with reads>0
        msj['maxcnt'] = msj[snames].max(axis=1) # max reads
        self.allsj = msj
        UT.write_pandas(msj, fn.allsj_txt(), 'h')
        UT.write_pandas(msj[['locus','#detected','maxcnt']], fn.allsj_stats(), 'h')

    def select_sj(self):
        """Select aggregated junctions according to several metrics"""
        fn = self.fnobj
        pr = self.params

        # calc self intersection to find overlapping junctions
        a = b = fn.sj0_bed()
        c = fn.fname('sj0.ovl.txt.gz')
        c = BT.bedtoolintersect(a,b,c,wao=True)
        # calc ratio
        cols0 = GGB.SJCOLS
        cols = cols0+['b_'+x for x in cols0]+['ovl']
        sjovl = UT.read_pandas(c, names=cols)

        sjovl = sjovl[sjovl['strand']==sjovl['b_strand']] # same strand
        LOG.debug('len(sjovl)={0}'.format(len(sjovl)))

        sjgr = sjovl.groupby(['chr','st','ed','strand'])
        sj2 = sjgr[['ucnt','mcnt','name']].first()
        sj2['ucnt_sum'] = sjgr['b_ucnt'].sum()
        sj2['mcnt_sum'] = sjgr['b_mcnt'].sum()
        sj2['sum'] = sj2['ucnt_sum']+sj2['mcnt_sum']
        sj2['cnt'] = sj2['ucnt']+sj2['mcnt']
        self.sj2 = sj2 = sj2.reset_index() # need chr,st,ed,strand at next step
        sj2['locus'] = UT.calc_locus_strand(sj2)
        sj2['ratio'] = sj2['ucnt']/sj2['ucnt_sum']
        sj2['ratio_m'] = sj2['mcnt']/sj2['mcnt_sum']
        sj2['ratio_a'] = sj2['cnt']/sj2['sum']

        # add #detected, maxcnt
        if hasattr(self, 'allsj'):
            allsj = self.allsj
        else:
            self.allsj = allsj = UT.read_pandas(fn.allsj_txt())

        l2d = UT.df2dict(allsj, 'locus', '#detected')
        l2m = UT.df2dict(allsj, 'locus', 'maxcnt')

        sj2['#detected'] = [l2d[x] for x in sj2['locus']]
        sj2['maxcnt'] = [l2m[x] for x in sj2['locus']]
        UT.write_pandas(sj2, fn.sj2_txt(), 'h')

        # select 
        idx1 = (sj2['ratio']>=pr['th_ratio'])|(sj2['ratio_a']>=pr['th_ratio'])
        idx2 = sj2['#detected']>pr['th_detected']
        idx3 = sj2['maxcnt']>pr['th_maxcnt1']
        idx4 = sj2['maxcnt']>pr['th_maxcnt2']

        self.sj4 = sj4 = sj2[idx1&((idx2&idx3)|idx4)]
        LOG.info('selectsj: in {0}'.format(len(sj2)))
        LOG.info('selectsj: {0} smaller than th_ratio({1})'.format(N.sum(~idx1),pr['th_ratio']))
        LOG.info('selectsj: {0} smaller than th_detected({1})'.format(N.sum(~idx2),pr['th_detected']))
        LOG.info('selectsj: {0} smaller than th_maxcnt1({1})'.format(N.sum(~idx3),pr['th_maxcnt1']))
        LOG.info('selectsj: {0} larger than th_maxcnt2({1})'.format(N.sum(idx4),pr['th_maxcnt2']))
        LOG.info('#selected SJ:{0}<={1}'.format(len(sj4),len(sj2)))
        cols = GGB.SJCOLS
        sjg = allsj.set_index('locus')[cols]
        self.sj1 = sjg.ix[sj4['locus'].values]
        # # sjgp = sjg.ix[sj4[sj4['strand']=='+']['locus'].values]
        # # sjgn = sjg.ix[sj4[sj4['strand']=='-']['locus'].values]
        # sjgp = sjg.ix[sj4[sj4['strand'].isin(['+','.'])]['locus'].values]
        # sjgn = sjg.ix[sj4[sj4['strand'].isin(['-','.'])]['locus'].values]        
        # UT.write_pandas(sjgp, fn.sj_bed('p'), '')
        # UT.write_pandas(sjgn, fn.sj_bed('n'), '')

    def write_sj(self):
        fn = self.fnobj
        if hasattr(self, 'sj1'):
            sjg = self.sj1
        else:
            sjg = self.sj0
        cols = GGB.SJCOLS
        sjgp = sjg[sjg['strand'].isin(['+','.'])][cols]
        sjgn = sjg[sjg['strand'].isin(['-','.'])][cols]
        UT.write_pandas(sjgp, fn.sj_bed('p'), '')
        UT.write_pandas(sjgn, fn.sj_bed('n'), '')

    def aggregate_bigwigs(self):
        fn = self.fnobj
        pr = self.params
        bwfiles = [x[1] for x in fn.bwpaths()] # [(name,bwfile),...]
        dstpath = fn.agg_bw()
        # scale = 1./len(bwfiles) # average
        # BW.merge_bigwigs_mp(bwfiles, pr['genome'], dstpath, scale=scale, np=pr['np'])
        BW.merge_bigwigs_mp(bwfiles, pr['genome'], dstpath, scale=None, np=pr['np'])

def make_ex_bed_chr(expaths, dstpre, chrom, covth, covdelta, tgts):
    withcov = True
    paths = {k:dstpre+chrom+k for k in tgts}
    dst = {k:open(paths[k],'w') for k in tgts }
    for i,(name, path) in enumerate(expaths):
        ex = UT.read_pandas(path)
        if withcov:
            ex['dup'] = ((ex['cov']-covth)/covdelta).astype(int)+1
        ex = ex[(ex['chr']==chrom)&(ex['cov']>covth)]
        me = ex[ex['cat']!='s']
        se = ex[ex['cat']=='s'].copy()
        comb = {'mep':(me,('+','.')),
                'men':(me,('-','.')),
                'sep':(se,('+','.')),
                'sen':(se,('-','.')),
                'se':(se,['.'])}
        comb = {t:comb[t] for t in tgts}
        for k,(tgt,strand) in comb.items():
            tgt = tgt[tgt['strand'].isin(strand)]
            if withcov:
                def _gen():
                    for chrom,st,ed,dup in UT.izipcols(tgt, ['chr','st','ed','dup']):
                        for i in range(dup):
                            yield '{0}\t{1}\t{2}\n'.format(chrom,st,ed)
                recs = [x for x in _gen()]
                txt = ''.join(recs)
            else:
                if k=='se':# save srcname, cov
                    #tgt['sname'] = path.split('/')[-1].replace('.gz','').replace('.txt','').replace('.ex2','').replace('.ex','')
                    tgt['sname'] = name
                    tmp = tgt[['chr','st','ed','sname','cov']].apply(lambda x: '\t'.join(map(str, x)),axis=1)
                    txt = '\n'.join(tmp.values)+'\n'
                else:   
                    txt = '\n'.join((tgt['chr']+'\t'+tgt['st'].astype(str)+'\t'+tgt['ed'].astype(str)).values)+'\n'
            dst[k].write(txt)
    for k, v in dst.items():
        v.close()
        UT.compress(paths[k])

def make_sj_bed_chr(sjpaths,dstpath,chrom,uth,mth):
    cols = ['chr','st','ed','strand','src','ucnt','mcnt']
    wpath = dstpath+chrom
    # n = len(sjpaths) # how many files?
    # scale = 1/float(n)
    with open(wpath,'w') as dst:
        for i, (name,spath) in enumerate(sjpaths):
            sj = GGB.read_sj(spath)
            #scale = 1e6/float(aligned)
            sj['_id'] = N.arange(len(sj))
            #name = os.path.basename(spath)[:-len('sj.txt.gz')]
            sj['src'] = name+':'+sj['_id'].astype(str)
            #sj['ucnt'] = sj['ucnt']*scale # average over all samples
            #sj['mcnt'] = sj['mcnt']*scale
            sj0 = sj[(sj['chr']==chrom)&((sj['ucnt']>=uth)|(sj['mcnt']>=mth))][cols]
            txt = '\n'.join(['\t'.join(map(str, x)) for x in sj0.values])
            dst.write(txt+'\n')
    return UT.compress(wpath)


class MergeAssemble(object):
    """Merge multiple assemblies into one. 

    1. run assembler for each strand
    2. combine outputs from each strand (ME: multi-exons)
    3. detect SE (single-exons)
    4. combine ME and SE
    5. calculate coverage based on averaged bigwig
    6. calculate junction coverage based on aggregated junctions

    """

    def __init__(self, fni, fna, saveintermediates=False, **kw):
        """
        Args:
            fni: MergeInputNames object
            fna: MergeAssemblyNames object

        Keywords:
            can be used to modify assembly parameters

        """
        self.fni = fni
        self.fna = fna
        self.params = AS.MPARAMS.copy()
        self.params.update(MERGEASMPARAM)
        self.params.update(kw)
        self.kw = kw
        self.saveintermediates = saveintermediates

        self.make_assemblers()

    def make_assemblers(self):
        fni = self.fni
        fna = self.fna
        pr = self.params
        sjexdic = {'mep': {'bw':fni.ex_bw('mep'),
                           'sj':fni.sj_bed('p')},
                   'men': {'bw':fni.ex_bw('men'),
                           'sj':fni.sj_bed('n')}}
        # FileNames objects for assemblers
        fns = {k:FN.FileNames(sname = '{0}.{1}'.format(fna.code,k),
                              bwfile = sjexdic[k]['bw'],
                              sjfile = sjexdic[k]['sj'],
                              outdir = fna.outdir,
                              refgtf = fna.refgtf) for k in sjexdic}
        savei = self.saveintermediates
        self.asms = asms = {k: AS.Assembler(fns[k], saveintermediates=savei, **pr) for k in fns}
        asms['men'].params['binstrand']='-'
        asms['mep'].params['binstrand']='+'


    def assemble(self):
        self.assemble_me1()
        self.assemble_me2()
        self.assemble_se()
        self.assemble_combine()
        self.assemble_writefiles()
        self.calc_merged_covs()
        self.assign_sjcnt()
        if not self.saveintermediates:
            self.fna.delete(delete=[],protect=['output'])
            # also delete outputs of mep, men assemblies
            self.asms['mep'].fnobj.delete(['output'])
            self.asms['men'].fnobj.delete(['output'])

    def assemble_me1(self):
        """do assembly separately for each strand"""
        asms = self.asms
        LOG.info('#########  START +strand assembly ######################################')
        asms['mep'].assemble()
        LOG.info('#########  START -strand assembly ######################################')
        asms['men'].assemble()        
        LOG.info('########################################################################')
        LOG.info('mep:{0}, men:{1}'.format(len(asms['mep'].ae),len(asms['men'].ae)))

    def _remove_se_from_me(self):
        fna = self.fna
        fnp = self.asms['mep'].fnobj
        fnn = self.asms['men'].fnobj
        exp = UT.read_pandas(fnp.ex_out())
        exn = UT.read_pandas(fnn.ex_out())
        mep,sep = UT.mese(exp)
        men,sen = UT.mese(exn)
        # FIX _gidx! 2016-03-22 
        #exp['_gidx'] = N.abs(exp['_gidx'])
        #exn['_gidx'] = -N.abs(exn['_gidx'])
        # ====> just remove all SE from this stage
        # remove SE overlapping opposite strand
        # def _select(se1,me2):
        #     cols0 = ['chr','st','ed','name','_id','strand']
        #     a = fna.fname('setmp.bed.gz')
        #     b = fna.fname('metmp.bed.gz')
        #     c = fna.fname('semetmp.bed.gz')
        #     a = UT.write_pandas(se1[cols0], a, '')
        #     b = UT.write_pandas(me2[cols0], b, '')
        #     c = PO.BT.bedtoolintersect(a,b,c,wao=True) # get original entry
        #     cols = cols0+['b_'+x for x in cols0]+['ovl']
        #     cdf = UT.read_pandas(c,names=cols)
        #     cdfg = cdf.groupby('_id')
        #     ovl = cdfg['ovl'].sum()
        #     #siz = cdfg.size()
        #     #nonovl_se = ovl[(ovl==0)|(siz>1)].index.values 
        #     nonovl_se = ovl[ovl==0].index.values 
        #     se2 = se1.set_index('_id').ix[nonovl_se].reset_index()
        #     return se2
        #sep2 = _select(sep,men)
        #sen2 = _select(sen,mep)
        #exp2 = PD.concat([mep,sep2],ignore_index=True)
        #exn2 = PD.concat([men,sen2],ignore_index=True)
        #return exp2,exn2
        return mep, men

    def assemble_me2(self):
        """combine strands for ME"""
        fna = self.fna
        fnp = self.asms['mep'].fnobj
        fnn = self.asms['men'].fnobj

        exp,exn = self._remove_se_from_me()
        sjp = UT.read_pandas(fnp.sj_out())
        sjn = UT.read_pandas(fnn.sj_out())
        self.expn = expn = PD.concat([exp,exn], ignore_index=True)
        self.sjpn = sjpn = PD.concat([sjp,sjn], ignore_index=True)
        expn['_id'] = N.arange(len(expn))
        sjpn['_id'] = N.arange(len(sjpn))
        
        # stats
        n0 = len(set(expn['_gidx'].values))
        np = len(set(exp['_gidx'].values))
        nn = len(set(exn['_gidx'].values))
        LOG.info('n0:{0}, np:{1}, nn:{2}, np+nn:{3}'.format(n0,np,nn,np+nn))
        
        # write EX,SJ
        self.ecols = ecols = ['chr','st','ed','name','sc1','strand',
                 '_id','_gidx','gname','cat','ptyp','cov','len',
                 'a_id','d_id','a_degree','d_degree','a_pos','d_pos']
        self.scols = scols = ['chr','st','ed','name','sc1','strand',
                 '_id','_gidx','gname','st-1',
                 'a_id','d_id','a_degree','d_degree','a_pos','d_pos',
                ]
        UT.write_pandas(expn[ecols], fna.fname('mepn.ex.txt.gz'), 'h')
        UT.write_pandas(sjpn[scols], fna.fname('mepn.sj.txt.gz'), 'h')
        
        # GENES BED
        gp = GGB.read_bed(fnp.genes_out())
        gn = GGB.read_bed(fnn.genes_out())
        gidp = set(exp['gname'].values)
        gidn = set(exn['gname'].values)
        gp = gp[[x in gidp for x in gp['name']]]
        gn = gn[[x in gidn for x in gn['name']]]
        self.genes = genes = PD.concat([gp,gn],ignore_index=True)
        GGB.write_bed(genes, fna.fname('mepn.genes.bed.gz'), ncols=12)

    def assemble_se(self):
        """ Calculate SE candidate (subtract ME) 
        Currently only filter with maxcov.

        This part needs to be improved.

        """
        # [TODO] also process sep, sen (stranded SEs)
        # [TODO] do power law fitting and adaptively find secov threshold
        # [Q] does power law still apply for aggregated coverages?

        fna = self.fna
        fni = self.fni
        pr = self.params

        if hasattr(self, 'expn'):
            expn = self.expn
        else:
            self.expn = expn = UT.read_pandas(fna.fname('mepn.ex.txt.gz'))
            self.sjpn = sjpn = UT.read_pandas(fna.fname('mepn.sj.txt.gz'))

        sebin = BW.bw2bed(
                    bwfile=fni.ex_bw('se'), 
                    bedfile=fna.fname('sebw0.bed.gz'), 
                    chroms=UT.chroms(pr['genome']), 
                    th=0)

        mefile = GGB.write_bed(expn, fna.fname('mepn.me.bed.gz'), ncols=3)
        sufile = BT.bedtoolintersect(sebin,mefile,fna.fname('mepn.se-me.bed.gz'),v=True) # -v subtract
        df = GGB.read_bed(sufile)

        # calculate SECOV, SEMAX
        self.secov = secov = CC.calc_cov_mp(
                                    bed=df, 
                                    bwname=fni.ex_bw('se'), 
                                    fname=fna.fname('secov.txt.gz'), 
                                    np=pr['np'], 
                                    which='cov')
        self.semax = semax = CC.calc_cov_mp(
                                    bed=df, 
                                    bwname=fni.ex_bw('se'), 
                                    fname=fna.fname('semax.txt.gz'), 
                                    np=pr['np'], 
                                    which='max')
        semax['cov'] = secov['cov']
        # threshold to get SE
        self.se0 = se0 = semax[semax['max']>pr['se_maxth']].copy()
        
        # save 
        gid0 = max(pr['se_gidstart'], N.max(N.abs(expn['_gidx'])))
        se0['_gidx'] = N.arange(gid0,gid0+len(se0))
        se0['name'] = ['S{0}'.format(x) for x in se0['_gidx']]
        se0['gname'] = se0['name']
        se0['sc1'] = se0['max']
        se0['strand'] = '.'
        sename = fna.fname('se.bed.gz',category='output')
        GGB.write_bed(se0, sename, ncols=6)

    def assemble_combine(self):
        """Combine ME/SE """
        fna = self.fna
        fni = self.fni
        pr = self.params
        ecols = self.ecols # from assemble_me2
        sjpn = self.sjpn
        expn = self.expn
        genes = self.genes
        se0 = self.se0

        # match SE to ME: ecols
        se0['_id'] = N.arange(len(expn),len(expn)+len(se0))
        se0['cat'] = 's'
        se0['ptyp'] = 's'
        se0['len'] = se0['ed']-se0['st']
        se0['a_id'] = 0
        se0['a_degree'] = 0
        se0['a_pos'] = se0['chr']+':'+se0['st'].astype(str)+':.'
        se0['d_id'] = 0
        se0['d_degree'] = 0
        se0['d_pos'] = se0['chr']+':'+se0['ed'].astype(str)+':.'
        
        self.sj0 = sj0 = sjpn
        self.ex0 = ex0 = PD.concat([expn[ecols],se0[ecols]], ignore_index=True)
        # fix id, ad info
        UT.set_ids(sj0)
        UT.set_ids(ex0)
        UT.set_ad_info(sj0,ex0)
        
        # make ci
        self.ci0 = ci0 = UT.chopintervals(ex0, fname=fna.ci_out())

        # calculate glen, tlen
        tlen = UT.calc_tlen(ex0, ci0)
        g2tlen = UT.df2dict(tlen, 'index', 'tlen')
        ex0['tlen'] = [g2tlen[x] for x in ex0['_gidx']]
        gr = ex0.groupby('_gidx')
        glen = gr['ed'].max() - gr['st'].min()
        g2glen = UT.series2dict(glen)
        ex0['glen'] = [g2glen[x] for x in ex0['_gidx']]

        # adjust genes bed
        se0['tst'] = se0['st']
        se0['ted'] = se0['ed']
        se0['sc2'] = se0['cov']
        se0['#exons'] = 1
        se0['esizes'] = se0['len'].astype(str)+','
        se0['estarts'] = '0,'
        bcols = GGB.BEDCOLS
        self.genes0 = genes0 = PD.concat([genes[bcols],se0[bcols]],ignore_index=True)

    def assemble_writefiles(self):
        """ Write EX, SJ, GENES output files """
        fna = self.fna
        scols = self.scols
        GGB.write_bed(self.sj0[scols], fna.sj_out('bed'), ncols=6)        
        UT.write_pandas(self.sj0[scols], fna.sj_out('txt'), 'h')
        
        GGB.write_bed(self.ex0, fna.ex_out('bed'), ncols=6)
        UT.write_pandas(self.ex0, fna.ex_out('txt'), 'h')
        
        GGB.write_bed(self.genes0, fna.genes_out('bed'), ncols=12)
        UT.write_pandas(self.genes0, fna.genes_out('txt'), 'h')

    def calc_merged_covs(self):
        """ calculate ecov, gcov against aggregated bigwig """
        fna = self.fna
        fni = self.fni
        pr = self.params

        expath = fna.ex_out('txt')
        cipath = fna.ci_out()
        bwpath = fni.agg_bw()
        dstpre = fna.fname('')
        covciname = fna.fname('covci.txt.gz') # register as temp files, delete later
        gcovname = fna.fname('gcov.txt.gz')
        ecovname = fna.fname('ecov.txt.gz')

        gcov = CC.calc_gcov(expath, cipath, bwpath, dstpre, override=True, np=pr['np'])
        ecov = CC.calc_ecov(expath, cipath, bwpath, dstpre, override=False, np=pr['np'])
        
        # set ecov, gcov columns
        i2g = UT.df2dict(gcov, '_gidx','gcov')
        i2e = UT.df2dict(ecov, 'eid','ecov')
        ex0 = self.ex0
        ex0['ecov'] = [i2e[x] for x in ex0['_id']]
        ex0['gcov'] = [i2g[x] for x in ex0['_gidx']]

        # set cov column for genes
        genes0 = self.genes0
        # _gidx not in genes0
        def name2gidx(s):
            if s[0]=='N':
                return -int(s[2:])
            if s[0]=='P':
                return int(s[2:])
            return int(s[1:])
        genes0['_gidx'] = [name2gidx(x) for x in genes0['name']]
        genes0['cov'] = [i2g[x] for x in genes0['_gidx']]
        
        # overwrite ex0, genes0
        UT.write_pandas(ex0, fna.ex_out('txt'), 'h')
        UT.write_pandas(genes0, fna.genes_out('txt'), 'h')

    def assign_sjcnt(self):
        """ calculate junction counts """
        fna = self.fna
        fni = self.fni

        sjg = GGB.read_sj(fni.sj0_bed())
        sjg['locus'] = UT.calc_locus_strand(sjg)

        l2u = UT.df2dict(sjg, 'locus','ucnt')
        l2m = UT.df2dict(sjg, 'locus','mcnt')
        if hasattr(self, 'sj0'):
            sj0 = self.sj0 
        else:
            self.sj0 = sj0 = UT.read_pandas(fna.sj_out('txt'))
        if 'locus' not in sj0.columns:
            sj0['locus'] = UT.calc_locus_strand(sj0)
        sj0['ucnt'] = [l2u.get(x,0) for x in sj0['locus']]
        sj0['mcnt'] = [l2m.get(x,0) for x in sj0['locus']]
        sj0['jcnt'] = [x or y for x,y in sj0[['ucnt','mcnt']].values]
        # overwrite sj0
        UT.write_pandas(sj0, fna.sj_out('txt'), 'h')





