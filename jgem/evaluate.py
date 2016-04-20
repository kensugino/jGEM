"""
.. module:: evaluate
    :synopsis: evaluate performance by comparing to a reference annotation

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""

# system imports
import gzip
import os
import subprocess
from collections import Counter
from operator import iadd
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

# 3rd party libraries
import pandas as PD
import numpy as N
import matplotlib.pyplot as P

# library imports
from jgem import utils as UT
from jgem import gtfgffbed as GGB
from jgem import bedtools as BT
from jgem import bigwig as BW
from jgem import filenames as FN
from jgem import calccov as CC

class EvalNames(FN.FileNamesBase):
    """Filename manager for evaluation process.

    Attributes:
        sjexprefix: path prefix to junction, exon files (*.sj.txt.gz and *.ex.txt.gz)
        code: assembly identifier
        outdir: output directory

    All output and temporary files are prefixed by **outdir/code**

    """

    def __init__(self, sjexprefix, code, outdir):
        self.sjexbase = sjexprefix
        self.code = code
        self.outdir = outdir
        for x in ['sj','ex','ci']:
            setattr(self, x+'path', '{0}.{1}.txt.gz'.format(sjexprefix,x))

        prefix = os.path.join(outdir, code)
        super(EvalNames, self).__init__(prefix)

    def fname2(self, suffix, code2, category='temp'):
        """Generate filenames furthre prefixed by code2.

        Args:
            suffix: (str)
            code2: (str) identifier (comparison target)
            category: (str)

        Returns:
            (outdir)/(code).(code2).(suffix)

        """
        suf = '{0}.{1}'.format(code2, suffix)
        return self.fname(suf, category)

    def modelpath(self, which):
        """Returns path to junction(sj)/exon(ex)/choppedinterval(ci) file.

        Args:
            which: one of 'sj','ex','ci'

        """
        return '{0}.{1}.txt.gz'.format(self.sjexbase, which)

    def model(self, which):
        """Returns model dataframe (junction/exon/chopped intervals).

        Args:
            which: one of 'sj','ex', 'ci'

        """
        if hasattr(self, which): # cached
            return getattr(self, which)

        path = self.modelpath(which)
        if os.path.exists(path): # file exists
            df = UT.read_pandas(path)
            setattr(self, which, df)
            return df
        # file does not exists, if ci then make from ex
        if which=='ci':
            expath = self.modelpath['ex']
            if os.path.exists(expath):
                self.ci = UT.chopintervals(self.model['ex'], path)
            else:
                raise RuntimeError('file {0} does not exist'.format(expath))
        else:
            raise RuntimeError('file {0} does not exist'.format(path))
        
    def savemodel(self, which, code2=None, category='temp'):
        """Save model. If code2 is None, overwrite original, if code2 is provided,
        writes to outdir/(code).(code2).(which).txt.gz. 

        Args:
            which: 'sj','ex','ci'
            code2: 2nd identifier
            category: filename category (default 'temp')

        Returns:
            file path or None (if model is not loaded)

        """
        if hasattr(self, which):
            if code2 is None:
                path = self.modelpath(which)
            else:
                path = self.fname2('{0}.txt.gz'.format(which),code2)
            return UT.write_pandas(getattr(self, which), path, 'h')
        return None

class EvalMatch(object):
    """Compare two models against a genome coverage (bigwig) 
    and junction counts (sjfile).

    Usage:
        >>> en1 = EvalNames(sjexpre_to_ref, 'ref', outdir)
        >>> en2 = EvalNames(sjexpre_to_target, 'tgt', outdir)
        >>> em = EvalMatch(en1,en2,bwfile,sjfile,datacode)
        >>> figs = em.calculate(outdir)

    """
    
    abbr = {'i':'internal exons',
            '5':"5' exons",
            '3':"3' exons",
            's':'single exons',
            'j':'junctions'}

    def __init__(self, en1, en2, bigwig, sjfile, datacode, binsize=1000):
        """
        Args:
            en1: EvalNames object, reference
            en2: EvalNames object, sensitivity of this model against en1 is calculated
            bigwig: path to normalized bigwig coverage file
            sjfile: path to normalized junction counts file
            datacode: code indicating data (bigwig & sjfile)
            binsize (int): for sensitivity plot (default 1000)

        """
        self.en1 = en1
        self.en2 = en2
        self.bigwig = bigwig
        self.sjfile = sjfile
        self.datacode = datacode
        self.closest = {}       
        self.stats = {}     
        self.binsize = binsize

    def calculate(self, np=1):
        """Calculate necessary data.

        1. for en1 and en2 calculate ecov,gcov,jcnt (prep_sjex)
        2. calculate match between en1 and en2 (find_match)
        3. calculate length ratio, detected numbers, sensitivity, etc. (calc_stats)

        """
        # calc exon, junction, gene coverage
        self.prep_sjex(self.en1, np)
        self.prep_sjex(self.en2, np)
        self.find_match()
        self.calc_stats()

    def save_stats(self):
        pass


    def colname(self, x):
        return '{0}_{1}'.format(x, self.datacode)

    def colname2(self, x, code):
        return '{0}_{1}_{2}'.format(x, self.datacode, code)


    def prep_sjex(self, en, np=1):
        """ Assign ecov, gcov, jcnt """
        sj = en.model('sj')
        ex = en.model('ex')
        savesj = False
        saveex = False
        # length
        if 'len' not in sj.columns:
            sj['len'] = sj['ed'] - sj['st']
            savesj = True
        if 'len' not in ex.columns:
            ex['len'] = ex['ed'] - ex['st']
            saveex = True
        # ecov
        ecovname = self.colname('ecov')
        if ecovname not in ex.columns:
            ecov = CC.calc_ecov(
                expath=en.modelpath('ex'), 
                cipath=en.modelpath('ci'), 
                bwpath=self.bigwig, 
                dstprefix=en.sjexbase, 
                override=False, np=np)
            ex[ecovname] = ecov.set_index('eid').ix[ex['_id'].values]['ecov'].values
            saveex = True
        # gcov, glen
        gcovname = self.colname('gcov')
        if gcovname not in ex.columns:
            gcov = CC.calc_gcov(
                expath=en.modelpath('ex'), 
                cipath=en.modelpath('ci'), 
                bwpath=self.bigwig, 
                dstprefix=en.sjexbase, 
                override=False, np=np)
            tmp = gcov.set_index('_gidx').ix[ex['_gidx'].values]
            ex[gcovname] = tmp['gcov'].values
            ex['glen'] = tmp['glen'].values # glen is only dependent on model not data
            saveex = True
        # sjcnt
        ucntname = self.colname('ucnt')
        mcntname = self.colname('mcnt')
        jcntname = self.colname('jcnt')
        sjfile = self.sjfile
        if ucntname not in sj.columns:
            if sjfile.endswith('.bed') or sjfile.endswith('.bed.gz'): # no header
                dsj = UT.read_pandas(sjfile, names=['chr','st','ed','name','ucnt','strand','mcnt'])
            else: # assume txt file with header
                dsj = UT.read_pandas(sjfile) 
            # locus based matching
            dsj['locus'] = UT.calc_locus_strand(dsj)
            sj['locus'] = UT.calc_locus_strand(sj)
            l2u = UT.df2dict(dsj, 'locus', 'ucnt')
            l2m = UT.df2dict(dsj, 'locus', 'mcnt')
            sj[ucntname] = [l2u.get(x,0) for x in sj['locus']]
            sj[mcntname] = [l2m.get(x,0) for x in sj['locus']]
            sj[jcntname] = [x or y for x,y in sj[[ucntname,mcntname]].values]
            savesj = True
        if saveex:
            en.savemodel('ex',self.datacode)
        if savesj:
            en.savemodel('sj',self.datacode)

    def find_match(self):
        en1 = self.en1
        en2 = self.en2
        # write internal,3,5,se exons separately for finding match
        a = en1.fname('ex.bed.gz')
        b = en2.fname('ex.bed.gz')
        c = en1.fname2('ex.ovl.txt.gz', en2.code)
        self.e1 = e1 = en1.model('ex')
        self.e2 = e2 = en2.model('ex')
        ecovname = self.colname('ecov')
        cols = ['chr','st','ed','cat','_id',ecovname,'_gidx','len','strand']
        a = UT.write_pandas(e1[cols],a,'')
        b = UT.write_pandas(e2[cols],b,'')
        c = BT.bedtoolintersect(a,b,c,wao=True)
        ocols = cols + ['b_'+x for x in cols] + ['ovl']
        self.ov = ov = UT.read_pandas(c, names=ocols) # overlaps of exons
        
        idxchr = ov['chr']==ov['b_chr'] 
        idxstrand = ov['strand']==ov['b_strand']
        idxp = (ov['strand']=='+')&idxstrand
        idxn = (ov['strand']=='-')&idxstrand
        idxst = ov['st']==ov['b_st']
        idxed = ov['ed']==ov['b_ed']
        idxcat = ov['cat']==ov['b_cat']
        idxcov = ov[ecovname]>0 # exons with reads
        LOG.info('='*10 + 'calculating match between {0} and {1}'.format(en1.code, en2.code))
        LOG.info('len(ov):{0}'.format(len(ov)))
        for k in ['idxchr','idxstrand','idxp','idxn','idxst','idxed','idxcat','idxcov']:
            v = locals()[k]
            LOG.info('#{0}:{1}'.format(k, N.sum(v)))
        
        # internal exon cat='i' and chr,st,ed,strand match
        self.ei = ei = ov[idxchr&idxstrand&idxst&idxed&idxcat&(ov['cat']=='i')].copy()
        # 5' cat='5' and chr,donor (+,ed)|(-,st) match, find closest
        self.e5 = e5 = ov[idxchr&((idxp&idxed)|(idxn&idxst))&idxcat&(ov['cat']=='5')].copy()
        # 3' cat='3' and chr,acceptor (+,st)|(-,ed) match
        self.e3 = e3 = ov[idxchr&((idxn&idxed)|(idxp&idxst))&idxcat&(ov['cat']=='3')] .copy()
        # se cat='s' and chr,
        self.es = es = ov[idxchr&(ov['cat']=='s')&idxcat].copy()
        
        # splice junction
        self.s1 = s1 = en1.model('sj')
        self.s2 = s2 = en2.model('sj')
        jcntname = self.colname('jcnt')
        l2c = UT.df2dict(s2, 'locus',jcntname)
        jhitname = self.colname2('jhit', en2.code)
        s1[jhitname] = [l2c.get(x,0) for x in s1['locus']] # corresponding s2 count
        self.sj= sj = s1[s1[jhitname]>0].copy() # only consider s2 count > 0
        
        # for batch processing
        self.e = {'i':ei,'5':e5,'3':e3,'s':es, 'j':sj}
        
    def calc_stats(self):

        ecovname = self.colname('ecov')
        jcntname = self.colname('jcnt')
        jhitname = self.colname2('jhit', self.en2.code)

        def _findclosest(e, which):
            e['dlen'] = N.abs(e['len']-e['b_len'].astype(int))
            e['ratio'] = e['b_len'].astype(float)/e['len']
            e = e.sort_values(['_id','dlen'],ascending=True)
            f = e.groupby('_id',sort=False).first().reset_index()
            self.closest[which] = f
            return f

        def _calc(x0,y0):
            avgx,avgy,minx,maxx,cnt = UT.calc_binned(x0, y0, num=self.binsize, returnminmax=True)
            LOG.debug('len(avgx)={0},len(avgy)={1},len(minx)={2},len(maxx)={3}'.
                format(len(avgx),len(avgy),len(minx),len(maxx)))
            avgy1 = N.concatenate([avgy,[0]])
            delta = maxx - minx
            hight = (avgy+avgy1[1:])/2.
            if len(maxx)>0:
                auc = N.sum(delta*hight)/(maxx[0]-minx[-1])
            else:
                auc = 0.
            LOG.debug('len(x0)={0},len(y0)={1}, auc={2:.3f}'.format(len(x0),len(y0),auc))
            return auc,maxx,avgy,x0,y0

        def _count(dw, da, which):
            hit = set(dw['_id'].values)
            if which != 'j':
                da = da[da[ecovname]>0]
            else:
                da = da[da[jcntname]>0]
            pop = set(da['_id'].values)
            dif = pop.difference(hit)
            if len(pop)==0:
                LOG.warning('no elements in {0} for population'.format(self.abbr[which]))
            if len(hit)==0:
                LOG.warning('no elements in {0} for hits'.format(self.abbr[which]))
            LOG.info( '[{4}] pop:{0},hit:{1},dif:{2},hit/pop:{3:.2f}'.format(
                len(pop),len(hit),len(dif), float(len(hit))/max(1,len(pop)), which) )
            return hit, pop, dif


        for which in ['i','5','3','s','j']:
            LOG.info(which+'='*10)
            cn = 'hit{0}'.format(which)
            if which != 'j':
                e1 = self.e1
                # use exons with reads
                ea = e1[(e1['cat']==which)][['_id',ecovname]].copy() # all exons
                ew = self.e[which] # matched exons
                hit, pop, dif = _count(ew, ea, which)
                ew2 = _findclosest(ew, which) # calculate ratio
                i2r = UT.df2dict(ew2,'_id','ratio')
                ea[cn] = [i2r.get(x,0) for x in ea['_id']]
                x = N.log2(ea[ecovname].values+1) # log coverage
                y = ea[cn].values
            else:
                sa = self.s1 
                sw = self.e['j']
                hit, pop, dif = _count(sw, sa, which)
                sa[cn] = [1 if x>0 else 0 for x in sa[jhitname]] # in case of NaN
                x = N.log2(sa[jcntname].values+1)
                y = sa[cn].values
                
            # # positive ratio
            # idx = y>0
            # y0 = y[idx]
            # x0 = x[idx]                
            # auc,maxx,avgy,x0,y0 = _calc(x0,y0)
            # # gen4 ecov>0, ratio
            # idx2 = x>0 
            # y2 = y[idx2]
            # x2 = x[idx2]
            # auc2,maxx2,avgy2,x2,y2 = _calc(x2,y2)

            # gen4 ecov>0, detected or not
            # idx2 = x>0
            # x2 = x[idx2]
            # y4 = N.array(y[idx2]>0, dtype=int)
            x2 = x
            y4 = N.array(y>0, dtype=int)
            auc4,maxx4,avgy4,x4,y4 = _calc(x2,y4)

            #self.stats[which] = (pop,hit,dif,auc,maxx,avgy,x,y,auc2,maxx2,avgy2,x2,y2,auc4,maxx4,avgy4,x4,y4)
            self.stats[which] = {'pop':pop,'hit':hit,'dif':dif,
                                'x':x,'y':y,'xpd':x4,'ypd':y4,
                                'auc':auc4,'maxx':maxx4,'avgy':avgy4,}

    def calc_completeness(self):
        """Completeness measures how much of the reference structure is recovered.

        1. ELC: exon length completeness = max(ratio of exon length covered by overlapping target gene)
        2. GLC: gene length completeness = max(ratio of gene length covered by overlapping target gene)
        3. ECC: exon count completeness = max(ratio of overlapping exon counts)
        4. JCC: junction count completeness = max(ratio of overlapping junction counts)


        """

    def plot_sensitivity(self, color='b', ypos=0, axr=None, lineonly=False):
        ws = ['i','5','3','s','j']
        p1c = self.en1.code # gen4
        p2c = self.en2.code

        def _plot_one(ax, which, label):
            s = self.stats[which]
            ax.plot(s['maxx'],100*s['avgy'],color+'.-',ms=5, label=label)
            ma = int(N.max(s['maxx']))
            ax.set_xlim([-0.5,ma+0.5])
            ax.text(ma/2,10*(1+ypos),'{0}:{1:.2f}'.format(label,s['auc']))

        if axr is None:
            fig,axr = P.subplots(1,len(ws),figsize=(3*len(ws),3),sharey=True)
            P.subplots_adjust(wspace=0.07, top=0.85)
        else:
            assert len(axr)==len(ws)
            fig = None

        for i,w in enumerate(ws):
            ax = axr[i]
            _plot_one(ax, w, p2c)
            if not lineonly:
                ax.set_title(self.abbr[w])
                if w!='j':
                    ax.set_xlabel('log2({0}_ecov+1)'.format(p1c))
                else:
                    ax.set_xlabel('log2({0}_jcnt+1)'.format(p1c))
        if not lineonly:
            axr[0].set_ylim([-5,105])
            axr[0].set_ylabel('%detected')
            axr[len(ws)-1].legend(loc='center left', bbox_to_anchor=(1.0,0.5))

        if fig is not None:
            fig.suptitle('{1}/{0}'.format(p1c,p2c))
        return axr

    def plot_ratio(self,axr=None,plotxlabel=True,label='',disp='both'):
        """Plot length ratios of best matching exons """
        p1c = self.en1.code # gen4
        p2c = self.en2.code
        tgts = ['5','3','s']
        if axr is None:
            fig,axr = P.subplots(1,len(tgts),figsize=(3*len(tgts),3),sharex=True,sharey=True)
            P.subplots_adjust(wspace=0.07,hspace=0.15,top=0.85)
        else:
            fig = None
        for i,w in enumerate(tgts):
            ax = axr[i]
            #pop,hit,dif,auc,maxx,avgy,x,y = self.stats[w]
            st = self.stats[w]
            auc,maxx,avgy,x,y = st['auc'],st['maxx'],st['avgy'],st['x'],st['y']
            if disp!='pdf':
                ax.plot(x,y,'.',ms=3, alpha=0.3)
            #ax.plot(maxx,avgy,'ro-',ms=3,alpha=0.3)
            ax.set_yscale('log')
            if disp!='png':
                if i==0:
                    ax.set_ylabel('{1}_len/{0}_len+1'.format(p1c,p2c))
                if plotxlabel:
                    ax.set_xlabel('log2({0}_ecov+1)'.format(p1c))
                ax.set_title(label+self.abbr[w])
                m = 10**(N.nanmean(N.log10(y[(x>0)&(y>0)])))
                ax.text(5,10**2,'avg:{0:.2f}'.format(m))
            else:
                ax.set_yticks([])
                ax.set_xticks([])
        if fig is not None:
            fig.suptitle('{1}/{0}'.format(p1c,p2c))
        return axr
        


class CompareMatches(object):

    def __init__(self, base, targets, bigwig, sjfilel, datacode, binsize=1000):
        """
        Args:
            base: EvalNames object
            targets: listl of EvalNames objects

        """
        self.base = base
        self.targets = targets
        self.matches = {x.code:EvalMatch(base,x,bigwig,sjfile,datacode,binsize) for x in targets}
        self.codes = self.matches.keys()
        for v in self.matches.values():
            v.calculate()



    def make_ratio_suppfig(self,prefix,base='gen4',tgts=['pndr','cmst','cmst2']):
        ws = ['5','3','s']
        # pdf w/o dots
        fig1,axr = P.subplots(len(tgts),len(ws),figsize=(3*len(ws),3*len(tgts)),sharey=True,sharex=True)
        P.subplots_adjust(wspace=0.07,hspace=0.07)
        for i,t in enumerate(tgts):
            em = self.ems[base][t]
            em.plot_ratio(axr[i],plotxlabel=False,label=t+' ',disp='pdf')
        axr[len(tgts)-1][1].set_xlabel('log2(gen4_ecov+1)')
        # png 
        fig2,axr = P.subplots(len(tgts),len(ws),figsize=(3*len(ws),3*len(tgts)),sharey=True,sharex=True)
        P.subplots_adjust(wspace=0.07,hspace=0.07)
        for i,t in enumerate(tgts):
            em = self.ems[base][t]
            em.plot_ratio(axr[i],plotxlabel=False,label=t+' ',disp='png')
        fig1.savefig(prefix+'.pdf')
        fig2.savefig(prefix+'.png', dpi=300)
        return fig1,fig2
        



