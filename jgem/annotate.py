"""

.. module:: compare
    :synopsis: compare to a reference and annotate known genes

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
# system imports
import os
from collections import Counter
from operator import iadd
from functools import reduce

# 3rd party libraries
import pandas as PD
import numpy as N
#import matplotlib.pyplot as P

# library imports
from jgem import utils as UT
from jgem import gtfgffbed as GGB
from jgem import bedtools as BT
from jgem import filenames as FN 
from jgem import evaluate as EV

class ComparatorNames(EV.EvalNames):
    """Filename manager for evaluation process. Same as EvalNames

    Attributes:
        sjexbase: path prefix to junction, exon files (\*.sj.txt.gz and \*.ex.txt.gz)
        code: assembly identifier
        outdir: output directory

    All outputs and temporary files are prefixed by **outdir/code**

    """
    def __init__(self, sjexbase, code, outdir):
        super(ComparatorNames, self).__init__(sjexbase, code, outdir)

    
class Comparator(object):
    """ Compare to reference and annotate. 

    """
    def __init__(self, cn_ref, cn_tgt, gnamecol='gene_name', gidxcol='gene_id'):
        """
        Args:
            cn_ref: Reference ComparatorNames object
            cn_tgt: Target ComparatorNames object
            gnamecol: name of the column containing gene names (gene_name for Gencode, symbol for RefSeq)
             (default gene_name)
            gidxcol: name of the column containing gene id (gene_id for Gencode, 
             _gidx if using connected components as genes)

        """
        self.cn_ref = cn_ref
        self.cn_tgt = cn_tgt
        self.refgnamecol = gnamecol
        self.refgidxcol = gidxcol

    def calc_overlaps(self):
        cref = self.cn_ref
        ctgt = self.cn_tgt
        a = ctgt.fname('cptmp.ex.bed.gz')
        b = cref.fname('cptmp.ex.bed.gz')
        c = ctgt.fname2('cptmp.ex.ovl.txt.gz', cref.code)
        cols = ['chr','st','ed','cat','_id','_gidx','len','strand']
        self.ex_tgt = etgt = ctgt.model('ex') #UT.read_pandas(p1.ex)
        self.ex_ref = eref = cref.model('ex') #UT.read_pandas(p2.ex)
        
        eref['_gidx'] = eref[self.refgidxcol]

        if 'len' not in etgt.columns:
            etgt['len'] = etgt['ed']-etgt['st']
        if 'len' not in eref.columns:
            eref['len'] = eref['ed']-eref['st']
        a = UT.write_pandas(etgt[cols],a,'')
        b = UT.write_pandas(eref[cols],b,'')
        c = BT.bedtoolintersect(a,b,c,wao=True)
        ocols = cols + ['b_'+x for x in cols] + ['ovl']
        self.ov = UT.read_pandas(c, names=ocols)

        # gene overlap
        gcols = ['chr','st','ed','strand']
        def _gbed(ex):
            gr = ex.groupby('_gidx')
            g = gr[gcols].first()
            g['st'] = gr['st'].min()
            g['ed'] = gr['ed'].max()
            return g.reset_index()
        gtgt = _gbed(etgt)
        gref = _gbed(eref)
        gcols2 = gcols+['_gidx']
        a2 = ctgt.fname('cptmp.gene.bed.gz')
        b2 = cref.fname('cptmp.gene.bed.gz')
        c2 = ctgt.fname2('gene.ovl.txt.gz', cref.code)
        a2 = UT.write_pandas(gtgt[gcols2],a2,'')
        b2 = UT.write_pandas(gref[gcols2],b2,'')
        c2 = BT.bedtoolintersect(a2,b2,c2,wao=True)
        gocols = gcols2 + ['b_'+x for x in gcols2] + ['ovl']
        self.gov = UT.read_pandas(c2,names=gocols)
        
    def assign_tcode_ex(self):
        ov = self.ov
        # restrict to proper ones (e.g. strand)
        idxchr = ov['chr']==ov['b_chr'] 
        idxstrand = ov['strand']==ov['b_strand'] # strand should be same
        # above will remove all known SE (which have '.'), add them back
        idxstrand = idxstrand|(ov['strand']=='.') 
        self.eo = eo = ov[idxchr&idxstrand].copy() # exon overlap subset
        eo['dlen'] = N.abs(eo['len']-eo['b_len'].astype(int)) # diff(length)
        eo = eo.sort_values(['_id','dlen'],ascending=True) # chose closest in length
        self.eog = eog = eo.groupby('_id',sort=False).first().reset_index()
        # exons: eknown, etcode
        ex = self.ex_tgt
        e2k = dict([(x,'k') for x in eog['_id']]) # overlapping exons
        rcode = self.cn_ref.code
        ekfld = 'eknown_'+rcode
        etfld = 'etcode_'+rcode
        ex[ekfld] = [e2k.get(x,'u') for x in ex['_id']]
        idxse = ex['cat']=='s'
        idxek = ex[ekfld]=='k'
        ex.loc[idxek&idxse,etfld] = 'k.se'
        ex.loc[idxek&(~idxse),etfld] = 'k.me'
        ex.loc[(~idxek)&idxse,etfld] = 'u.se'
        ex.loc[(~idxek)&(~idxse),etfld] = 'u.me'

    def assign_antisense(self):
        """Check whether exons overlap in anti-sense direction.
        Assign code (ex_as_ovl_(refcode)) to exon-wise and then gene-wise (gene_as_ovl_(refcode)).
        """
        ov = self.ov
        # restrict to anti-sense
        idxchr = ov['chr']==ov['b_chr']
        idxstrand = ov['strand']!=ov['b_strand'] # opposite strand 
        # above will include all exons with strand '.' => exclude
        idxstrand = idxstrand&(ov['strand']!='.')
        self.eoas = eoas = ov[idxchr&idxstrand].copy() 
        # take closest in length
        eoas['dlen'] = N.abs(eoas['len']-eoas['b_len'].astype(int))
        eoas = eoas.sort_values(['_id','dlen'],ascending=True)
        self.eoasg = eoasg = eoas.groupby('_id',sort=False).first().reset_index()
        # exons: eknown, etcode
        ex = self.ex_tgt
        e2k = dict([(x,'y') for x in eoasg['_id']])
        rcode = self.cn_ref.code
        easfld = 'ex_as_ovl_'+rcode
        gasfld = 'gene_as_ovl_'+rcode
        ex[easfld] = [e2k.get(x,'n') for x in ex['_id'].values]
        gasovl = ex.groupby('_gidx')[easfld].apply(lambda x: 'y' in set(x))
        g2k = dict([(x,'y') for x in gasovl[gasovl].index.values])
        ex[gasfld] = [g2k.get(x,'n') for x in ex['_gidx'].values]

    def assign_intergenic(self):
        """Check whether gene is intergenic or not """
        gov = self.gov  
        ex = self.ex_tgt
        gov2 = gov[gov['ovl']>0]
        # ['chr','st','ed','strand','_gidx','b_chr','b_st','b_ed','b_strand','b__gidx','ovl']
        ovlgids = set(gov2['_gidx'].values)
        g2o = dict([(x,'n') for x in ovlgids])
        rcode = self.cn_ref.code
        itgfld = 'intergenic_'+rcode
        ex[itgfld] = [g2o.get(x,'y') for x in ex['_gidx'].values]
        
    def assign_tcode_gene(self):
        eog = self.eog # overlapping exons
        ex = self.ex_tgt
        # genes: gknown, gtcode
        g2k = dict([(x,'k') for x in set(eog['_gidx'].values)])
        rcode = self.cn_ref.code
        gkfld = 'gknown_'+rcode
        gtfld = 'gtcode_'+rcode
        ex[gkfld] = [g2k.get(x,'u') for x in ex['_gidx']]
        idxse = ex['cat']=='s'
        idxgk = ex[gkfld]=='k'
        ex.loc[idxgk&idxse,gtfld] = 'k.se'
        ex.loc[idxgk&(~idxse),gtfld] = 'k.me'
        ex.loc[(~idxgk)&idxse,gtfld] = 'u.se'
        ex.loc[(~idxgk)&(~idxse),gtfld] = 'u.me'

        # ref_gidxs, ref_gidx0
        self.eogg = eogg = eog.groupby('_gidx')['b__gidx'] # overlapping ref gene _gidx's
        self.g2cnt = g2cnt = {k:Counter(list(v)) for k,v in eogg}
        self.g2gi = g2gi = {k:','.join(map(str,v.keys())) for k,v in g2cnt.items()} # all of ovl _gidx
        self.g2gi0 = g2gi0 = {k:v.most_common()[0][0] for k,v in g2cnt.items()} # most common _gidx
        ex[rcode+'_gidxs'] = [g2gi.get(x, N.nan) for x in ex['_gidx']]
        ex[rcode+'_gidx0'] = [g2gi0.get(x, N.nan) for x in ex['_gidx']]

        # ref_symbols, ref_symbol0
        sfld = self.refgnamecol
        # self.e2g = e2g = self.ex_ref.groupby('_gidx')[sfld].apply(lambda x: list(set(x))).reset_index()
        # self.g2s = g2s = UT.df2dict(e2g, '_gidx', sfld) # gidx => list of syms
        self.g2s = g2s = self.ex_ref.groupby('_gidx')[sfld].apply(lambda x: list(set(x)))
        # self.g2s = g2s = UT.df2dict(e2g, '_gidx', sfld) # gidx => list of syms
        # convert _gidx => sym
        try:
            self.g2gs = g2gs = {k:','.join(set(reduce(iadd,[g2s[int(z)] for z in v.keys()],[]))) for k,v in g2cnt.items()}
            self.g2gs0 = g2gs0 = {k:','.join(g2s[int(v.most_common()[0][0])]) for k,v in g2cnt.items()}
        except:
            self.g2gs = g2gs = {k:','.join(set(reduce(iadd,[g2s[z] for z in v.keys()],[]))) for k,v in g2cnt.items()}
            self.g2gs0 = g2gs0 = {k:','.join(g2s[v.most_common()[0][0]]) for k,v in g2cnt.items()}
        #g2gs0 = {k:g2s[int(g2gi0[k])] for k,v in g2cnt.items()}
        ex[rcode+'_syms'] = [g2gs.get(x, N.nan) for x in ex['_gidx']]
        ex[rcode+'_sym0'] = [g2gs0.get(x,N.nan) for x in ex['_gidx']]
        
    def assign_tcode_sj(self):
        self.sj_tgt = stgt = self.cn_tgt.model('sj') #UT.read_pandas(self.p1.sj)
        self.sj_ref = sref = self.cn_ref.model('sj') #UT.read_pandas(self.p2.sj)
        if 'locus' not in stgt.columns:
            stgt['locus'] = GGB.calc_locus_strand(stgt)
        if 'locus' not in sref.columns:
            sref['locus'] = GGB.calc_locus_strand(sref)
        l2c = dict([(x,'k.me') for x in sref['locus']])
        rcode = self.cn_ref.code
        setfld = 'etcode_'+rcode
        sgtfld = 'gtcode_'+rcode
        stgt[setfld] = [l2c.get(x,'u.me') for x in stgt['locus']]
        g2c = UT.df2dict(self.ex_tgt, '_gidx', 'gtcode_'+rcode)
        stgt[sgtfld] = [g2c.get(x,'u.me') for x in stgt['_gidx']]
                
    def annotate(self, save=True, overwrite=True):
        """Annotate target by overlaps to reference.

        Args:
            overwrite (bool): overwrite original sj,ex files (default True)

        """
        self.calc_overlaps()
        self.assign_tcode_ex()
        self.assign_intergenic()
        self.assign_antisense()
        self.assign_tcode_gene()
        self.assign_tcode_sj()
        cntgt = self.cn_tgt
        cnref = self.cn_ref
        if save:
            if overwrite:
                cntgt.savemodel('ex', category='output')
                cntgt.savemodel('sj', category='output')
            else:
                refcode = cnref.code
                cntgt.savemodel('ex', code2=refcode, category='output')
                cntgt.savemodel('sj', code2=refcode, category='output')
        cntgt.delete(['temp'],protect=['output'])
        cnref.delete(['temp'],protect=['output'])


        
   
#### REFEFENCE ANNOTATION RELATED #######################################################

def gencode_remove_badSE(ex):
    """Gencode put some single exons (SEs) together with overlapping gene, then annotate these 
    SEs with annotation of the overlapping gene, for example, gene_type as "protein_coding", 
    when they are not protein_coding at all. This messes things up when calculating threshold
    to distinguish protein_coding and non_coding using PhyloCSF, etc. 

    This function detect these SEs and returns a set of _gidx's which does not contain these
    bad SEs.

    """
    tmp = ex.groupby(['gene_name','_gidx']).first().reset_index()
    tmpsize = tmp.groupby('gene_name').size()

    gids1 = tmpsize[tmpsize==1].index.values # set of unique ones
    gids2 = tmpsize[tmpsize>1].index.values # have unconnected subcomponents=>clean up these

    tgts = tmp[tmp['gene_name'].isin(gids2)] 
    tgts2 = tgts[tgts['cat']!='s'] # remove SE which belong to another gene

    gidx1 = tmp[tmp['gene_name'].isin(gids1)]['_gidx'].values
    gidx2 = tgts2['_gidx'].values
    gidx0 = sorted(set(list(gidx1)+list(gidx2)))
    return gidx0   