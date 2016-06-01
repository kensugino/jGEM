"""

.. module:: assembler
    :synopsis: assemble genes from RNASeq data (normalized genome coverage (bigwig) and junctions)

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
# system imports
import gzip
import os
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
# class
from ponder import plotutils as PU

def make_dict(df,f1,f2):
    dic = {}
    for k,v in df[[f1,f2]].values:
        dic.setdefault(k,set()).add(v)
    dic = {k:list(dic[k]) for k in dic}
    return dic

class SpecificByDM(object):
    """
    Args:
        si: sample info dataframe
        gcov: gene coverage (normalized)
        gcovlevel: column name in sample info dataframe containing 
         gcov column entries (usually name)
        maxeth: min maxexpression required (default 1)
        zth1: normed logdiff smaller than this will be ignored for pairs whose minmin < minth
        zth2: normed logdiff smaller than this will be ignored for pairs whose minmin >= minth
        mmth: threshold to decide a pair is both expressed (min of min)
        rdth: threshold for calculating RD (reproducible detectable ratio)
        gbed: for annotation
        
    """
    def __init__(self, si, gcov, gbed, gcovlevel='name',maxeth=1,zth1=0.2,zth2=0.6,mmth=1,rdth=0.5):
        # first make sure si is limited to gcov
        self.snames = snames = [x for x in si[gcovlevel] if x in gcov.columns]
        print('#sname={0}, #si={1}, #gcov={2}'.format(len(snames),len(si),gcov.shape[1]))
        self.si = si[si[gcovlevel].isin(snames)].copy()

        self.gcov = gcov[snames] # use intersection
        self.gcovlevel = gcovlevel # sampleinfo column name for gcov columns
        self.gbed = gbed

        self.maxeth = maxeth
        self.zth1 = zth1
        self.zth2 = zth2
        self.mmth = mmth
        self.rdth = rdth
        
        self.dms = {} # holds dm
        
    
    def make_dm(self, targetlevel):
        """calculate 2 DMs (logdiff and minmin) at specified level """
        # first make gcovlevel <=> targetlevel mapping
        si = self.si
        gl = self.gcovlevel
        gc = self.gcov
        ts = si.groupby(targetlevel, sort=False).first().index.values
        
        g2t = UT.df2dict(si, gl, targetlevel)
        t2g = make_dict(si, targetlevel, gl)
        
        lgc = N.log2(gc+1)
        v0 = lgc.groupby(g2t, axis=1).mean() # target level
        maxe = v0.max(axis=1)
        gids = maxe[maxe>N.log2(self.maxeth+1)].index.values
        v = v0.ix[gids][ts] # restrict to expressed
        # do the math in numpy to get normalized logdiff DM
        m = v.values
        logdiff = N.abs(m[:,:,N.newaxis]-m[:,N.newaxis,:])
        maxdiff = logdiff.max(axis=2).max(axis=1)
        normdiff = logdiff/maxdiff[:,N.newaxis,N.newaxis] # normalized
        dm = PD.Panel(normdiff, v.index, ts, ts) 
        # calculate minmin DM
        gmin = gc.ix[gids].groupby(g2t, axis=1).min()[ts].values
        a = gmin[:,:,N.newaxis] # i
        b = gmin[:,N.newaxis,:] # j
        minmin = N.minimum(a,b)
        mm = PD.Panel(minmin, v.index, ts, ts)
        self.dms[targetlevel] = dict(ts=ts,g2t=g2t,t2g=t2g,dm=dm,mm=mm,v=v)

    def get_snames(self, targetlevel, levelnames):
        si = self.si
        # [levelnames,...] in targetlevel column (group or cg1) => sample names (name)
        return si[si[targetlevel].isin(levelnames)]['name'].values
    
    def calc_one_specific(self, targetlevel, levelnames, gids=None):
        """
        Args:
            targetlevel: a column name in sampleinfo dataframe
            g1: a set of names in the targetlevel column
        """
        dm = self.dms[targetlevel]['dm']
        mm = self.dms[targetlevel]['mm']
        v = self.dms[targetlevel]['v']
        if gids is not None:
            dm = dm.ix[gids]
            mm = mm.ix[gids]
            v = v.ix[gids]
        
        c = dm.major_axis
        i1 = c.isin(levelnames)
        i2 = ~i1
        n1 = N.sum(i1)
        n2 = N.sum(i2)
        print(n1,n2)
        mask1 = i1[:,N.newaxis]*i1[N.newaxis,:]+i2[:,N.newaxis]*i2[N.newaxis,:]
        mask2 = N.ones(len(i1)) - mask1
        sum1 = mask1.sum() - (n1+n2)
        sum2 = mask2.sum()
        dmv = dm.values
        dmv1 = dmv*mask1[N.newaxis,:,:]
        idx1 = dmv1<self.zth1
        idx2 = dmv1<self.zth2
        mmv = mm.values
        mmmask = mmv>=self.mmth
        idx1 = idx1*(~mmmask)
        idx2 = idx2*mmmask
        dmv1[idx1] = 0. # ignore those pairs within groups and min of the pair > minminth
        dmv1[idx2] = 0.
        dmv2 = dmv*mask2[N.newaxis,:,:] # don't apply zeroth to between groups
        sc1 = dmv1.sum(axis=2).sum(axis=1)
        sc2 = dmv2.sum(axis=2).sum(axis=1)
        df = PD.DataFrame({'_gidx':dm.items, 'sc1':sc1,'sc2':sc2,})
        df['score'] = (df['sc2']-df['sc1'])/(2.*n1*n2)
        df['di1'] = df['sc1']/float(sum1)
        df['di2'] = df['sc2']/float(sum2)
        df['didiff'] = df['di2']-df['di1']
        # calculate rd <=== should do this at sample level!
        #v1 = v.ix[df['_gidx'].values]
        #df['rd'] = (v1[levelnames]>self.rdth).mean(axis=1).values
        #df['gcov'] = v1[levelnames].mean(axis=1).values
        #cmpl = [x for x in v1.columns if x not in levelnames]
        #df['gcov2'] = v1[cmpl].mean(axis=1).values
        gc = self.gcov.ix[df['_gidx'].values] # align
        sn = self.get_snames(targetlevel, levelnames) # corresponding samples
        df['rd'] = (gc[sn]>self.rdth).mean(axis=1).values
        df['gcov'] = gc[sn].mean(axis=1).values
        cmpl = [x for x in gc.columns if x not in sn]
        df['gcov2'] = gc[cmpl].mean(axis=1).values
        
        df = df.sort_values('score',ascending=False)  
        
        self._locals = locals()
        return df
        
    
    def calc_many_specific(self, targetlevel, key2names, scoreth=None, rdratioth=0.6):
        """
        Args:
            targetlevel: name or cg1
            key2names: dict groupname (key) to names in targetlevel
            
        """
        dfs = []
        for k, ln in key2names.items():
            print('{0}...'.format(k))
            df = self.calc_one_specific(targetlevel, ln)
            cols = list(df.columns)
            if scoreth is not None:
                df = df[df['score']>scoreth].copy()
                print('scoreth{0}:{1}'.format(scoreth,len(df)))
            if rdratioth is not None:
                idx1 = (df['gcov']>df['gcov2'])&(df['rd']>rdratioth)
                idx2 = (df['gcov']<=df['gcov2'])&((1-df['rd'])>rdratioth)
                df = df[idx1|idx2].copy()
                print('rdratioth{0}:{1}'.format(rdratioth,len(df)))
            df['key'] = k
            df = df.sort_values('score',ascending=False)        
            df['rank'] = N.arange(len(df))
            df['id'] = df['key']+'.'+df['rank'].astype(str)
            dfs.append(df)
        df0 = PD.concat(dfs, ignore_index=True)
        g2cg1 = UT.df2dict(self.si, 'group', 'cg1')
        df0['region'] = [g2cg1.get(x,x) for x in df0['key']]
        df0 = df0[['region','key','id']+cols]
        df0 = self.annotate(df0)
        return df0

    def annotate(self, df):
        gbed = self.gbed
        gbedcols=['gen9_sym0','glocus','#uexons','#junc','pCSF','p60','gknown']
        for c in gbedcols:
            df[c] = gbed.ix[df['_gidx']][c].values
        df = df.rename(columns={'gen9_sym0':'symbol'}).set_index('_gidx')
        idx = df['symbol'].isnull()
        df.loc[idx, 'symbol'] = gbed.ix[df[idx].index]['gname']        
        return df
    
    def plot_dm(self, targetlevel, gid, key2pos, lvmax, title=None, fontsize=6,
                ymax=None, plotylim=True, g1=None, ax=None):
        if ax is None:
            fig, ax = P.subplots(1,1,figsize=(3,4))
        # make axes
        dmrect = (0.,0.25,1.,0.750) # main
        exrect = (0.,0.00,1.,0.145) # expression
        sbrect = (0.,0.15,1.,0.090) # legend
        args = dict(axisbg='w',frameon=True, xticks=[], yticks=[])
        ax_dm = PU.add_subplot_axes(ax,dmrect,args)
        ax_sb = PU.add_subplot_axes(ax,sbrect,args)
        args['frameon'] = False
        ax_ex = PU.add_subplot_axes(ax,exrect,args)
        
        dm1 = self.dms[targetlevel]['dm'].ix[gid]
        mm1 = self.dms[targetlevel]['mm'].ix[gid]
        v1 = self.dms[targetlevel]['v'].ix[gid]
        gc1 = self.gcov.ix[gid]
        
        cols = [x[1] for x in sorted([(key2pos[x],x) for x in dm1.columns])]
        dm1 = dm1[cols].ix[cols]
        mm1 = mm1[cols].ix[cols]
        v1 = v1.ix[cols]
        snames = [y for x in cols for y in self.get_snames(targetlevel, [x])]
        gc1 = gc1.ix[snames]
        
        # main
        dmv = dm1.values.copy()
        if g1 is not None:
            mmv = mm1.values
            c = dm1.columns
            i1 = c.isin(g1)
            i2 = ~i1
            mask = i1[:,N.newaxis]*i1[N.newaxis,:]+i2[:,N.newaxis]*i2[N.newaxis,:]
            th1 = dmv.max()*self.zth1
            th2 = dmv.max()*self.zth2
            idx1 = (mmv<self.mmth)&(dmv<th1)
            idx2 = (mmv>=self.mmth)&(dmv<th2)
            dmv[idx1*mask] = 0.
            dmv[idx1*mask] = 0.
            print(N.sum(i1),N.sum(i2))
        ax_dm.imshow(dmv, aspect='auto', interpolation='none', cmap='gray')
        # expression
        #ax_ex.plot(v1.values, '-', color='k', lw=0.3)
        ax_ex.plot(gc1.values, '-', color='k', lw=0.3)
        
        ymin = min(0, N.floor(N.min(v1)))
        if ymax is None:
            ymax = max(1, N.ceil(N.max(v1)))
        ax_ex.set_ylim((ymin,ymax))
        xmax = ax_ex.get_xlim()[1]
        if plotylim:
            ax_ex.text(xmax,ymax,'{:.0f}'.format(ymax), va='top',ha='left',fontsize=fontsize-2)
            ax_ex.text(xmax,ymin,'{:.0f}'.format(ymin), va='bottom',ha='left',fontsize=fontsize-2)
        # legend
        self._plot_legend(cols, key2pos, lvmax, ax_sb)
        # title
        if title is None:
            title='gid:{0}'.format(gid)
        ax.set_title(title,fontsize=fontsize)
        P.setp(ax, xticks=[], yticks=[], frame_on=False)        
        return ax

    def plot_all_gcov(self, df, targetlevel, key2pos, key2level, lvmax, 
                      figsize=(9,6), normrow=True, proportional=False):
        # df: output of calc_many_specific
        gcov = self.dms[targetlevel]['v'] # log2(ug1kn+1).groupby(...)
        cols = [x[1] for x in sorted([(key2pos[x],x) for x in gcov.columns])]
        gcov = gcov[cols]
        # make sure highest score is selected
        if df.index.name == '_gidx':
            df = df.reset_index()
        df = df.sort_values(['_gidx','score'],ascending=False)
        df = df.groupby('_gidx', sort=False).first()        
        # set pos & level
        df['_pos_'] = [key2pos[x] for x in df['key']]
        df['_lvl_'] = [key2level[x] for x in df['key']]

        # make axes
        fig,ax = P.subplots(1,1,figsize=figsize)
        args = dict(axisbg='w', frameon=True, xticks=[], yticks=[])
        axl = PU.add_subplot_axes(ax, [0,0,0.9,0.07],args) # legends at the bottom
        if proportional:
            nsub = N.array([N.sum(df['_lvl_']==3-i) for i in range(3)])
            hs = 0.9*nsub/N.sum(nsub)
            ys = [0.08, 0.08+0.01+hs[0], 0.08+0.02+hs[0]+hs[1]]
        else:
            hs = [0.3 for i in range(3)]
            ys = [0.08+0.31*i for i in range(3)]

        axm = [PU.add_subplot_axes(ax, [0.00,ys[i],0.90,hs[i]],args) for i in range(3)]
        axu = [PU.add_subplot_axes(ax, [0.91,ys[i],0.03,hs[i]],args) for i in range(3)]
        axc = [PU.add_subplot_axes(ax, [0.96,ys[i],0.03,hs[i]/2.],args) for i in range(3)]

        
        def _plot_main(sub,gcov,ax):
            gids = [x for x in sub.index.values]
            m = gcov.ix[gids].values
            mmin = m.min(axis=1)
            mmax = m.max(axis=1)
            #mn = (m-mmin[:,N.newaxis])/(mmax-mmin)[:,N.newaxis]
            if normrow:
                mn = m/mmax[:,N.newaxis]
                zmin,zmax = 0,1
            else:
                mn = m
                zmin = 0 #mmin.min()
                zmax = mmax.max()
            ax.imshow(mn, aspect='auto', interpolation='nearest', cmap='gray_r')
            ax.text(0,0,str(len(sub)),horizontalalignment='right', transform=ax.transAxes)
            #ax.text(0,1,'0',horizontalalignment='right', transform=ax.transAxes)
            return zmin,zmax
            
        def _plot_ukbar(sub, ax):
            uk = N.array([[x.startswith('u.') for x in sub['gknown']]]).T
            ax.imshow(uk, aspect='auto', interpolation='nearest', cmap='Reds')
            
        def _plot_cbar(zmin,zmax,ax):
            cb = N.array([N.arange(0,zmax,zmax/32.)[::-1]]).T
            ax.imshow(cb,aspect='auto', interpolation='nearest', cmap='gray_r', vmax=zmax)
            ax.text(1,32,str(zmin))
            ax.text(1,0,str(zmax))
            
        
        # main&known/unknwon&colorbar
        for i in range(3):
            sub = df[df['_lvl_']==3-i].sort_values('_pos_')
            zmin,zmax = _plot_main(sub, gcov, axm[i])
            _plot_ukbar(sub, axu[i])
            if normrow:
                if i==0:
                    _plot_cbar(zmin, zmax, axc[i])
                else:
                    P.setp(axc[i], xticks=[], yticks=[], frame_on=False)
            else:
                _plot_cbar(zmin, zmax, axc[i])
        # legends
        self._plot_legend(cols, key2pos, lvmax, axl)
        
        P.setp(ax, xticks=[], yticks=[], frame_on=False)
        return fig
    
    def _plot_legend(self, cols, key2pos, lvmax, ax):
        si = self.si
        h = 1./3
        args = dict(axisbg='w', frameon=False, xticks=[], yticks=[])
        #args['frameon'] = False
        ax_sb = [PU.add_subplot_axes(ax, [0,h*i,1,h],args) for i in range(3)]
        #flds = ['gtrans3','region','type']
        clrs = ['jet','nipy_spectral','gist_rainbow']
        #for i,(a,fld,cm) in enumerate(zip(ax_sb,flds,clrs)):
        for i,(a,cm) in enumerate(zip(ax_sb,clrs)):
            m = N.array([[key2pos[x][2-i] for x in cols]])
            vmax = lvmax[2-i]
            a.imshow(m, aspect='auto', interpolation='none', cmap=cm, vmin=0, vmax=vmax)
    
    
    
    
    
CLS = dict(gtrans3='jet',region='nipy_spectral',type='gist_rainbow')
FLDS = ['gtrans3','region','type']
def make_legend_bars(orders,flds=FLDS, cls=CLS):
    nfld = len(flds)
    nums = N.array([len(orders[x]) for x in flds])
    tot = N.sum(nums)
    spaces = 0.01*(nfld-1)
    dy = (1-spaces)/tot
    hs = [dy*x for x in nums]
    ys = [0.01*i+N.sum(hs[:i]) for i in range(nfld)]
    
    fig, ax = P.subplots(1,1,figsize=(3,0.25*tot+0.1*(nfld-1)))
    args = dict(axisbg='w',frameon=True, xticks=[], yticks=[])
    for i in range(nfld):
        axi = PU.add_subplot_axes(ax, [0,ys[i],0.15, hs[i]], args)
        m = N.array([range(nums[i])]).T
        axi.imshow(m, aspect='auto', interpolation='nearest', cmap=cls[flds[i]], vmin=0, vmax=nums[i])
        axi.set_yticks(range(nums[i]))
        axi.set_yticklabels(orders[flds[i]])
        axi.tick_params(axis='both', which='both',length=0)
        
    P.setp(ax, frame_on=False, xticks=[], yticks=[])
    return fig