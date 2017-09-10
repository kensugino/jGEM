# => moved to ponder/plotgene.py

import pandas as PD
import numpy as N
import matplotlib.pyplot as P

#import bedtools as BT
#import gtfgffbed as GGB
import jgem.utils as UT

import random
from matplotlib.patches import Rectangle, PathPatch, Circle, Ellipse
from matplotlib.lines import Line2D
from matplotlib.path import Path
import matplotlib.colors as C
import matplotlib.cm as CM

class SpliceFig(object):
    
    
    def __init__(self, ex, sj, xmargin=None, ymargin=0.25, compress=True, ecov='ecov', 
            ucnt='ucnt',mcnt='mcnt',minlw=1,drawscalebar=True,ecovth=None,jcntth=None,
            origin=None, sortexby=None,fontsize=7):
        self.ymargin = ymargin
        self.ecov = ecov
        self.ucnt = ucnt
        self.jcnt = jcnt = 'jcnt'
        self.mcnt = mcnt
        self.minlw = minlw
        self.drawscalebar = drawscalebar
        self.ecovth = ecovth
        self.jcntth = jcntth
        self.ex = ex = ex.copy()
        self.sj = sj = sj.copy()
        self.compress = compress
        self.fontsize=fontsize
        if sortexby is None:
            self.sortexby = ecov
        else:
            self.sortexby = sortexby # when plotting multiple and comparing, you want to use same sorting
        # start and end, strand
        
        if ex.iloc[0]['strand']=='+':
            if origin is None:
                origin = ex['st'].min()
            ex['xst'] = ex['st']-origin
            ex['xed'] = ex['ed']-origin
            self.strand='+'
            self.origin=origin
        else:
            if origin is None:
                origin = ex['ed'].max()
            ex['xst'] = origin - ex['ed']
            ex['xed'] = origin - ex['st']
            self.strand='-'
            self.origin=origin
        # fix old a_id null 
        if (ex['a_id'].min()==-1) and (N.sum(ex['a_id']==0)==0):
            ex.loc[ex['a_id']==-1, 'a_id'] = 0
            ex.loc[ex['d_id']==-1, 'd_id'] = 0
            sj.loc[sj['a_id']==-1, 'a_id'] = 0
            sj.loc[sj['d_id']==-1, 'd_id'] = 0
        
        ex['len'] = ex['xed'] - ex['xst']
        if xmargin is None:
            xmargin = int(ex['len'].mean())
        self.xmargin = xmargin

        if ecov not in ex.columns:
            ex[ecov] = 1
        if (ucnt not in sj.columns) or (mcnt not in sj.columns):
            sj[jcnt] = 1
            sj[jcnt+'_ls'] = 'solid'
        else:
            # sj uniq, mult
            sj[jcnt] = [x or y for x,y in sj[[ucnt,mcnt]].values]
            sj[jcnt+'_ls'] = ['solid' if x else 'dashed' for x in sj[ucnt]]

        if ecovth is not None:
            self.ex = ex = ex[ex[ecov]>ecovth].copy()
        if jcntth is not None:
            self.sj = sj = sj[sj[jcnt]>jcntth].copy()
        if len(ex)==0:
            return

        # find exon groups
        if 'asize' not in ex.columns:
            a2size = dict(UT.izipcols(ex.groupby('a_id').size().reset_index(), ['a_id',0]))
            d2size = dict(UT.izipcols(ex.groupby('d_id').size().reset_index(), ['d_id',0]))
            a2size[0]=0
            d2size[0]=0
            ex['asize'] = [a2size[x] for x in ex['a_id']]
            ex['dsize'] = [d2size[x] for x in ex['d_id']]
        ex['group'] = ['a{0}'.format(ai) if (a!=0 and a>d) else'd{0}'.format(di) for a,ai,d,di in ex[['asize','a_id','dsize','d_id']].values]
        # find exon group st, ed
        exg = ex.groupby('group')
        g2st = dict(UT.izipcols(exg['xst'].min().reset_index(), ['group','xst']))
        g2ed = dict(UT.izipcols(exg['xed'].max().reset_index(), ['group','xed']))
        g2size = dict(UT.izipcols(exg.size().reset_index(), ['group',0]))
        ex['gst'] = [g2st[x] for x in ex['group']]
        ex['ged'] = [g2ed[x] for x in ex['group']]
        ex['gsize'] = [g2size[x] for x in ex['group']]
        #self.ex = ex = ex.sort_values(['group',ecov]) #'gst','ged','xst','xed'])
        self.ex = ex = ex.sort_values(['group',self.sortexby]) #'gst','ged','xst','xed'])
        # find exon y pos within group
        def _eypos(gs):
            g0,s0 = gs[0] # first g
            cnt = 0
            yield cnt - (s0-1)/2.
            for g1,s1 in gs[1:]:
                if g1==g0:
                    cnt +=1
                else:
                    cnt = 0
                yield cnt - (s1-1)/2.
                g0 = g1
        ex['eypos'] = [x for x in _eypos(ex[['group','gsize']].values)]
        # find group y center pos
        self.gr = gr = ex.groupby('group')[['gst','ged','gsize']].first().sort_values(['gst','ged'])
        gr['len'] = gr['ged']-gr['gst']
        def _gypos(gr):
            side = 1
            r0 = gr.iloc[0]
            h = r0['gsize']/2.
            ged0 = r0['ged']
            gy0 = {1:h, -1:-h} # remember filled height both side (1,-1)
            yield 0 # first one gets center
            for gst1,ged1,gsiz1 in gr[['gst','ged','gsize']].values[1:]:
                h = gsiz1/2.
                if ged0<=gst1: # no overlap
                    gy0 = {1:h, -1:-h}
                    yield 0
                else:
                    gy1 = gy0[side] + side*gsiz1/2.
                    gy0[side] = gy0[side] + side*gsiz1
                    side = -1*side # flip side
                    yield gy1
                gst0 = gst1
                ged0 = max(ged0, ged1)
        gr['gypos'] = [x for x in _gypos(gr)]
        # compress x coord
        if compress:
            def _gxst(gr):
                r0 = gr.iloc[0]
                delta = 0
                yield r0['gst'] - delta # 0
                ged0 = r0['ged']
                for i, r1 in gr.iloc[1:].iterrows():
                    gst1 = r1['gst']
                    if gst1-ged0>self.xmargin:
                        delta += (gst1-ged0-self.xmargin)
                    yield gst1 - delta
                    ged0 = r1['ged']
            gr['cst'] = [x for x in _gxst(gr)]
        else:
            gr['cst'] = gr['gst']
        #gr['ced'] = gr['cst']+gr['len']
        ex['cst0'] = [gr['cst'].ix[g]+(xst-gst) for g,xst,gst in ex[['group','xst','gst']].values]
        ex['ced0'] = ex['cst0']+ex['len']
        if self.strand=='+':
            ex['cst'] = origin + ex['cst0']
            ex['ced'] = origin + ex['ced0']
        else:
            ex['cst'] = origin - ex['ced0']
            ex['ced'] = origin - ex['cst0']
        ex['ey'] = [ey+gr['gypos'].ix[g] for ey,g in ex[['eypos','group']].values]

    def draw(self, ax=None, cm='R', cm2='G', ecov2=None, jcnt2=None, xlim=None):
        if len(self.ex)==0:
            return
        if ax is None:
            fig,ax = P.subplots(1,1,figsize=(12,5))
        self.draw_junctions(ax, cm=cm, jcnt=self.jcnt, cm2=cm2, jcnt2=jcnt2, ecov=self.ecov, ecov2=ecov2,xlim=xlim)
        self.draw_exons(ax, cm=cm, ecov=self.ecov, cm2=cm2, ecov2=ecov2, xlim=xlim)
        if self.drawscalebar:
            self.draw_bar(ax)

    def draw_bar(self, ax):
        # 1kb scale bar at right top
        ylim = ax.get_ylim()
        xlim = ax.get_xlim()
        yd = ylim[1]-ylim[0]
        xd = xlim[1]-xlim[0]
        ypos = ylim[0]+yd*0.8
        ypos2 = ylim[0]+yd*0.9
        xpos = xlim[0]+xd*0.9
        xpos1 = max(xlim[0], xpos-1000)
        xpos2 = max(xlim[0], xpos-500)
        ax.plot([xpos1, xpos], [ypos, ypos], 'k-', lw=2)
        suf = '({0})'.format(self.strand)
        if self.compress:
            suf = '({0},c)'.format(self.strand) #'(intron compressed)'
        if xpos1==xpos-1000:
            ax.text(xpos2,ypos2,'1kb{0}'.format(suf),fontsize=self.fontsize)
        else:
            ax.text(xpos2,ypos2,'{0:d}bp{1}'.format(xpos-xpos1,suf),fontsize=self.fontsize)
        
    def draw_exons(self, ax, cm='Reds', ecov='ecov', cm2='Greens', ecov2=None, xlim=None):
        ex = self.ex
        hh = 0.5-self.ymargin
        emax = N.log2(ex[ecov].max()+1)
        #print 'emax=', ex[ecov].max()
        sm = self.get_scalarmap(emax,mapname=cm)
        if ecov2 is None:
            for x0,w,y,ec in ex[['cst','len','ey',ecov]].values:
                #print 'ec=',ec
                ec = N.log2(ec+1)
                fc = sm.to_rgba(ec)
                #eclr = sm2.to_rgba(ec)
                lw = max(0.3, ec/emax)
                ax.add_patch(Rectangle((x0,y-hh),w,2*hh,facecolor=fc,lw=lw,alpha=0.9))#edgecolor=ec))
        else:
            emax2 = N.log2(ex[ecov2].max()+1)
            sm2 = self.get_scalarmap(emax2,mapname=cm2)
            for x0,w,y,ec,ec2 in ex[['cst','len','ey',ecov,ecov2]].values:
                ec = N.log2(ec+1)
                ec2 = N.log2(ec2+1)
                fc = sm.to_rgba(ec)
                fc2 = sm2.to_rgba(ec2)
                fc3 = self._addcolor(fc,fc2)
                lw = N.max([0.3, ec/emax, ec2/emax2])
                ax.add_patch(Rectangle((x0,y-hh),w,2*hh,facecolor=fc3,lw=lw,alpha=0.9))#edgecolor=ec))

        xmi = ex['cst'].min()
        xma = ex['ced'].max()
        ymi = ex['ey'].min()
        yma = ex['ey'].max()
        if xlim is None:
            ax.set_xlim([xmi-self.xmargin, xma+self.xmargin])
        else:
            ax.set_xlim(xlim)
        d = 2+self.ymargin
        ax.set_ylim([ymi-d, yma+d])

    def _addcolor(self,fc,fc2):
        #return tuple(map(lambda x: min(x, 1.), N.array(fc)+N.array(fc2)))
        #return tuple([max(x,y) for x,y in zip(fc,fc2)])
        return tuple([min(x,y) for x,y in zip(fc,fc2)])
        #return tuple((N.array(fc)+N.array(fc2))/2.)
                        
    def draw_junctions(self, ax, cm='Reds', jcnt='jcnt', cm2='Greens', jcnt2=None, ecov='ecov', ecov2=None, xlim=None):
        ex = self.ex
        sj = self.sj.sort_values(['d_id',jcnt],ascending=False)
        if len(sj)==0:
            return
        maxjcnt = N.log2(sj[jcnt].max()+1)
        print('sjmax:{0},sjmaxidx:{1},locus: {2}:{3}-{4}'.format(sj[jcnt].max(), sj[jcnt].argmax(), ex['chr'].iloc[0], ex['st'].min(), ex['ed'].max()))
        sm = self.get_scalarmap(maxjcnt,mapname=cm)
        lwmax = 2
        if xlim is None:
            xw = ex['ced'].max()-ex['cst'].min()
        else:
            xw = xlim[1]-xlim[0]
        yw = ex['ey'].max()-ex['ey'].min()
        # cw = 0.01*xw
        # ch = 0.1*(1-2*self.ymargin)
        minlw = self.minlw
        self._side=1
        print('xw=',xw, 'xlim=', xlim, 'yw=', yw)
        def _draw(x0,y0,x1,y1,ec,lw,ls):
            global dcnt
            if (len(ex)>2):
                pts = self.bezierpts(x0,y0,x1,y1,xw,yw)
                a = Path(pts, [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4])
                ax.add_patch(PathPatch(a, fc='none', alpha=0.8, lw=lw, ec=ec, linestyle=ls))
            else:
                ax.add_line(Line2D((x0,x1),(y0,y1),alpha=0.9, lw=lw, color=ec,linestyle=ls))
            # place dot in acceptor/donor
            #ax.add_patch(Ellipse((x0,y0), width=cw, height=ch, facecolor='k'))
            #ax.add_patch(Ellipse((x1,y1), width=cw, height=ch, facecolor='k'))
            #ax.plot([x0],[y0],'k.',ms=1.5)
            #ax.plot([x1],[y1],'k.',ms=1.5)

        d0 = sj.iloc[0]['d_id']
        afld = 'cst' if self.strand=='+' else 'ced'
        dfld = 'ced' if self.strand=='+' else 'cst'
        #print('afld:{0},dfld:{1},strand:{2}'.format(afld,dfld,self.strand))

        if jcnt2 is None:
            for ai,di,jc,ls in sj[['a_id','d_id',jcnt,jcnt+'_ls']].values:
                if di!=d0:
                    self._side=1 # initialize
                d0 = di
                # connect donor to acceptor
                jc = N.log2(jc+1)
                ea = ex[ex['a_id']==ai]
                ed = ex[ex['d_id']==di]
                c1m = N.log2(ea[ecov].max()+1)
                c0m = N.log2(ed[ecov].max()+1)
                for x1,y1,c1,i1 in ea[[afld,'ey',ecov,'_id']].values:
                    c1 = N.log2(c1+1)
                    for x0,y0,c0,i0 in ed[[dfld,'ey',ecov,'_id']].values:
                        c0 = N.log2(c0+1)
                        jc1 = jc*max(0.3, (c1/c1m)*(c0/c0m))
                        if jc1==0:
                            ec = sm.to_rgba(0.7*maxjcnt)
                            _draw(x0,y0,x1,y1,ec,2,'dotted')
                        else:
                            ec = sm.to_rgba(jc1)
                            lw = max(minlw,lwmax*jc1/maxjcnt)
                            _draw(x0,y0,x1,y1,ec,lw,ls)
                        if y0==0 and y1==0:
                            self._side *= -1 # flip
        else:
            maxjcnt2 = N.log2(sj[jcnt2].max()+1)
            sm2 = self.get_scalarmap(maxjcnt2,mapname=cm2)
            for ai,di,jc,jc2,ls in sj[['a_id','d_id',jcnt,jcnt2,jcnt+'_ls']].values:
                if di!=d0:
                    self._side=1
                d0=di
                # connect donor to acceptor
                jc = N.log2(jc+1)
                jc2 = N.log2(jc2+1)
                ea = ex[ex['a_id']==ai]
                ed = ex[ex['d_id']==di]
                c1m = N.log2(ea[ecov].max()+1)
                c0m = N.log2(ed[ecov].max()+1)
                c1m2 = N.log2(ea[ecov2].max()+1)
                c0m2 = N.log2(ed[ecov2].max()+1)
                for x1,y1,c1,c12 in ea[[afld,'ey',ecov,ecov2]].values:
                    c1 = N.log2(c1+1)
                    c12 = N.log2(c12+1)
                    for x0,y0,c0,c02 in ed[[dfld,'ey',ecov,ecov2]].values:
                        c0 = N.log2(c0+1)
                        c02 = N.log2(c02+1)
                        jc1 = jc*max(0.3, (c1/c1m)*(c0/c0m))
                        jc21 = jc2*max(0.3, (c12/c1m2)*(c02/c0m2))
                        lw = max(minlw, lwmax*N.max([jc1/maxjcnt,jc21/maxjcnt2]))
                        if (jc1==0) and (jc21==0):
                            ec = sm.to_rgba(0.5*maxjcnt)
                            ec2 = sm2.to_rgba(0.5*maxjcnt2)
                            ec3 = self._addcolor(ec,ec2)
                            _draw(x0,y0,x1,y1,ec3,2,'dotted')
                        else:
                            ec = sm.to_rgba(jc1)
                            ec2 = sm2.to_rgba(jc21)
                            ec3 = self._addcolor(ec,ec2)
                            _draw(x0,y0,x1,y1,ec3,lw,ls)
                        if y0==0 and y1==0:
                            self._side *=-1
                        
                                    
    def get_scalarmap(self,vmax,vmin=0, mapname='Reds'):
        return Colors(mapname, vmax, vmin)
                        
    def bezierpts(self,x0,y0,x1,y1,xw,yw):
        width = x1-x0
        off = 0.1
        cx0 = x0+width*off
        cx1 = x0+width*(1-off)
        coef = 1
        coef2 = (width/float(xw)-0.1)*4
        if y0==0:
            if y1==0:
                d = self._side
            elif y1<0:
                d = -1
            else:
                d = 1
        elif (y0>0):
            d = 1
        else:
            d = -1
        cy0 = y0+d*coef
        cy1 = y1+d*coef
        if y0<y1:
            if d==1:
                cy0 += (y1-y0)*coef2
            else:
                cy1 -= (y1-y0)*coef2
        elif y0>y1:
            if d==1:
                cy1 += (y0-y1)*coef2
            else:
                cy0 -= (y0-y1)*coef2
        return [(x0,y0),(cx0,cy0),(cx1,cy1),(x1,y1)]

        
        
class Colors(object):
    
    def __init__(self, mapname, vmax, vmin=0, nl=32):
        self.mn = mapname
        self.vmin = vmin
        self.vmax = vmax
        self.d = d = 1./nl
        if mapname=='C':
            self.rgba = [(1.-x,1.,1.,1.) for x in N.arange(0,1+d,d)]
        elif mapname=='M':
            self.rgba = [(1.,1.-x,1.,1.) for x in N.arange(0,1+d,d)]
        elif mapname=='Y':
            self.rgba = [(1.,1.,1.-x,1.) for x in N.arange(0,1+d,d)]
        elif mapname=='R':
            self.rgba = [(1.,1.-x,1.-x,1.) for x in N.arange(0,1+d,d)]
        elif mapname=='G':
            self.rgba = [(1.-x,1.,1.-x,1.) for x in N.arange(0,1+d,d)]
        elif mapname=='B':
            self.rgba = [(1.-x,1.-x,1.,1.) for x in N.arange(0,1+d,d)]
        else:
            cm = P.get_cmap(mapname)
            cnorm = C.Normalize(vmin=0,vmax=1.)
            self.sm = sm = CM.ScalarMappable(norm=cnorm,cmap=cm)
            self.rgba = [sm.to_rgba(x) for x in N.arange(0,1+d,d)]
            
    def to_rgba(self, v):
        d = self.d
        if self.mn in ['R','G','B','C','M','Y']:
            vn = max(0., (v-self.vmin)/self.vmax)
            vn = min(1., vn)
            vni = int(vn/d)
            return self.rgba[vni]
        return self.sm.to_rgba(v)
        
        
