"""
    GO Overrepresentation Analysis
    
"""

from random import Random
R = Random()

import sys
import os
import numpy as N
import pandas as PD
from matplotlib import pylab as P
from time import time
from scipy.stats import chi2

from bioinfo3.go.graphviz import Graphviz
from scipy.stats import hypergeom
from rpy2.robjects import r, FloatVector

pmf = hypergeom.pmf
def hyperg(x,n1,n2,num):
    return pmf(x,n1+n2,n1,num)
def rhyperg(x,n1,n2,num):    
    return r.dhyper(x,n1,n2,num)[0]
def rpchisq(x,df):    
    return r.pchisq(x, df)[0]

"""
# BUG in scipy.stats.hypergeom?
from scipy import stats
>>> stats.hypergeom.pmf(13,12854,581,149)
# array(nan)  # was this years ago
0.0099366208035666682  # looks like it's fixed (as of 2014-01-02)

from rpy2.robjects import r
>>> r.dhyper(13,581,12854-581,149)[0]
# 0.0099366208036166248  # years ago
0.009936620803616614  # as of 2014-01-02 

import timeit
def compare_scipy_rpy():
    t0 = timeit.Timer(stmt='hypergeom.pmf([13]*1000,12854,581,149)',
                      setup='from scipy.stats import hypergeom')
    t1 = timeit.Timer(stmt='r.dhyper(numpy.array([13]*1000),581,12854-581,149)',
                      setup='from rpy2.robjects import r;import numpy')
    rep = 1
    while t0.timeit(rep) < .25:
        rep *=10
    print('Going for %s reps' % rep)
    print('scipy.stats.hypergeom.pmf:',  t0.repeat(3, rep))
    print('r.dhyper:', t1.repeat(3,rep))

# as of 2014-01-02 
# Going for 1000 reps
# scipy.stats.hypergeom.pmf: [0.25429932100814767, 0.25074021398904733, 0.24888762301998213]
# r.dhyper: [0.16580171600799076, 0.15831547800917178, 0.15673348799464293]

# final winner is GSL with Cython : about 10~11 times faster than rpy    

#     
"""

class EnrichmentAnalysis(object):

    def __init__(self, population_ids, id2term):
        self.pop_ids = sorted(list(set(population_ids)))
        self.id2term = id2term
        # pop_id => term
        self.i2t = {k: set(id2term[k]) for k in population_ids if k in id2term} 
        t2i = {}
        for i in self.i2t:
            for t in self.i2t[i]:
                t2i.setdefault(t,set()).add(i)
        self.t2i = t2i

    def calculate(self, selected_ids):
        self.sel_ids = sorted(list(set(selected_ids)))
        
        i2t = self.i2t
        t2i = self.t2i
        tot = len(i2t)  # population total (with term association)
        s2t = {s:i2t[s] for s in self.sel_ids if s in i2t}  # id => term in selected
        num = len(s2t)
        t2s = {}
        for s in s2t:
            for t in s2t[s]:
                t2s.setdefault(t,set()).add(s)
        terms = sorted(list(t2s.keys()))
        er = EnrichmentResult(tot, num, terms)
        df = er.df
        p = float(num)/float(tot)
        for t in terms:
            if t not in t2s:
                continue
            n0 = len(t2s[t]) # number of selected for term t
            n1 = len(t2i[t]) # number of term t for population
            pval = N.sum(pmf(range(n0,num+1),tot,n1,num))
            odds = (float(n0)/float(n1))/p
            ids = ','.join(t2s[t])
            df.loc[t] = [n1,float(n1)/tot,n0,float(n0)/num,odds,pval,ids]
        return er


class EnrichmentResult(object):

    def __init__(self, pop_total, sel_total, terms):
        self.pop_total = pop_total
        self.sel_total = sel_total
        index = PD.Index(terms)  # term id as index
        columns = ['total', 'total %' ,'selected', 'selected %', 'odds', 'p.values', 'ids']
        self.df = PD.DataFrame({}, index=index, columns=columns)

    def __repr__(self):
        lines = 'Population: %d\n' % self.pop_total
        lines += 'Selected: %d\n' % self.sel_total
        lines += self.df.__str__()
        return lines


class GOSlimEnrichment(object):

    FLDS = ['selected %', 'odds', 'p.values']

    def __init__(self, populations, id2slim):
        self.populations= populations # dict
        self.id2slim = id2slim # ['gob','gom','goc'] => id2term

    def calculate(self, selected, conditions):
        e = self.populations
        s = self.selected = selected
        self.conditions = conditions
        A = EnrichmentAnalysis
        which = ['gob','gom','goc']
        i2s = self.id2slim
        self.rslts = {c: {w: A(e[c],i2s[w]).calculate(s[c]) for w in which} for c in conditions}

    def plot_goslim_bar(self, conds, flds = ['odds', '-log10(p)','selected %'],gos = ['goc','gom','gob'], **kw):
        # flds 'selected %', 
        # some prep
        xlims = {'selected %':(0,100),'odds':(0,3),'-log10(p)':(0,150)}
        xticks = {'selected %':[0,50,100],'odds':[0,2,4],'-log10(p)':[0,2,4]}
        th = {'odds':2,'-log10(p)':-N.log10(0.05)}
        titles = {}#{'selected %':'%','odds':'odds','p.values':'-log10(p)'}
        # merge
        rslts = self.rslts
        cps = {}
        dfs = {w:{f:PD.DataFrame({c:rslts[c][w].df[f] for c in conds},columns=conds, dtype='float') for f in self.FLDS} for w in gos}
        for w in gos:
            dfs[w]['selected %'] = 100*(dfs[w]['selected %'])
            #x = (-2*N.log(dfs[w]['p.values'])).sum(axis=1)
            #cp = PD.Series(1. - chi2.cdf(x, 2*len(conds)), index=dfs[w]['p.values'].index)
            #cp = (-N.log10(dfs[w]['p.values'])).mean(axis=1)
            #cp = -N.log10((dfs[w]['p.values']).mean(axis=1))
            #print cp
            #cp.sort(ascending=False)
            #cps[w] = cp
            dfs[w]['-log10(p)'] = -N.log10(dfs[w]['p.values'])
            tmp = rslts[conds[0]][w].df
            idx = [(not x.startswith('other')) and (y>10) for x,y in zip(tmp.index.values, tmp['total'].values)]
            for f in self.FLDS+['-log10(p)']: # sort according to cp value
                dfs[w][f] = dfs[w][f][idx]
            cp = dfs[w]['odds']
            cp = cp.sort(ascending=False)
            for f in self.FLDS+['-log10(p)']: # sort according to cp value
                dfs[w][f] = dfs[w][f].loc[cp.index]
            print(dfs[w])
        #self.cps = cps
        # calculate sizes in pts
        labelw = kw.get('labelw', 250.)
        margin = kw.get('margin', 25.)
        spacer = kw.get('spacer', 10.)
        panelw = kw.get('panelw', 100.)
        lineh = kw.get('lineh', 20.)
        panelh = {w: len(dfs[w][flds[0]])*lineh for w in gos}
        figw = float(margin+labelw+(panelw+spacer)*len(flds)+margin)
        figh = float(margin+sum(panelh.values())+spacer*len(gos)+margin)
        dpi=80
        fig = P.figure(figsize=(figw/dpi,figh/dpi))
        #dpi = float(fig.get_dpi())
        #fig.set_size_inches(figw/dpi,figh/dpi)
        bottom = margin/figh
        width = panelw/figw
        for i in range(len(gos)):
            left = (margin+labelw)/figw
            height = panelh[gos[i]]/figh
            w = gos[i]
            for j in range(len(flds)):
                f = flds[j]
                ax = fig.add_axes([left,bottom,width,height])
                df = dfs[w][f]
                df.plot(kind='barh', ax=ax, legend=False, grid=False)
                ax.set_xlim(xlims[f])
                ax.set_xticks(xticks[f])
                ax.yaxis.set_tick_params(which='both', width=0, size=0)
                if f in th:
                    ax.plot([th[f],th[f]], ax.get_ylim(), 'r--')
                if j==0 and i==0:
                    patches, labels = ax.get_legend_handles_labels()
                    P.legend(patches, labels, loc='best', prop={'size':8})
                if j>0:
                    ax.set_yticklabels([])
                if i>0:
                    ax.set_xticklabels([])
                else:
                    for t in ax.xaxis.get_ticklabels():
                        t.set_fontsize(8)
                if i == len(gos)-1:
                    ax.set_title(titles.get(f,f))
                    ax.title.set_fontsize(10)
                left += (panelw+spacer)/figw
            bottom += (panelh[gos[i]]+spacer)/figh  
        fig.show()
        return fig


    def plot_goslim_mat(self, conds, xlabels=None, flds = ['odds', '-log10(p)','cp'],gos = ['goc','gom','gob'], **kw):
        if xlabels is None:
            xlabels = conds
        # flds 'selected %'
        # some prep
        # xlims = {'selected %':(0,100),'odds':(0,4),'-log10(p)':(0,4),'cp':(0,5)}
        # xticks = {'selected %':[0,50,100],'odds':[0,2,4],'-log10(p)':[0,2,4],'cp':[0,5]}
        xlims = {'selected %':(0,100),'odds':(0,4),'-log10(p)':(0,4),'cp':(0,21)}
        xticks = {'selected %':[0,50,100],'odds':[0,2,4],'-log10(p)':[0,2,4],'cp':[0,10,20]}
        th = {'odds':2,'-log10(p)':-N.log10(0.05)}
        zlims = {'selected %':(0,1),'odds':(0,3),'-log10(p)':(0,3),'cp':(0,10)}
        #titles = {'cp':'combined p', 'selected %':'%'}#,'odds':'odds','-log10(p)':'-log10(p)'}
        titles = {'cp':'-log10(cp)*odds', 'selected %':'%'}#,'odds':'odds','-log10(p)':'-log10(p)'}
        modes = {'selected %':'imshow','odds':'imshow','-log10(p)':'imshow','cp':'barh'}
        # merge
        rslts = self.rslts
        dfs = {w:{f:PD.DataFrame({c:rslts[c][w].df[f] for c in conds},columns=conds, dtype='float') for f in self.FLDS} for w in gos}
        for w in gos:
            x = (-2*N.log(dfs[w]['p.values'])).sum(axis=1)
            cp = -N.log10(PD.Series(1. - chi2.cdf(x, 2*len(conds)), index=dfs[w]['p.values'].index))
            #cp = (-N.log10(dfs[w]['p.values'])).mean(axis=1)
            #print(w,cp)
            p = dfs[w]['p.values']
            o = dfs[w]['odds']
            #pmean = -N.log10(p.mean(axis=1))
            #pmean = (-N.log10(p)).mean(axis=1)
            omean = o.mean(axis=1)
            #omean[o.isnull().sum(axis=1)>0] = 0.
            #cp = -N.log10(pmean) * omean
            #cp = pmean*omean
            cp = cp*omean
            #cp = ((-N.log10(p))*o).mean(axis=1)
            print(w, cp)
            cp.sort(ascending=True)
            cp[N.isinf(cp)] = 100.
            dfs[w]['cp'] = cp
            dfs[w]['-log10(p)'] = -N.log10(dfs[w]['p.values'])
            for f in self.FLDS+['-log10(p)']:
                dfs[w][f] = dfs[w][f].loc[cp.index]
        self.dfs = dfs
        # calculate sizes in pts
        labelh = kw.get('labelh',50.)
        labelw = kw.get('labelw', 250.)
        margin = kw.get('margin', 30.)
        spacer = kw.get('spacer', 10.)
        #panelw = kw.get('panelw', 100.)
        colorbarh = kw.get('colorbarh',10.)
        lineh = kw.get('lineh', 15.)
        panelw = lineh*len(conds)
        cm = getattr(P.cm, kw.get('cm', 'gist_yarg'))
        panelh = {w: len(dfs[w][flds[0]])*lineh for w in gos}
        figw = float(margin+labelw+(panelw+spacer)*len(flds)+margin)
        figh = float(margin+labelh+sum(panelh.values())+spacer*len(gos)+margin+colorbarh+2*spacer)
        dpi=80
        fig = P.figure(figsize=(figw/dpi,figh/dpi))
        #dpi = float(fig.get_dpi())
        #fig.set_size_inches(figw/dpi,figh/dpi)
        bottom = (margin+labelh)/figh
        width = panelw/figw
        cbh = colorbarh/figh
        for i in range(len(gos)):
            left = (margin+labelw)/figw
            height = panelh[gos[i]]/figh
            w = gos[i]
            for j in range(len(flds)):
                f = flds[j]
                ax = fig.add_axes([left,bottom,width,height])
                df = dfs[w][f]
                vmin, vmax = zlims[f]
                if modes[f]=='imshow':
                    im = ax.imshow(df.values, cmap=cm, aspect='auto', vmin=vmin, vmax=vmax, origin='lower',interpolation='nearest')
                    ax.set_xticks(N.arange(len(conds)))
                    #imshow does not work on saving vector (pdf) if interpolation set to 'none'
                    #im = ax.pcolormesh(df.values, cmap=cm, vmin=vmin, vmax=vmax, shading='flat')
                    # ax.set_xticks(N.arange(len(conds))+0.5)
                    ax.yaxis.set_tick_params(which='both', width=0, size=0)
                    ax.xaxis.set_tick_params(which='both', width=0, size=0)
                else:
                    if f=='cp':
                        #df.plot(kind='barh', ax=ax, legend=False, colormap='Blues', grid=False)
                        ax.barh(N.arange(0,len(df))+0.2, df.values, height=0.6, color='#eaeaea')
                    else:
                        df.plot(kind='barh', ax=ax, legend=False, grid=False)
                    ax.set_xlim(xlims[f])
                    ax.set_xticks(xticks[f])
                    ax.yaxis.set_tick_params(which='both', width=0, size=0)
                    if f in th:
                        ax.plot([th[f],th[f]], ax.get_ylim(), 'r--')
                ax.set_yticks(N.arange(len(df))+0.2)
                # ax.set_yticks(N.arange(len(df))+0.5)
                if j>0:
                    ax.set_yticklabels([])
                else:
                    ax.set_yticklabels(list(df.index))
                if i>0:
                    ax.set_xticklabels([])
                elif modes[f]=='imshow':
                    ax.set_xticklabels(xlabels)
                    for t in ax.xaxis.get_ticklabels():
                        t.set_fontsize(10)
                        t.set_rotation(90)
                if i == len(gos)-1:
                    if modes[f] == 'imshow':
                        cax = fig.add_axes([left,bottom+(panelh[gos[i]]+2*spacer)/figh,width,cbh])
                        cbar = P.colorbar(im, cax=cax, orientation='horizontal')
                        #cax.xaxis.tick_top()
                        cbar.set_ticks(zlims[f])
                        for t in cax.xaxis.get_ticklabels():
                            t.set_fontsize(10)
                        cax.set_title(titles.get(f,f))
                        cax.title.set_fontsize(12)
                    else:
                        ax.set_title(titles.get(f,f))
                        ax.title.set_fontsize(12)
                left += (panelw+spacer)/figw
            bottom += (panelh[gos[i]]+spacer)/figh
        fig.show()
        return fig



class GOEnrichmentAnalysis(object):

    def __init__(self, population_ids, id2gos, go):
        """
        @param population_ids
        @param id2gos (dict): gob,gom,goc -> (i2g (dict): id -> go term list)
        @param go: bioinfo.go.GO object

        """
        self.pop_ids = list(set(population_ids))
        self.id2gos = id2gos
        self.go = go
        self.parsed = {}
        # parse GO
        for w in ['gob','gom','goc']:
            print("parsing %s... " % (w,))
            self.parsed[w] = self._parseGO(w)

    def _parseGO(self, which):
        i2g = self.id2gos[which]
        g2x = self.go.id2idx
        ids = self.go.ids
        s2g = {} # id -> [GO] unique item list no duplications
        for k in self.pop_ids:
            v = [g2x[x] for x in i2g.get(k,[])]
            if len(v)>0: # if it has GO annotation
                s2g[k] = list(set(v))
        g2s = {} # GO -> [ids]
        for s in s2g: # reverse previous map
            for g in s2g[s]:
                g2s.setdefault(g,[]).append(s)
        # make an "envelope" of the selected gos, i.e. add all the parents        
        # gather all parents
        parents = set()
        go = self.go
        for g in g2s:
            #parents = parents.union(set(go.get_ancestors(g)))
            parents = parents.union({g2x[x] for x in go.get_ancestors(ids[g])})
        for p in parents:
            if p not in g2s:
                g2s[p] = []
        for g, v in g2s.items():
            g2s[g] = list(set(v)) # this part is unnecessary
        return s2g, g2s
        
    def calculate(self, selected_ids, which):
        sel_ids = list(set(selected_ids))
        rslt = GOResult(sel_ids, which, self)
        ids = self.go.ids
        s2g,g2s = self.parsed[which] # id->go, go->id
        # first calculate g2s0 the map from GO term -> selected id
        num = 0
        tmp = {}
        g2s0 = {}
        for s in sel_ids:
            if s in s2g:
                num += 1
                for g in s2g[s]:
                    tmp.setdefault(g,[]).append(s)
        for g in tmp:
            g2s0[g] = list(set(tmp[g]))
            rslt.g2i_node[ids[g]] = g2s0[g]
        tot = len(s2g)
        rslt.num = num # number of selected genes with GO associated (total drawn)
        rslt.tot = tot # total number of genes with GO associated  (total in urn)
        p = float(num)/float(tot) # average ratio of selected genes with GO
        # now calculate node counts (has to do for all of the nodes to accomodate descendant counts)
        sum = N.sum
        #hyperg = statsext.hyperg_pdf  # legacy speedup
        for g in g2s: 
            n0 = len(g2s0.get(g,[])) # number of selected with specific GO node
            des = self._descendants(g2s0, g) 
            t0 = len(des) # number of selected with specific GO tree (white drawn)
            if n0==0 and t0==0:
                continue
            idg = ids[g]
            rslt.g2i_tree[idg] = des
            if n0==0:
                rslt.g2i_node[idg] = []
            n1 = len(g2s[g]) # number of total with specific GO node
            t1 = len(self._descendants(g2s, g)) # t1: number of total with specific GO tree (total white)
            rslt.g2cnt[idg] = (n0,n1,t0,t1)
            #n2 = tot-n1 # number of total not specific GO node 
            rslt.g2p_node[idg] = sum(pmf(range(n0,num+1),tot,n1,num))
            #t2 = tot-t1 # number of total not with specific GO tree (total black)
            rslt.g2p_tree[idg] = sum(pmf(range(t0,num+1),tot,t1,num))
            # dhyper(x,m,n,k)
            # x: white drawn  = t0->num
            # m: total white  = t1 
            # n: total black  = t2
            # k: number drawn = num
            rslt.g2over[idg] = float(t0)/float(t1) > p
        rslt.make_table()
        return rslt


    def _descendants(self, g2s, g):
        des = self.go.descendants[g] # use precalculated descendants
        count = {}
        for g in des:
            if g in g2s:
                for s in g2s[g]:
                    count[s] = 1
        return count.keys()    


class GOResult(object):
    
    def __init__(self,sel_ids, which, analysis_obj):
        self.g2cnt = {}# go-> counts (n0,n1,m0,m1) (node,node_total, des, des_total)
        self.g2i_node = {}# go->id for selected
        self.g2i_tree = {}# go->id including descendants
        self.g2p_node = {}# go->p-value for node
        self.g2p_tree = {}# go->p-value for tree
        self.g2over = {} # go->overrepresented? (true) underpresented(false)
        self.sel_ids = sel_ids
        self.num = 0
        self.tot = 0
        self.which = which
        self.analysis = analysis_obj

    def make_table(self):
        columns = ['GO_id','GO_term','tot.node','cnt.node','p.value.node','ids.node',
                                     'tot.tree','cnt.tree','p.value.tree','ids.tree','odds','over']
        rows = []
        ratio = float(self.num)/self.tot
        for g in self.g2cnt:
            n0,n1,t0,t1 = self.g2cnt[g]
            np = self.g2p_node[g]
            tp = self.g2p_tree[g]
            ni = self.g2i_node[g]
            ti = self.g2i_tree[g]
            over = self.g2over[g]
            odds = (float(t0)/t1)/ratio
            term = self.analysis.go.get_name(g)
            rows.append([g,term,n1,n0,np,ni,t1,t0,tp,ti,odds,over])
        self.df = PD.DataFrame(rows, columns=columns).set_index('GO_id')

    def report(self, th_node=1., th_tree=1., num_node=1, num_tree=1, \
               show=1, out=sys.stdout,mode='tree'):
        goids = [(v,k) for k,v in self.g2p_tree.items() if \
                 (v<=th_tree and len(self.g2i_tree[k])>=num_tree)]
        ids_node = [(v,k) for k,v in self.g2p_node.items() if \
                    (v<=th_node and len(self.g2i_node[k])>=num_node)]
        goids.sort()
        ids_node.sort()
        goids = [k for v,k in goids]
        ids_node = [k for v,k in ids_node]
        if mode=='tree':
            ids = goids
        else:
            ids = ids_node
        if show:
            go = self.analysis.go
            out.write( "GO(%s) OverRep Analysis: sel=%d, num=%d, total=%d, th_tree=%g, num_tree=%d, \
            th_node=%g, num_node=%d, chosen from %d\n" % (\
                self.which,len(goids), self.num, self.tot, th_tree, num_tree, \
                th_node, num_node, len(self.g2i_tree)))
            out.write( "GO id\tnode\ttotal\tp-value\ttree\ttotal\tp-value\todds\tname\tover_under\n")
            ratio = float(self.num)/self.tot
            for g in ids:
                n0,n1,t0,t1 = self.g2cnt[g]
                odds = (float(t0)/t1)/ratio
                out.write( "%s\t" % (g,) + \
                           "%d\t%d\t%1.2e\t" % (n0,n1,self.g2p_node[g]) +\
                           "%d\t%d\t%1.2e\t" % (t0,t1,self.g2p_tree[g]) +\
                           "%.3f\t" % (odds,)+\
                           "%s\t" % (go.get_name(g),)+\
                           "%s\n" % (['Under','Over'][self.g2over[g]],)
                          )
        return ids
    
    def write(self,th_tree,num_tree,fname):
        f = open(fname, 'w')
        self.report(th_tree=th_tree,num_tree=num_tree,out=f)

    def make_dot(self, th_tree, num_tree, dotcfg={}):
        seed = self.report(th_tree=th_tree,num_tree=num_tree, show=0)
        if len(seed)==0:
            print("NOTHING TO DRAW")
            return None
        pvalues = self.g2p_tree
        over = self.g2over
        labels = {}
        for g in seed:
            labels[g] = '%s\\n(%d/%d/%1.1e)\\n(%d/%d/%1.1e)' % (
                self.analysis.go.get_name(g),
                self.g2cnt[g][0],self.g2cnt[g][1],self.g2p_node[g],
                self.g2cnt[g][2],self.g2cnt[g][3],self.g2p_tree[g])
        g = Graphviz(self.analysis.go)
        return g.make_dot(seed, pvalues, over, labels,dotcfg=dotcfg)


    def make_dot2(self, root=None, which='gob', depth=2, dotcfg={},cutoff=True):
        go = self.analysis.go
        if root is None:
            root = go.get_root_id(which)
        pvalues = self.g2p_tree
        over = self.g2over
        g = Graphviz(go)
        return g.make_dot2(pvalues, over, root=root, depth=depth, 
                            dotcfg=dotcfg,cutoff=cutoff)
        