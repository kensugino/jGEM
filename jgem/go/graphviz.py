"""
    Copyright (c) by Ken Sugino, 2008


"""

import math
import numpy as N
#from sets import set
#from rpy import r
#from bioinfo.lib.configattributes import *

class Graphviz(object):
    
    def __init__(self,go):
        self.go = go
        self.id2idx = go.id2idx
        self.ids = go.ids
        self.parents = go.parents
        self.children = go.children
        self.roots = go.roots
            
    def graph1(self, seed):
        """
        Create tree structure from the set:

        @type seed: list or tuple or set
        @param seed: a seed of the tree (list of GO id's)
        @rtype : edges:list of (GO id, GO id) 
        """
        i2x = self.id2idx
        #x2i = self.ids
        x2p = self.parents
        d2g = {}
        edges = set()
        def _calc(nodeidx):
            for x in x2p[nodeidx]:
                edges.add('%d;%d' % (x, nodeidx))
                _calc(x)
        for goid in seed:
            _calc(i2x[goid])
        edges = [tuple(x.split(';')) for x in edges]
        for x,y in edges:
            d2g[x] = [int(x)]
            d2g[y] = [int(y)]
        return d2g, edges
        
    def calc_levels(self):
        #if hasattr(self, 'levels'):
        #    return self.levels
        x2c = self.children
        i2x = self.id2idx
        #x2l = [{} for x in xrange(len(x2i))] # (parent->level map)
        x2l = [-1]*len(x2c) # minimum level
        def _calc(node, level):
            if x2l[node]>0:
                minlevel = min(x2l[node], level)
            else:
                minlevel = level
            x2l[node] = minlevel
            for x in x2c[node]:
                _calc(x, minlevel+1)
        for root in self.roots.values():
            _calc(i2x[root], 0)
        #self.levels = x2l
        return x2l
    
    def graph2(self, root='GO:0008150', depth=2):
        """
            RECOVER UNIQUENESS!
            all "root,node" -> rename to "node"
            all "root,(any node chain of length<depth),node" -> rename to "node"
             
        """
        x2c, i2x = self.children, self.id2idx        
        d2g = {} # map from deconvoluted node -> [go node]
        edges = set() # dcnode;dcnode
        #x2l = self.calc_levels()
        def _calc(gonode, dcnode, level):
            d2g.setdefault(dcnode, []).append(gonode)
            for x in x2c[gonode]:
                if level>=depth:
                    ndcnode = dcnode+','
                else:
                    ndcnode = str(x)
                edges.add('%s;%s' % (dcnode, ndcnode))
                _calc(x, ndcnode,level+1)

        _calc(i2x[root],str(i2x[root]),0)
        edges = [tuple(x.split(';')) for x in edges] 
        return d2g, edges
    
    
    def graph3(self, root='GO:0008150', depth=2):
        x2c, i2x = self.children, self.id2idx        
        d2g = {} # map from deconvoluted node -> [go node]
        edges = set() # dcnode;dcnode
        x2l = self.calc_levels()
        def _calc(gonode, dcnode):
            d2g.setdefault(dcnode, []).append(gonode)
            for x in x2c[gonode]:
                if x2l[x]<=depth:
                    #ndcnode = dcnode+','
                    #do not add edges
                    #else:
                    ndcnode = str(x)
                    edges.add('%s;%s' % (dcnode, ndcnode))
                    _calc(x, ndcnode)

        _calc(i2x[root],str(i2x[root]))
        edges = [tuple(x.split(';')) for x in edges] 
        return d2g, edges
        
        
        
    def make_dot2(self, pvalues={}, over={}, colors=[], root='GO:0008150', depth=2, 
                              plog=True, dotcfg={}, usename=True, cutoff=True):
        """" 
            pvalues: (dict) go id -> pvalue : pvalue will be mapped to color
            over: (dict) go id -> True (overrepresented) or False (underrepresented)
            colors: (over_colors, under_colors): colors for over or under represented pvalues
            root: root node go id
            depth: where to start representing all items in one level as one node
            plog: whether to take log10 to pvalues
            dotcfg: config for dot
            
        """
        pcolors, ncolors, npcolors, nncolors = _dot_colors(colors)
        dotheader, cfg = _dot_config(dotcfg,{'nodesep':'0.2','ranksep':'1.5',
                                             'nodeheight':'0.5'})

        if cutoff:
            d2g, edges = self.graph3(root, depth)
        else:
            d2g, edges = self.graph2(root, depth)
        
        i2x = self.id2idx
        x2i = self.ids
        pmax, pmin = _dot_pvalues(pvalues,plog,d2g,i2x)
        dotnodes = d2g.keys()
        for i, dcnode in enumerate(dotnodes):
            dcname = dcnode.split(',')[-1]
            if dcname:
                dcname = x2i[int(dcname)]
            nodecolor = ''
            shape = 'ellipse'
            goids = [x2i[x] for x in d2g[dcnode]]
            pvals = [pvalues[x] for x in goids if x in pvalues]
            if len(pvals)>0:
                gmin = min(pvals)
                goid = [x for x in goids if ((x in pvalues) and (pvalues[x]==gmin))][0]
                gmin = min(gmin,1.0)
                if plog:
                    gmin = math.log10(gmin)
                overs = [over[x] for x in goids if x in over]
                mixed = (True in overs) and (False in overs)
                if mixed:
                    shape = 'circle'
                if over[goid]: 
                    nodecolor = pcolors[int(math.floor(npcolors*(gmin-pmax)/(pmin-pmax)))]
                else:
                    nodecolor = ncolors[int(math.floor(nncolors*(gmin-pmax)/(pmin-pmax)))]
            if dcname and usename:
                dcname = self.go.get_name(dcname)
            if nodecolor:
                opt = '[label="%s",shape="%s",style="filled", fillcolor="%s"]' % (dcname,shape,nodecolor)
            else:
                opt = '[label="%s",shape="%s"]' % (dcname,shape)
            dotnodes[i] = '"%s" %s;' % (dcnode, opt)
        dotnodes = '\n'.join(dotnodes)+'\n'
        dotedges = '\n'.join(['"%s" -- "%s";' % x for x in edges]) 
        return dotheader+dotnodes+dotedges+"\n}\n"
    
    
    def make_dot(self, seed, pvalues={}, over={}, labels={}, colors=[], 
                 plog=True, dotcfg={}, usename=True, breaklen=10):
        
        pcolors, ncolors, npcolors, nncolors = _dot_colors(colors)
        dotheader, cfg = _dot_config(dotcfg,{'nodesep':'0.1','ranksep':'1'})

        d2g, edges = self.graph1(seed)
        i2x = self.id2idx
        x2i = self.ids
        pmax, pmin = _dot_pvalues(pvalues,plog,d2g,i2x)
        dotnodes = d2g.keys()
        for i, dcnode in enumerate(dotnodes):
            goid = x2i[d2g[dcnode][0]]
            nodecolor = ''
            shape = cfg['shape']
            if goid in pvalues:
                pval = min(pvalues[goid],1.0)
                if plog:
                    pval = math.log10(pval)
                if over.get(goid, True):
                    nodecolor = pcolors[int(math.floor(npcolors*(pval-pmax)/(pmin-pmax)))]
                else:
                    nodecolor = ncolors[int(math.floor(nncolors*(pval-pmax)/(pmin-pmax)))]
            name = goid
            if usename:
                name = self.go.get_name(goid)
            if goid in labels:
                name = labels[goid]
            if breaklen>0:
                name = _dot_divide_line(name,breaklen)
            if nodecolor:
                opt = '[label="%s",shape="%s",style="filled", fillcolor="%s"]' % (name,shape,nodecolor)
            else:
                opt = '[label="%s",shape="%s"]' % (name,shape)
            dotnodes[i] = '"%s" %s;' % (dcnode, opt)
        dotnodes = '\n'.join(dotnodes)+"\n"
        dotedges = '\n'.join(['"%s" -- "%s";' % x for x in edges]) 
        return dotheader+dotnodes+dotedges+"\n}\n"

            
def _dot_colors(colors):
    # Calculate colors
    if not colors:
        #seq = r.seq(255,0,-1)
        seq = N.arange(255,-1,-1)
        #pcolors = r.rgb(255,seq,seq,maxColorValue=255)  # red 
        pcolors = ['#%02x%02x%02x' % x for x in zip([255]*len(seq),seq,seq)]
        #ncolors = r.rgb(seq,255,seq,maxColorValue=255)  # green
        ncolors = ['#%02x%02x%02x' % x for x in zip(seq,[255]*len(seq),seq)]
    else:
        pcolors, ncolors = colors
    npcolors,nncolors = len(pcolors)-1, len(ncolors)-1
    return pcolors,ncolors,npcolors,nncolors

def _dot_pvalues(pvalues,plog,d2g,i2x):
    pmax, pmin = 1.0, 0.0
    if len(pvalues)>0:
        goidx = [x for y in d2g.values() for x in y]
        pvals = [pvalues[x] for x in pvalues if i2x[x] in goidx]
        pmax, pmin = max(pvals), min(pvals)
        print("pvalue max:%g, min:%g" % (pmax, pmin))
        pmax = 1.0
        if plog:
            pmax,pmin = map(math.log10, (pmax,pmin))
    return pmax, pmin


"""
       style = 'filled' | 'invisible' | 'diagonals' | 'rounded'
       shape = 'box' | 'ellipse' | 'circle' | 'point' | 'triangle'

       style     = 'dashed' | 'dotted' | 'solid' | 'invis' | 'bold'
       arrowhead = 'box' | 'crow' | 'diamond' | 'dot' | 'inv' | 'none' | 'tee' | 'vee'
       weight    = number (the larger the number the closer the nodes will be)
"""
def _dot_config(dotcfg,cfg2={}):
    # dot config
    cfg = dict(
               nodesep="0.2",
               ranksep="3",
               shape="box",
               nodeheight="0.5",
               nodewidth="1",
               nodestyle="filled",
               edgestyle="solid",
               size="4,6",
               fontsize="10",
               rankdir="LR",
               splines="true",
               graph="graph",
               arrowhead="normal",
               )
    cfg.update(cfg2)
    cfg.update(dotcfg)
    dotheader = """%(graph)s go{
    size = "%(size)s";
    nodesep = "%(nodesep)s";
    ranksep = "%(ranksep)s";
    rankdir = "%(rankdir)s";
    splines = "%(splines)s";
    node [fontsize="%(fontsize)s",height="%(nodeheight)s",width="%(nodewidth)s",fixedsize="false"];
    edge [arrowhead="%(arrowhead)s"];\n""" % cfg
    return dotheader, cfg

def _dot_divide_line(txt,breaklen=10):
    words = txt.split(' ')
    tmplist = []
    tmpline = words[0]
    for w in words[1:]:
        if len(tmpline+w)>breaklen:
            tmplist.append(tmpline)
            tmpline = w
        else:
            tmpline += ' ' + w
    tmplist.append(tmpline)
    return '\\n'.join(tmplist)


