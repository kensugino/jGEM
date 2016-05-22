"""

.. module:: phylo
    :synopsis: a collection of phylo reloated stuffs (PhyloCSF, Phylo60)

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
import subprocess
import gzip
import os
import errno
import math
import sys
import csv
import time
import collections
import functools
from functools import partial, reduce
from operator import iadd
from itertools import repeat, groupby
try:
    from itertools import izip_longest as zip_longest
except:
    from itertools import zip_longest
import json
import shutil
import uuid

import pandas as PD
import numpy as N
import scipy.optimize as SO

import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

from jgem import utils as UT
from jgem import bigwig as BW

class PhyloCSF(object):
    """ Wrapper to access PhyloCSF bigwig. 

    """
    FRAMES = [0,1,2]
    STRANDS = ['+','-']
    FNAMETMPL = 'PhyloCSFwind01{strand}{frame}.bw'
    SMFNAMETMPL = 'PhyloCSF{strand}{frame}.bw'

    def __init__(self, phylocsfdir):
        self.pcdir = phylocsfdir
        self.fobjs = {'+':{},'-':{}}
        self.bws = {'+':{},'-':{}}

    def __enter__(self):
        # self.fobjs = {'+':{},'-':{}}
        # self.bws = {'+':{},'-':{}}
        for s in self.STRANDS:
            for f in self.FRAMES:
                bwname = self.FNAMETMPL.format(strand=s,frame=f)
                bwpath = os.path.join(self.pcdir, bwname)
                self.fobjs[s][f] = fp = open(bwpath, 'rb')
                self.bws[s][f] = BW.BigWigFile(fp)
        return self
        
    def __exit__(self, exc_type, exc_value, traceback):
        for s in self.STRANDS:
            for f in self.FRAMES:
                self.fobjs[s][f].close()


    def get(self, s, f, chrom, st, ed):
        bw = self.bws[s][f]
        a = bw.get_as_array(chrom,st,ed)
        if a is None:
            a = N.array([]) # null array
        else:
            a[N.isnan(a)]=0.
        return a

    def score(self, s, f, chrom, st, ed):
        a = self.get(s,f,chrom,st,ed)
        # a[a<0]=0
        return N.sum(a)

    def maxscore(self,s,chrom,st,ed):
        return N.max([self.score(s,f,chrom,st,ed) for f in self.FRAMES])

    def calc_scores(self, gbed):
        # gbed : beddf subset containing uexons in one gene
        # adds score+, score- columns to gbed
        for s in self.STRANDS:
            colname = 'score{0}'.format(s)
            gbed[colname] = [self.maxscore(s,chrom,st,ed) for chrom,st,ed in gbed[['chr','st','ed']].values]
        return gbed

    def calculate(self, unionexbed, addcols=['_id','_gidx'], np=10):
        """ Calculate PhyloCSF score.

        Args:
            unionexbed: Pandas DataFrame bed containing unioned exons
            addcols: additional cols to copy from unionexbed (other than chr,st,ed)
            np: number of CPU to use

        Returns:
            A dictionary: gene_id => score

        """
        # process chrom wise
        args = []
        for chrom in unionexbed['chr'].unique():
            uechr = unionexbed[unionexbed['chr']==chrom][['chr','st','ed']+addcols].copy()
            mycopy = PhyloCSF(self.pcdir) # looks like MP doesn't get invoked unless you use a copy...
            args.append((uechr, mycopy))

        rslts = UT.process_mp(calc_worker, args, np=np, doreduce=False)

        df = PD.concat(rslts, ignore_index=True)
        return df

def calc_worker(uex, pcobj):
    with pcobj:
        df = pcobj.calc_scores(uex)
    return df




class Phylo60(object):

    def __init__(self, path):
        self.path = path

    def __enter__(self):
        self.fobj = fp = open(self.path, 'rb')
        self.bw = BW.BigWigFile(fp)
        return self
        
    def __exit__(self, exc_type, exc_value, traceback):
        self.fobj.close()

    def get(self, chrom, st, ed):
        bw = self.bw
        a = bw.get_as_array(chrom,st,ed)
        if a is None:
            a = N.array([]) # null array
        else:
            a[N.isnan(a)]=0.
        return a

    def score(self, chrom, st, ed):
        a = self.get(chrom,st,ed)
        return N.sum(a)

    def calc_scores(self, gbed):
        # gbed : beddf subset containing uexons in one gene
        # adds phylo60score column to gbed
        colname = 'phylo60score'
        gbed[colname] = [self.score(chrom,st,ed) for chrom,st,ed in gbed[['chr','st','ed']].values]
        return gbed

    def calculate(self, unionexbed, addcols=['_id','_gidx'], np=10):
        """ Calculate PhyloCSF score.

        Args:
            unionexbed: Pandas DataFrame bed containing unioned exons
            addcols: additional cols to copy from unionexbed (other than chr,st,ed)
            np: number of CPU to use

        """
        # process chrom wise
        args = []
        for chrom in unionexbed['chr'].unique():
            uechr = unionexbed[unionexbed['chr']==chrom][['chr','st','ed']+addcols].copy()
            mycopy = Phylo60(self.path)
            args.append((uechr, mycopy))

        rslts = UT.process_mp(calc_worker, args, np=np, doreduce=False)

        df = PD.concat(rslts, ignore_index=True)
        return df

