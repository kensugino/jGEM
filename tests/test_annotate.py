import os
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

import pytest
import pandas as PD
import numpy as N

from jgem import utils as UT
from jgem import annotate as AN


def test_comparator(datadir, outdir, mergedsjexbase, g4sjexprefix):
	#g4sjexprefix = os.path.join(datadir, 'assemblies/gencode.vM4.chr21119')
	ctgt = AN.ComparatorNames(mergedsjexbase, 'merged', outdir)
	cref = AN.ComparatorNames(g4sjexprefix, 'gen4', outdir)
	cp = AN.Comparator(cref, ctgt, 'gene_name')
	cp.annotate(overwrite=False)
	expath = os.path.join(outdir, 'merged.gen4.ex.txt.gz')
	sjpath = os.path.join(outdir, 'merged.gen4.sj.txt.gz')
	assert os.path.exists(expath)
	assert os.path.exists(sjpath)
	assert N.sum([',' in x for x in cp.g2gs.values()]) == 2 # genes 
	ex0 = cp.ex_tgt
	ex0 = ex0[ex0['gen4_gidx0'].notnull()]
	ex1 = UT.read_pandas(expath)
	ex1 = ex1[ex1['gen4_gidx0'].notnull()]
	# ex0 == ex1 ?
	assert N.sum(ex0['gen4_gidxs'].str.contains(',')) == 27 # exons
	assert N.sum(ex0['gen4_gidx0'] != ex1['gen4_gidx0']) == 0
	assert N.sum(ex0['gen4_sym0'] != ex1['gen4_sym0']) == 0
	assert N.sum(ex0['gen4_syms'] != ex1['gen4_syms']) == 0
	# check gidx0
	# gidx0 column should be all single (i.e. does not contain ',')
	idx = ex0['gen4_gidx0'].astype(str).str.contains(',')
	assert N.sum(idx)==0
	idx = ex1['gen4_gidx0'].astype(str).str.contains(',')
	assert N.sum(idx)==0
	# check sym0
	idx = ex0['gen4_sym0'].astype(str).str.contains(',')
	assert N.sum(idx)==0 
	idx = ex1['gen4_sym0'].astype(str).str.contains(',')
	assert N.sum(idx)==0 

