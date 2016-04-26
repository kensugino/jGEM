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
	g4sjexprefix = os.path.join(datadir, 'assemblies/gencode.vM4.chr21119')
	ctgt = AN.ComparatorNames(mergedsjexbase, 'merged', outdir)
	cref = AN.ComparatorNames(g4sjexprefix, 'gen4', outdir)
	cp = AN.Comparator(cref, ctgt, outdir, 'gene_name')
	cp.annotate(overwrite=False)
	expath = os.path.join(outdir, 'merged.gen4.ex.txt.gz')
	sjpath = os.path.join(outdir, 'merged.gen4.sj.txt.gz')
	assert os.path.exists(expath)
	assert os.path.exists(sjpath)