import os
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

import pytest
import pandas as PD
import numpy as N
import subprocess

from jgem import utils as UT
from jgem import calccov as CC


# def test_chopintervals():
# 	pass

# def test_calc_cov_ovl_mp():
# 	pass

# def test_calc_ecov_chrom():
# 	pass

# def test_calc_ecov_mp():
# 	pass

# def test_calc_cov_mp():
# 	pass

# def test_calc_max_chrom():
# 	pass

# def test_calc_cov_chrom():
# 	pass

# def test_calc_sjcnt():
# 	pass

def test_calc_ecov(g4sjexprefix, bigwig, outdir):
	expath = g4sjexprefix+'.ex.txt.gz'
	cipath = g4sjexprefix+'.ci.txt.gz'
	dstpre = os.path.join(outdir, 'g4')
	ecov = CC.calc_ecov(expath,cipath,bigwig,dstpre,True,1)
	assert 'ecov' in ecov.columns
	LOG.debug(ecov.head())
	# 117     chr1    4857813 4857976 3.0008723189128688
	ecov.set_index('eid',inplace=True)
	assert ecov.ix[117]['chr'] == 'chr1'
	assert ecov.ix[117]['st'] == 4857813 
	assert ecov.ix[117]['ed'] == 4857976 
	assert ecov.ix[117]['ecov'] == 3.0008723189128688

# def test_calc_gcov():
# 	pass

# def test_calc_gcov2000():
# 	pass

# def test_calc_gene_cov():
# 	pass

# def test_calc_glen():
# 	pass


