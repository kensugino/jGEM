import os
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

import pytest
import pandas as PD
import numpy as N

from jgem import utils as UT
from jgem import collector as CL


# def test_collect_jcnt(sampleinfo, mergedsjexbase, datadir, outdir):
# 	c = CL.Collector('acode', mergedsjexbase, 'scode', sampleinfo, outdir)
# 	c.collect_jcnt()

# def test_collect_ecov(sampleinfo, mergedsjexbase, datadir, outdir):
# 	c = CL.Collector('acode', mergedsjexbase, 'scode', sampleinfo, outdir)
# 	c.collect_ecov()	

# def test_collect_gcov(sampleinfo, mergedsjexbase, datadir, outdir):
# 	c = CL.Collector('acode', mergedsjexbase, 'scode', sampleinfo, outdir)
# 	c.collect_gcov()		

# def test_collect_gcov1000(sampleinfo, mergedsjexbase, datadir, outdir):
# 	c = CL.Collector('acode', mergedsjexbase, 'scode', sampleinfo, outdir)
# 	c.collect_gcov1000(unique=False)	



