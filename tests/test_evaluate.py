import os
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

import pytest
import pandas as PD
import numpy as N

from jgem import utils as UT
from jgem import evaluate as EV


def test_evalnames(g4gtfpath, outdir):
	sjexbase = g4gtfpath[:-len('.gtf.gz')]
	en = EV.EvalNames(sjexbase, 'gen4', outdir)
	p = os.path.join(str(outdir), 'gen4')
	assert en._prefix == p
	assert en.modelpath('sj') == sjexbase+'.sj.txt.gz'
	assert len(en.model('sj')) == 17631
	assert en.fname2('test.txt.gz','pndr') == p+'.pndr.test.txt.gz'
	assert en.fname('test2.bed.gz') == p+'.test2.bed.gz'

def test_evalmatch(sjexprefix, g4sjexprefix, outdir, bigwig, sjbed):
	en1 = EV.EvalNames(g4sjexprefix, 'gen4', outdir)
	en2 = EV.EvalNames(sjexprefix, 'jgem', outdir)
	em = EV.EvalMatch(en1, en2, bigwig, sjbed, 'FevDR1623',100)
	em.calculate(np=1)
	s1 = en1.model('sj')
	s2 = en2.model('sj')
	e1 = en1.model('ex')
	e2 = en2.model('ex')
	assert hasattr(em, 'e1')
	assert hasattr(em, 's1')
	# prep_sjex
	assert 'len' in s1.columns
	assert 'len' in e2.columns
	assert em.colname('ecov') == 'ecov_FevDR1623'
	assert em.colname2('jhit','jgem') == 'jhit_FevDR1623_jgem'
	assert 'ecov_FevDR1623' in e1.columns
	assert 'ecov_FevDR1623' in e2.columns
	assert 'gcov_FevDR1623' in e1.columns
	assert 'gcov_FevDR1623' in e2.columns
	assert 'jcnt_FevDR1623' in s1.columns
	assert 'jcnt_FevDR1623' in s2.columns
	# find_match
	for k,v in em.e.items():
		LOG.info('{0}:{1}'.format(k,len(v)))
	for k in ['i','5','3','s']:
		v = em.e[k]
		assert len(v['cat'].unique())==1
		assert v['cat'].iloc[0] == k
	assert 'jhit_FevDR1623_jgem' in s1.columns
	# calc_stats
	assert len(em.closest) == 4
	assert len(em.stats) == 5
	assert len(em.ov) == 30449
	assert len(em.e['i']) == 5635
	# sensitivity plot
	axr = em.plot_sensitivity()
	# ratio fig
	fig = em.plot_ratio()
	# number fig
	# completeness calc, fig



