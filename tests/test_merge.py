import os
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

import pytest
import pandas as PD
import numpy as N

from jgem import utils as UT
from jgem import merge as MG
from jgem import gtfgffbed as GGB

def test_MergeInputNames(testsampleinfo, outdir, datadir):
	with pytest.raises(ValueError):
		fn = MG.MergeInputNames(testsampleinfo, 'test-mergeinputnames', outdir)
	fn = MG.MergeInputNames(testsampleinfo, 'test-mergeinputnames', outdir, checkfiles=False)
	p = os.path.join(outdir, 'test-mergeinputnames.')
	expaths = fn.expaths()
	assert ('sample1', 'assemblies/s1.ex.txt.gz') == expaths[0]
	assert ('sample2', 'assemblies/s2.sj.txt.gz') == fn.sjopaths()[1]
	LOG.debug(fn.sjpaths()[2])
	assert ('sample3', os.path.join(datadir,'SJ/test3.sj.bed.gz')) == tuple(fn.sjpaths()[2])
	LOG.debug(fn.bwpaths()[3])
	assert ('sample4', os.path.join(datadir,'bedtools/test.bw')) == tuple(fn.bwpaths()[3])
	assert fn.sj0_bed() == p + 'sj0.bed.gz'
	assert fn.sj_bed('p') == p + 'sj.p.bed.gz'
	assert fn.allsj_txt() == p + 'allsj.txt.gz'
	assert fn.allsj_stats() == p + 'allsj.stats.txt.gz'
	assert fn.ex_bw('mep') == p + 'ex.mep.bw'
	assert fn.ex_bed('men') == p + 'ex.men.bed.gz'
	assert fn.agg_bw() == p + 'allsample.bw'
	assert fn.snames() == ['sample1','sample2','sample3','sample4']


def test_agg_bw2(sampleinfo, outdir):
	fn = MG.MergeInputNames(sampleinfo, 'Fev_merge_test', outdir)
	mi = MG.MergeInputs(fn, genome='mm10', np=3)
	mi.aggregate_bigwigs()
	assert os.path.exists(fn.agg_bw())

def test_make_sj_bed(sampleinfo, outdir):
	fni = MG.MergeInputNames(sampleinfo, 'Fev_merge_test', outdir)
	mi = MG.MergeInputs(fni, genome='mm10', np=3)
	mi.make_sj_bed()
	assert os.path.exists(fni.sj0_bed())
	assert os.path.exists(fni.allsj_txt())
	assert os.path.exists(fni.allsj_stats())
	assert os.path.exists(fni.sj_bed('p'))
	assert os.path.exists(fni.sj_bed('n'))
	# sj include unstranded junction
	sjp = GGB.read_sj(fni.sj_bed('p'))
	sjp['locus'] = UT.calc_locus_strand(sjp)
	assert 'chr1:118545460-118547864:.' in sjp['locus'].values


def test_make_ex_bigwigs(sampleinfo, outdir):
	fn = MG.MergeInputNames(sampleinfo, 'Fev_merge_test', outdir)
	mi = MG.MergeInputs(fn, genome='mm10', np=3)
	mi.make_ex_bigwigs()
	assert os.path.exists(fn.ex_bw('mep'))
	assert os.path.exists(fn.ex_bw('men'))
	assert os.path.exists(fn.ex_bw('se'))
	#assert os.path.exists(fn.ex_bw('sep'))
	#assert os.path.exists(fn.ex_bw('sen'))

# def test_aggregate_bigwigs(testsampleinfo, outdir):
# 	fn = MG.MergeInputNames(testsampleinfo, 'test-mergeinputnames', outdir)
# 	mi = MG.MergeInputs(fn, genome='dm3', np=1)
# 	mi.aggregate_bigwigs()
# 	assert os.path.exists(fn.agg_bw())

def test_MergeAssemblyNames(outdir):
	od = os.path.join(outdir,'_merge')
	fn = MG.MergeAssemblyNames('test.man', od)
	assert fn._prefix == od+'/test.man'
	assert fn.ex_out('bed') == od+'/test.man.ex.bed.gz'
	assert fn.sj_out('txt') == od+'/test.man.sj.txt.gz'
	assert fn.ci_out() == od+'/test.man.ci.txt.gz'
	assert fn.genes_out('bed') == od+'/test.man.genes.bed.gz'


def test_MergeAssemble(sampleinfo, outdir, datadir):
	fni = MG.MergeInputNames(sampleinfo, 'Fev_merge_test', outdir=outdir)
	mi = MG.MergeInputs(fni, genome='mm10', np=1, th_detected=0, th_maxcnt1=0, th_ratio=0.0005)#, jie_sjth=1)
	if not os.path.exists(fni.agg_bw()):
		mi.aggregate_bigwigs()
	if not os.path.exists(fni.ex_bw('men')):
		mi.make_ex_bigwigs()
	if not os.path.exists(fni.sj_bed('p')):
		mi.make_sj_bed()
	assert os.path.exists(fni.ex_bw('men'))
	assert os.path.exists(fni.ex_bw('mep'))
	assert os.path.exists(fni.ex_bw('se'))
	assert os.path.exists(fni.sj_bed('p'))
	assert os.path.exists(fni.sj_bed('n'))
	assert os.path.exists(fni.agg_bw())
	fna = MG.MergeAssemblyNames('Fev_merge_asm_test', os.path.join(outdir,'_merge'))
	ma = MG.MergeAssemble(fni, fna, np=1, se_maxth=2)
	ma.assemble()
	assert os.path.exists(fna.ex_out('txt'))
	assert os.path.exists(fna.sj_out('txt'))
	assert os.path.exists(fna.genes_out('txt'))
	assert os.path.exists(fna.ci_out())

