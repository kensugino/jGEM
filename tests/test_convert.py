import os
import pytest
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

from jgem import utils as UT
from jgem import convert as CV
from jgem import gtfgffbed as GGB


def test_bed2exonsj(testbed12):
	b12 = GGB.read_bed(testbed12)
	sj,ex = CV.bed2exonsj(b12)
	print(sj.iloc[:10])
	print(ex.iloc[:10])
	# assert 0



## [TODO] better testing on the followings?
def test_gtf2exonsj(g4Xkr4gtf, g4Xkr4sjex):
	sj0,ex0 = g4Xkr4sjex
	sj1,ex1 = CV.gtf2exonsj(g4Xkr4gtf)
	assert set(sj0['locus'].values)==set(sj1['locus'].values)
	assert set(ex0['locus'].values)==set(ex1['locus'].values)


# def test_bed2exonsj():
# 	pass

def test_gtf2sjex(g4gtfpath, outdir):
	g2sj = CV.GTF2SJEX(g4gtfpath, outdir)
	assert g2sj.exists()
	base = os.path.basename(g4gtfpath[:-7])
	pre = os.path.join(outdir, base)
	assert g2sj.fname('test') == pre+'.test'
	sjpath,expath = g2sj.sjexpaths()
	assert sjpath == pre+'.sj.txt.gz'
	assert expath == pre+'.ex.txt.gz'
	if os.path.exists(sjpath):
		os.unlink(sjpath)
	if os.path.exists(expath):
		os.unlink(expath)
	sj,ex = g2sj.make_sjex()
	assert os.path.exists(sjpath)
	assert os.path.exists(expath)
	assert '_id' in ex.columns
	assert '_gidx' in ex.columns # find_gene
	assert 'd_id' in ex.columns # set_ad_info
	assert 'cat' in ex.columns # category
	LOG.info('len(sj)={0}, len(ex)={1}'.format(len(sj),len(ex)))
	assert len(sj) == 17631
	assert len(ex) == 27151


def test_gtf2sjex2(g4gtfpath, outdir):
	g2sj = CV.GTF2SJEX(g4gtfpath, outdir)
	cipath = g2sj.cipath()
	if os.path.exists(cipath):
		os.unlink(cipath)
	ci = g2sj.ci()
	assert os.path.exists(cipath)
	assert len(ci) == 28776


