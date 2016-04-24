import pytest

import os
import pandas as PD

import jgem
from jgem import gtfgffbed as GGB
from jgem import filenames as FN
from jgem import assembler as AS
from jgem import convert as CV
from jgem import utils as UT

# def pytest_addoption(parser):
#     parser.addoption("--runslow", action="store_true",
#         help="run slow tests")

# slow = pytest.mark.skipif(
#     not pytest.config.getoption("--runslow"),
#     reason="need --runslow option to run"
# )        

@pytest.fixture(scope='session')
def datadir():
	"returns directory containing test data"
	return os.path.join(os.path.dirname(__file__), 'data')

@pytest.fixture(scope='session')
def outdir():
	"returns directory for test outputs"
	return os.path.join(os.path.dirname(__file__), 'out')

@pytest.fixture(scope='session')
def sampleinfo(datadir):
	"returns pandas dataframe containing sample info"
	sipath = os.path.join(datadir, 'sampleinfo.xlsx')
	si = PD.read_excel(sipath)
	si['bwfile'] = [os.path.join(datadir, 'bigwig', x) for x in si['bigwig']]
	si['sjfile'] = [os.path.join(datadir, 'SJ', x) for x in si['sjbed']]
	si['sjexpre'] = [os.path.join(datadir, 'assemblies', x) for x in si['name']]
	return si

@pytest.fixture(scope='session')
def sjtab(sampleinfo, datadir):
	"returns path to test SJ.out.tab file"
	name = sampleinfo.iloc[0]['sjtab']
	return os.path.join(datadir, 'SJ', name)

@pytest.fixture(scope='session')
def sname(sampleinfo, outdir):
	"returns name of sample"
	return sampleinfo.iloc[0]['name']

@pytest.fixture(scope='session')
def sjbed(sampleinfo, outdir, datadir):
	"returns path to test sjbed"
	name0 = sampleinfo.iloc[0]['sjbed']
	sjbed = os.path.join(datadir, 'SJ', name0)
	if not os.path.exists(sjbed):
		name1 = sampleinfo.iloc[0]['sjtab']
		sjtab = os.path.join(datadir, 'SJ', name1)
		aligned = sampleinfo.iloc[0]['aligned']
		sj = GGB.sjtab2sjbed(sjtab, sjbed, aligned)
	return sjbed

@pytest.fixture(scope='session')
def sjexprefix(sampleinfo, datadir):
	"returns prefix to SJ, EX"
	name = sampleinfo.iloc[0]['name']
	return os.path.join(datadir, 'assemblies', name)

@pytest.fixture(scope='session')
def bigwig(sampleinfo, datadir):
	"returns path to test bigwig"
	name = sampleinfo.iloc[0]['bigwig']
	return os.path.join(datadir, 'bigwig', name)

@pytest.fixture(scope='session')
def g4gtfpath(datadir):
	"returns path to gencode4 gtf"
	return os.path.join(datadir, 'assemblies','gencode.vM4.chr1.gtf.gz')

@pytest.fixture(scope='session')
def g4sjexprefix(datadir):
	"returns path to gencode4 gtf"
	pre = os.path.join(datadir, 'assemblies','gencode.vM4.chr1')
	# make sure SJ,EX exists\
	if (not os.path.exists(pre+'.ex.txt.gz')) or (not os.path.exists(pre+'.sj.txt.gz')):
		gtf2sjex = CV.GTF2SJEX(pre+'.gtf.gz')
		sj,ex = gtf2sjex.sjex()
	return pre

@pytest.fixture(scope='session')
def g4gtf(g4gtfpath):
	#GTFCOLS = ['chr','src','typ','st','ed','sc1','strand','sc2','extra']
	#return PD.read_table(g4gtfpath, names=GTFCOLS, compression='gzip')
	return GGB.read_gtf(g4gtfpath)

@pytest.fixture(scope='session')
def g4sjex(g4gtfpath):
	"returns gencode4 chr1 sj, ex dataframes"
	gtf2sjex = CV.GTF2SJEX(g4gtfpath)
	return gtf2sjex.sjex()
	# spath = os.path.join(datadir, 'gencode.vM4.chr1.sj.txt.gz')
	# epath = os.path.join(datadir, 'gencode.vM4.chr1.ex.txt.gz')
	# sj = PD.read_table(spath,compression='gzip')
	# ex = PD.read_table(epath,compression='gzip')
	# return sj,ex

@pytest.fixture(scope='session')
def g4Xkr4gtf(datadir):
	"returns gencod4 Xkr4 gtf exons"
	path = os.path.join(datadir, 'assemblies','gencode.vM4.Xkr4.gtf.gz')
	return GGB.read_gtf(path)

@pytest.fixture(scope='session')
def g4Xkr4sjex(datadir):
	"returns gencode4 Xkr4 sj, ex dataframes"
	path = os.path.join(datadir, 'assemblies','gencode.vM4.Xkr4.gtf.gz')
	g2se = CV.GTF2SJEX(path)
	return g2se.sjex()
	# spath = os.path.join(datadir, 'gencode.vM4.Xkr4.sj.txt.gz')
	# epath = os.path.join(datadir, 'gencode.vM4.Xkr4.ex.txt.gz')
	# sj = PD.read_table(spath,compression='gzip')
	# ex = PD.read_table(epath,compression='gzip')
	# return sj,ex

@pytest.fixture(scope='session')
def mm10chromsizes():
	"returns path to mm10.chrom.sizes"
	d = os.path.dirname(jgem.__file__)
	path = os.path.join(d,'data', 'mm10.chrom.sizes')
	return path

@pytest.fixture(scope='session')
def fnobj(sname, sjbed, bigwig, g4gtfpath, outdir):
	"returns FileNames object"
	return FN.FileNames(sname, bigwig, sjbed, str(outdir), g4gtfpath)

@pytest.fixture(scope='session')
def sj(sjbed):
	"returns sj datafrae"
	sj0 = GGB.read_bed(sjbed)
	return sj0.iloc[:5000]

@pytest.fixture(scope='session')
def asm(fnobj, sj):
	"returns Assembler object"
	a = AS.Assembler(fnobj, merging=False, saveintermediates=False)
	a.sj = sj
	a.params['np'] = 2
	a.params['override'] = True
	return a

@pytest.fixture(scope='session')
def testbam(datadir):
	return os.path.join(datadir, 'bedtools/test.bam')

@pytest.fixture(scope='session')
def testbed(datadir):
	return os.path.join(datadir, 'bedtools/test.bed.gz')

@pytest.fixture(scope='session')
def testbed12(datadir):
	return os.path.join(datadir, 'bedtools/test.bed12.gz')

@pytest.fixture(scope='session')
def testwig(datadir):
	return os.path.join(datadir, 'bedtools/test.wig')

@pytest.fixture(scope='session')
def testchromsizes():
	return UT.chromsizes('dm3')

@pytest.fixture(scope='session')
def testsampleinfo(datadir):
	si = UT.read_pandas(os.path.join(datadir, 'bedtools/test-si.txt'))
	si['bwfile'] = datadir + '/' + si['bwfile']
	return si


