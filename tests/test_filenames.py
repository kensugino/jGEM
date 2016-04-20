import os
import pytest
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

from jgem import filenames as FN
from jgem import gtfgffbed as GGB

def test_FileNamesBase(tmpdir, g4Xkr4sjex):
	sj,ex = g4Xkr4sjex # example DataFrame
	# instantiation
	p = str(tmpdir)
	fnb = FN.FileNamesBase(prefix=p)
	t1 = p+'.test1.txt.gz'
	assert fnb.fname('test1.txt.gz',category='temp') == t1
	t2 = p+'.test2.txt.gz'
	assert fnb.txtname('test2') == t2
	t3 = p+'.test3.bed.gz'
	assert fnb.bedname('test3') == t3
	sjpath = p+'.sj.txt.gz'
	assert fnb.write_txt(sj, 'sj') == sjpath
	assert os.path.exists(sjpath)
	expath = p+'.ex.bed.gz'
	#LOG.debug(ex.columns)
	assert fnb.write_bed(ex, 'ex', ncols=6) == expath
	assert os.path.exists(expath)
	sj1 = fnb.read_txt('sj')
	assert (sj1 == sj).all().all()
	ex1 = fnb.read_bed('ex')
	assert (ex1 == ex[GGB.BEDCOLS[:6]]).all().all()
	t4 = p+'.test4.txt.gz'
	assert fnb.txtname('test4', category='c1') == t4
	assert fnb._fnames['temp'] == set([t1,t2,t3,sjpath,expath])
	assert fnb._fnames['c1'] == set([t4])
	fnb.delete(['temp'])
	assert os.path.exists(sjpath) == False
	assert os.path.exists(expath) == False
	fnb.delete(['c1'])

def test_FileNames():
	fn = FN.FileNames('sname','/path/to/bwfile','/path/to/sjfile','/path/to/outdir')
	assert fn._prefix == '/path/to/outdir/sname'
	assert fn._fnames == {}
