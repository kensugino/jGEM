import os
import pytest
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)
import glob

import pandas as PD
import numpy as N

from jgem import assembler as AS
from jgem import gtfgffbed as GGB
from jgem import utils as UT

# TODO: rather than checking number of elements which can change with change in algorithm
#       check actual existence of element (e.g. Snap25, Gapdh exons that should always be there)


GBED = os.path.abspath(os.path.join(os.path.dirname(__file__), 'out/Fev_DR_m70_1623.genes.bed.gz'))

def test_cleanup(fnobj):
	fpat = fnobj.fname('*')
	flist = glob.glob(fpat)
	LOG.info('#files to delete {0}'.format(len(flist)))
	fnobj.delete_prefixed()
	flist2 = glob.glob(fpat)
	assert len(flist2)==0

def test_selectsj(asm):
	LOG.debug(asm.sj.head())
	f = AS.SELECTSJ(asm)
	assert f.fnobj == asm.fnobj
	assert f.params == asm.params
	chroms = ['chr1','chr2','chrX','chrM','chrZ']
	df = PD.DataFrame({'chr':chroms})
	assert f.chroms(df) == ['chr1','chr2','chrX']
	assert len(asm.sj) == 5000
	f()
	#LOG.info('{0}'.format(len(asm.sj)))
	assert N.sum(~asm.sj['strand'].isin(['+','-','.']))==0
	assert len(asm.sj) == 4991
	

def test_checksjsupport(asm):
	#asm.params['binth'] = 0.1
	f = AS.CHECKSJSUPPORT(asm)
	assert len(asm.sj) == 4991
	f()
	LOG.info('{0}'.format(len(asm.sj)))
	assert N.sum(~asm.sj['strand'].isin(['+','-','.']))==0
	assert len(asm.sj) == 4991
	assert os.path.exists(asm.fnobj.bedname('checksjsupport.sj'))
	

def test_removejie(asm):
	f = AS.REMOVEJIE(asm)
	f()
	LOG.info('{0}'.format(len(asm.sj)))
	assert N.sum(~asm.sj['strand'].isin(['+','-','.']))==0
	assert len(asm.sj) == 4986
	

def test_sj2ex(asm):
	f = AS.SJ2EX(asm)
	f()
	LOG.info('{0},{1}'.format(len(asm.me),len(asm.sj)))
	assert N.sum(~asm.me['strand'].isin(['+','-','.']))==0
	assert len(asm.me) == 5280
	assert len(asm.sj) == 4986
	

def test_mergeexons(asm):
	f = AS.MERGEEXONS(asm)
	f()
	LOG.info('{0},{1}'.format(len(asm.me),len(asm.sj)))
	assert N.sum(~asm.me['strand'].isin(['+','-','.']))==0
	assert len(asm.me) == 5341
	assert len(asm.sj) == 4986
	

# def test_addjie(asm):
# 	pass

def test_findedges2(asm):
	f = AS.FINDEDGES2(asm)
	assert len(asm.sj) == 4986
	assert len(asm.me) == 5341
	f()
	LOG.info('{0},{1}'.format(len(asm.me),len(asm.sj)))
	assert N.sum(~asm.me['strand'].isin(['+','-','.']))==0
	assert len(asm.me) == 6039
	assert len(asm.sj) == 4959
	

def test_findedges(asm):
	f = AS.FINDEDGES(asm)
	assert len(asm.sj) == 4959
	assert len(asm.me) == 6039
	f()
	LOG.info('{0},{1}'.format(len(asm.me),len(asm.sj)))
	assert N.sum(~asm.me['strand'].isin(['+','-','.']))==0
	assert len(asm.me) == 6039
	assert len(asm.sj) == 4959
	

def test_fixstrand(asm):
	f = AS.FIXSTRAND(asm)
	assert len(asm.sj) == 4959
	assert len(asm.me) == 6039
	f()
	LOG.info('{0},{1}'.format(len(asm.me),len(asm.sj)))
	assert N.sum(~asm.me['strand'].isin(['+','-','.']))==0
	assert len(asm.me) == 6039
	assert len(asm.sj) == 4959
	

def test_findirets(asm):
	f = AS.FINDIRETS(asm)
	assert len(asm.sj) == 4959
	assert len(asm.me) == 6039
	f()
	LOG.info('{0},{1}'.format(len(asm.me),len(asm.sj)))
	assert N.sum(~asm.me['strand'].isin(['+','-','.']))==0
	assert len(asm.me) == 6157
	assert len(asm.sj) == 4959
	

def test_edgefixer(asm):
	f = AS.EDGEFIXER(asm)
	assert len(asm.sj) == 4959
	assert len(asm.me) == 6157
	f()
	LOG.info('{0},{1}'.format(len(asm.me),len(asm.sj)))
	assert N.sum(~asm.me['strand'].isin(['+','-','.']))==0
	assert len(asm.me) == 5753
	assert len(asm.sj) == 4959
	

def test_findsecovth(asm):
	f = AS.FINDSECOVTH(asm)
	assert len(asm.sj) == 4959
	assert len(asm.me) == 5753
	f()
	LOG.info('{0},{1}'.format(len(asm.me),len(asm.sj)))
	assert N.sum(~asm.me['strand'].isin(['+','-','.']))==0
	assert len(asm.me) == 5753
	assert len(asm.sj) == 4959
	assert abs(asm.secovth - 0.4999999999999999)<1e-6
	

def test_findse(asm):
	f = AS.FINDSE(asm)
	assert len(asm.sj) == 4959
	assert len(asm.me) == 5753
	f()
	LOG.info('{0},{1},{2},{3}'.format(len(asm.me),len(asm.sj),len(asm.ae),len(asm.se)))
	assert N.sum(~asm.sj['strand'].isin(['+','-','.']))==0
	assert N.sum(~asm.ae['strand'].isin(['+','-','.']))==0
	assert len(asm.me) == 6083
	assert len(asm.sj) == 4959
	assert len(asm.ae) == 11055
	assert len(asm.se) == 4972
	

def test_find53ir(asm):
	f = AS.FIND53IR(asm)
	assert len(asm.sj) == 4959
	assert len(asm.me) == 6083
	f()
	LOG.info('{0},{1},{2},{3}'.format(len(asm.me),len(asm.sj),len(asm.ae),len(asm.se)))
	assert N.sum(~asm.ae['strand'].isin(['+','-','.']))==0
	assert len(asm.ae) == 7359 #7382
	assert len(asm.sj) == 4959
	

def test_calccov(asm):
	f = AS.CALCCOV(asm)
	f()
	LOG.info('{0},{1},{2},{3}'.format(len(asm.me),len(asm.sj),len(asm.ae),len(asm.se)))
	assert N.sum(~asm.ae['strand'].isin(['+','-','.']))==0
	assert 'cov' in asm.ae.columns
	

def test_setinfo(asm):
	f = AS.SETINFO(asm)
	f()
	LOG.info('{0},{1},{2},{3}'.format(len(asm.me),len(asm.sj),len(asm.ae),len(asm.se)))
	assert N.sum(~asm.ae['strand'].isin(['+','-','.']))==0
	assert 'd_id' in asm.ae.columns
	assert 'cat' in asm.ae.columns
	

def test_findgenes(asm):
	f = AS.FINDGENES(asm)
	f()
	LOG.info('{0},{1},{2},{3}'.format(len(asm.me),len(asm.sj),len(asm.ae),len(asm.se)))
	assert N.sum(~asm.ae['strand'].isin(['+','-','.']))==0
	assert '_gidx' in asm.ae.columns
	assert 'gname' in asm.ae.columns
	assert len(asm.genes) == 2008
	

def test_selectseme(asm):
	f = AS.SELECTSEME(asm)
	UT.set_exon_category(asm.sj,asm.ae)
	f()
	LOG.info('{0},{1},{2},{3}'.format(len(asm.me),len(asm.sj),len(asm.ae),len(asm.se)))
	assert len(asm.ae) == 6169 #6153
	

def test_fixedges2(asm):
	f = AS.FIXEDGES2(asm)
	f()
	LOG.info('{0},{1},{2},{3}'.format(len(asm.me),len(asm.sj),len(asm.ae),len(asm.se)))
	assert len(asm.ae) == 6338 #6311
	

def test_writesjex(asm):
	f = AS.WRITESJEX(asm)
	f()
	assert os.path.exists(asm.fnobj.txtname('sj'))
	assert os.path.exists(asm.fnobj.txtname('ex'))
	

def test_writegenes(asm):
	f = AS.WRITEGENES(asm)
	assert GBED not in asm.fnobj._fnames['temp']
	assert GBED not in asm.fnobj._fnames['output']
	f()
	assert os.path.exists(asm.fnobj.bedname('genes', category='read'))
	assert GBED not in asm.fnobj._fnames['temp']
	assert GBED in asm.fnobj._fnames['output']

def test_delete(asm):
	assert os.path.exists(asm.fnobj.txtname('sj'))==True
	asm.delete_intermediates()
	assert os.path.exists(asm.fnobj.txtname('sj'))==True
	asm.fnobj.delete(dcats=['output'])
	assert os.path.exists(asm.fnobj.txtname('sj'))==False

def test_assembler1(fnobj):
	asm1 = AS.Assembler(fnobj, False, False)
	asm1.assemble()
	assert os.path.exists(asm1.fnobj.bedname('genes', category='read'))

def test_assembler0(fnobj):
	# instantiation (__init__)
	asm1 = AS.Assembler(fnobj, False, False)
	fname = fnobj.txtname('assemble.params')
	if os.path.exists(fname):
		os.unlink(fname)
	assert asm1.params['merging'] == False
	# check_params
	asm1.check_params()
	assert asm1.params['checksjsupport'] == False
	assert asm1.params['override'] == False
	assert os.path.exists(fname) == True

	asm2 = AS.Assembler(fnobj, True, True)
	assert asm2.params['merging'] == True
	asm2.check_params()
	assert asm2.params['checksjsupport'] == True
	assert asm2.params['override'] == True

	asm3 = AS.Assembler(fnobj, False, True, mpth=0.5, binth=0.1)
	assert asm3.params['mpth'] == 0.5
	assert asm3.params['binth'] == 0.1
	asm3.check_params()
	assert asm3.params['checksjsupport'] == True
	assert asm3.params['override'] == True
	asm3.params['newparam'] = 1
	asm3.check_params()
	assert asm3.params['checksjsupport'] == True
	assert asm3.params['override'] == True







