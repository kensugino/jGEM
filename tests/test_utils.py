import os, time
import pytest
import pandas as PD
import numpy as N

from jgem import utils as UT

def test_calc_locus():
	df =  PD.DataFrame({'chr':['chr1','chr2','chr3'],'st':[0,1,2],'ed':[10,20,30]})
	l = UT.calc_locus(df)
	o = ['chr1:0-10','chr2:1-20','chr3:2-30']
	assert o==list(l.values)

def test_calc_locus_strand():
	df =  PD.DataFrame({'chr':['chr1','chr2','chr3'],'st':[0,1,2],'ed':[10,20,30],'strand':['+','-','+']})
	l = UT.calc_locus_strand(df)
	o = ['chr1:0-10:+','chr2:1-20:-','chr3:2-30:+']
	assert o==list(l.values)

def test_makedirs(tmpdir):
	path = os.path.join(str(tmpdir), 'a/b/c')
	UT.makedirs(path)
	assert os.path.exists(path)
	# should not raise
	UT.makedirs(path)
	
	# make a file
	path2 = os.path.join(str(tmpdir), 'a/b/c/d')
	open(path2,'w').write('test\n')
	# should raise
	with pytest.raises(OSError):
		UT.makedirs(path2)

def test_notstale(outdir):
	a = os.path.join(outdir, 'a')
	af = open(a,'w').write('a')
	b = os.path.join(outdir, 'b')
	open(b,'w').write('b')
	c = os.path.join(outdir, 'c')
	open(c,'w').write('c')
	d = os.path.join(outdir, 'd')
	# simple: a < b
	assert UT.notstale(a, b)==True
	# multiple: [a,b] < c
	assert UT.notstale([a,b], c)==True
	# non-existent cache
	assert UT.notstale(a, d)==False
	# non-existent input
	# with pytest.raises(RuntimeError): # changed
	# 	UT.notstale(d, a)

def test_exid():
	ex0 = ['chr1',0,1,'name',0,'+']
	ex1 = {'chr':'chr1','st':0,'ed':1,'name':'name','sc1':0,'strand':'+'}
	eid = 'chr1:0-1:+'
	assert UT.exid(ex0)==eid
	assert UT.exid(ex1)==eid

def test_sjid():
	ex0 = ['chr1',0,1,'name',0,'+']
	ex1 = {'chr':'chr1','st':0,'ed':1,'name':'name','sc1':0,'strand':'+'}
	eid = 'chr1:0^1:+'
	assert UT.sjid(ex0)==eid
	assert UT.sjid(ex1)==eid

def test_set_ids():
	df = PD.DataFrame({'chr':['chr1']*10, 'st':[2,5,4,6,7,8,9,1,0,3],'ed':[1,2,3,4,5,6,7,8,9,0]})
	UT.set_ids(df)
	assert '_id' in df.columns
	assert all(df['_id'] == N.arange(10))
	assert all(df.index == N.arange(10))


# def test_set_ad_info():
# 	pass

# def test_set_pos_info():
# 	pass

# def test_find_nullidx():
# 	pass

# def test_set_exon_category():
# 	pass

# def test_set_info():
# 	pass

# def test_me_se():
# 	pass

def test_chromsizes():
	assert UT.chroms('mm10') == ['chr{0}'.format(i+1) for i in range(19)]+['chrX','chrY']
	assert UT.chroms('dm3') == ['chr2L','chr2LHet','chr2R','chr2RHet','chr3L','chr3LHet',
								'chr3R','chr3RHet','chr4','chrX','chrXHet','chrYHet',
								'chrU','chrUextra']
	assert os.path.exists(UT.chromsizes('mm10'))
	assert os.path.exists(UT.chromsizes('dm3'))
	assert os.path.exists(UT.chromsizes('hg19'))
