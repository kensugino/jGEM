import os, time
import pytest
import pandas as PD
import numpy as N
try:
    from StringIO import StringIO
except:
    from io import StringIO

from jgem import utils as UT



def test_make_union_gene_bed():
	TESTDATA=StringIO("""chr,st,ed,name,sc1,strand
chr1,0,10,a,0,+
chr1,5,20,a,1,-
chr1,25,30,a,1,+
chr1,40,45,b,2,-
chr1,45,50,b,2,+
chr1,49,55,c,2,+
chr2,55,60,d,3,-
chr2,60,70,d,4,+
chr2,70,80,e,4,-
chr2,80,90,e,5,+
	""")
	df = PD.DataFrame.from_csv(TESTDATA, sep=",", index_col=False)
	print(df)
	udf1 = UT.make_unionex(df, gidx='name')
	print(udf1)
	assert len(udf1) == 6
	assert all(udf1.columns == ['chr','st','ed','name','sc1','strand'])
	assert list(udf1.iloc[0]) == ['chr1',0,20,'a',0,'+']
	assert list(udf1.iloc[-1]) == ['chr2',70,90,'e',4,'-']
	udf2 = UT.make_unionex(df, gidx='sc1')
	print(udf2)
	assert len(udf2) == 7


def test_union_contiguous():
	TESTDATA=StringIO("""chr,st,ed,name,sc1,strand
chr1,0,10,a,0,+
chr1,5,20,a,0,-
chr1,25,30,a,0,+
chr1,40,45,a,0,-
chr1,45,50,a,0,+
chr1,49,55,a,0,+
chr2,55,60,a,0,-
chr2,60,70,a,0,+
chr2,70,80,a,0,-
chr2,80,90,a,0,+
	""")

	df = PD.DataFrame.from_csv(TESTDATA, sep=",", index_col=False)
	print(df)
	udf = UT.union_contiguous(df)
	print(udf)
	assert len(udf) == 4
	assert all(udf.columns == ['chr','st','ed','name','sc1','strand'])
	assert list(udf.iloc[0]) == ['chr1',0,20,'a',0,'+']
	assert list(udf.iloc[-1]) == ['chr2',55,90,'a',0,'-']



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
