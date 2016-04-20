import os
import pytest
import pandas as PD
import numpy as N

from jgem import utils as UT
from jgem import trimex as TE


def test_trim_ex_worker(g4sjex):
	sj0,ex0 = g4sjex
	cols = list(ex0.columns)
	if 'len' not in cols:
		ex0['len'] = ex0['ed']-ex0['st']
	gidxs = ex0['gene_id'].unique()[:20]
	ex = ex0[ex0['gene_id'].isin(gidxs)] 

	tex = PD.DataFrame(TE.trim_ex_worker((ex,1000,'gene_id')),columns=cols+['len'])
	tex['len2'] = tex['ed']-tex['st']
	#g2s = UT.series2dict(tex.groupby('gene_id')['len2'].sum())
	#tex['glen'] = [g2s[x] for x in tex['gene_id']]
	ex = ex.set_index('gene_id')
	tex = tex.set_index('gene_id')

	# single exon <= 1000bp no change
	gid = 'ENSMUSG00000064842.1' # 110 bp
	assert len(ex.ix[[gid]])==1
	assert len(tex.ix[[gid]])==1
	assert tex.ix[gid]['len2']==ex.ix[gid]['len']

	# single exon > 1000bp
	gid = 'ENSMUSG00000102348.1' # 1634 bp
	assert len(ex.ix[[gid]])==1
	assert len(tex.ix[[gid]])==1
	assert ex.ix[gid]['len']>1000
	assert tex.ix[gid]['len2']==1000

	# multiple exon <= 1000bp
	gid = 'ENSMUSG00000089699.1' # 149+101 = 250 bp
	assert len(ex.ix[[gid]])>1
	assert ex.ix[gid]['len'].sum()<1000
	assert len(tex.ix[[gid]])>1
	assert tex.ix[gid]['len2'].sum()==ex.ix[gid]['len'].sum()

	# multipe exon > 1000bp but 3'-most UTR > 1000 bp
	gid = 'ENSMUSG00000025902.9'
	assert len(ex.ix[[gid]])>1
	assert len(tex.ix[[gid]])==1
	assert ex.ix[gid]['len'].sum()>1000
	assert tex.ix[gid]['len2']==1000

	# multiple exons > 1000bp 3'-most UTR < 1000 bp
	gid = 'ENSMUSG00000102948.1' # 194+1673 => 194+806
	assert len(ex.ix[[gid]])>1
	assert len(tex.ix[[gid]])>1
	assert ex.ix[gid]['len'].sum()>1000
	assert tex.ix[gid]['len2'].sum()==1000


def test_trim_ex(g4sjex, tmpdir):
	sj0,ex0 = g4sjex
	cols = list(ex0.columns)
	if 'len' not in cols:
		ex0['len'] = ex0['ed']-ex0['st']
	gidxs = ex0['gene_id'].unique()[:20]
	ex = ex0[ex0['gene_id'].isin(gidxs)].copy()
	ex.loc[ex['gene_id'].isin(gidxs[10:]),'chr'] = 'chr2'
	expath = os.path.join(str(tmpdir),'ex.txt.gz')
	UT.write_pandas(ex,expath,'h')
	dstpath = os.path.join(str(tmpdir),'tex.txt.gz')
	dstcipath = os.path.join(str(tmpdir),'texci.txt.gz')
	tex = TE.trim_ex(expath,dstpath,dstcipath,1000,'gene_id',2)
	assert len(ex) == 58
	assert len(tex) == 25
	assert os.path.exists(dstpath)
	assert os.path.exists(dstcipath)




