
import os
import pytest

from jgem import bedtools as BT
import gzip


def test_bedtoolintersect(tmpdir):
	adata = """chr1	10	20
chr1	30	40
"""

	bdata = """chr1	15	20
"""

	odata = """chr1	10	20
"""

	a = tmpdir.join('a.bed')
	a.write(adata)
	b = tmpdir.join('b.bed')
	b.write(bdata)
	c = tmpdir.join('c.bed')
	c1 = BT.bedtoolintersect(str(a),str(b),str(c), wa=True)
	print(c1)
	assert c1==str(c)+'.gz'
	cdata = gzip.open(c1).read().decode('ascii')
	assert cdata == odata

def test_calcovlratio(tmpdir):
	adata = """chr1	10	20
chr1	30	40
"""

	bdata = """chr1	15	20
chr1	38	48
"""

	odata = """chr	st	ed	len	ovl	ovlratio	notcovbp
chr1	10	20	10	5	0.5	5
chr1	30	40	10	2	0.2	8
"""
	a = tmpdir.join('a.bed')
	a.write(adata)
	b = tmpdir.join('b.bed')
	b.write(bdata)
	c = tmpdir.join('c.bed.gz')
	ovl = BT.calc_ovlratio(str(a),str(b),str(c),3,3,idcol=['chr','st','ed'])
	cdata = gzip.open(str(c)).read().decode('ascii')
	assert odata == cdata

def test_fillgap(tmpdir):
	adata = """chr1	10	20
chr1	30	40
chr1	100	110
chr1	120	150
"""
	odata = """chr1	10	40
chr1	100	150
"""

	a = tmpdir.join('a.bed')
	b = tmpdir.join('b.bed')
	a.write(adata)
	c = BT.fillgap(str(a),str(b), gap=50)
	cdata = gzip.open(str(c)).read().decode('ascii')
	assert odata == cdata

