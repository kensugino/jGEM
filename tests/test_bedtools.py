
import os
import pytest

from jgem import bedtools as BT
from jgem import gtfgffbed as GGB
import gzip


# def test_mapbed2bw(testbed7, outdir):
# 	bwpre = os.path.join(outdir, 'mapbed2bw.test')
# 	BT.mapbed2bw(testbed7, bwpre,'mm10')

def test_get_total_bp(testbed, testbed12):
	bed = GGB.read_bed(testbed)
	r1 = BT.get_total_bp(bed)
	r2 = BT.get_total_bp_bedfile(testbed, bed12=False)
	r3 = BT.get_total_bp_bedfile(testbed12, bed12=True)
	assert r1 == r2
	assert r1 == r3
	assert r1 == (1641348, 405746)

def test_wig2bw(testwig, testchromsizes, outdir):
	opath = os.path.join(outdir, 'testwig.bw')
	if os.path.exists(opath):
		os.unlink(opath)
	BT.wig2bw(testwig, testchromsizes, opath)
	assert os.path.exists(opath)

def test_bam2bw(testbam, testchromsizes, outdir):
	opath = os.path.join(outdir, 'testbam.bw')
	if os.path.exists(opath):
		os.unlink(opath)
	BT.bam2bw(testbam, testchromsizes, opath)
	assert os.path.exists(opath)

def test_bed2bw(testbed, testchromsizes, outdir):
	opath = os.path.join(outdir, 'testbed.bw')
	if os.path.exists(opath):
		os.unlink(opath)
	BT.bed2bw(testbed, testchromsizes, opath)
	assert os.path.exists(opath)


def test_bam2wig(testbam, testchromsizes, outdir):
	opath = os.path.join(outdir, 'test3.wig.gz')
	if os.path.exists(opath):
		os.unlink(opath)
	ret = BT.bam2wig(testbam, testchromsizes, opath)
	assert os.path.exists(ret)
	assert opath==ret

	opath = os.path.join(outdir, 'test4.wig')
	if os.path.exists(opath):
		os.unlink(opath)
	ret = BT.bam2wig(testbam, testchromsizes, opath, scale=2.)
	assert os.path.exists(ret)
	assert opath==ret
	assert len(open(opath).readlines())==60352

def test_bed2wig(testbed, testbed12, testchromsizes, outdir):
	opath = os.path.join(outdir, 'test.wig.gz')
	if os.path.exists(opath):
		os.unlink(opath)
	ret = BT.bed2wig(testbed, testchromsizes, opath)
	assert os.path.exists(ret)
	assert opath==ret

	opath = os.path.join(outdir, 'test2.wig')
	if os.path.exists(opath):
		os.unlink(opath)
	ret = BT.bed2wig(testbed, testchromsizes, wigpath=opath, scale=2.)
	assert os.path.exists(ret)
	assert opath==ret
	assert len(open(opath).readlines())==60352

	opath = os.path.join(outdir, 'test.bed12.wig')
	if os.path.exists(opath):
		os.unlink(opath)
	ret = BT.bed2wig(testbed12, testchromsizes, wigpath=opath)
	assert os.path.exists(ret)
	assert opath==ret
	assert len(open(opath).readlines())==60352

def test_bam2bed12(testbam, outdir):
	opath = os.path.join(outdir, 'test.bed12.gz')
	if os.path.exists(opath):
		os.unlink(opath)
	ret = BT.bam2bed12(testbam, opath)
	assert os.path.exists(ret)
	assert opath==ret

	opath = os.path.join(outdir, 'test.bed12')
	if os.path.exists(opath):
		os.unlink(opath)
	ret = BT.bam2bed12(testbam, opath)
	assert os.path.exists(ret)
	assert opath==ret
	assert len(open(opath).readlines())==45593

def test_bam2bed(testbam, outdir):
	opath = os.path.join(outdir, 'test.bed.gz')
	if os.path.exists(opath):
		os.unlink(opath)
	ret = BT.bam2bed(testbam, opath)
	assert os.path.exists(ret)
	assert opath==ret
	opath = os.path.join(outdir, 'test2.bed')
	if os.path.exists(opath):
		os.unlink(opath)
	ret = BT.bam2bed(testbam, opath)
	assert os.path.exists(ret)
	assert opath==ret
	assert len(open(opath).readlines())==46624


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
	assert c1==str(c)
	cdata = open(c1).read()
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
	cdata = open(str(c)).read()
	assert odata == cdata

