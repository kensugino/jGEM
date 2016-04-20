import os
import pytest
import pandas as PD
import numpy as N
import subprocess

from jgem import utils as UT
from jgem import bigwig as BW

# def test_cnt_bam():
# 	"needs samtools"
# 	pass

# def test_bam2bw():
# 	"needs genomeCoverageBed and wigToBigWig"
# 	pass

def test_apply_threshold(bigwig):
	it = BW.apply_threshold(bigwig, 1, ['chr1'])
	for i in range(10):
		x = next(it)
	assert x == ('chr1', 3203918, 3203926)

def test_block_iter(bigwig):
	it = BW.block_iter(bigwig, 'chr1', chunk=int(1e6))
	for i in range(10):
		x = next(it)
	assert x == (3118639, 3118733, 0.05421944335103035)

# def test_bw2bed():
# 	pass

def test_get_bigwig_as_array(bigwig):
	a = BW.get_bigwig_as_array(bigwig, 'chr1', 3.09e6, 3.10e6)
	assert N.sum(a>0.04) == 3

# def test_merge_bigwigs_chr():
# 	pass

def test_array2wiggle_chr(bigwig, tmpdir, mm10chromsizes):
	b = BW.get_bigwig_as_array(bigwig, 'chr1', 0, 4e6)
	testwig = os.path.join(str(tmpdir), 'test2.wig')
	testbw = os.path.join(str(tmpdir), 'test2.bw')
	BW.array2wiggle_chr(b, 'chr1', testwig)
	cmd = ['wigToBigWig', testwig, mm10chromsizes, testbw]
	subprocess.call(cmd)
	assert os.path.exists(testbw)
	c = BW.get_bigwig_as_array(testbw, 'chr1', 0, 4e6)
	assert all(b==c)

# def test_merge_bigwigs_mp():
# 	pass

