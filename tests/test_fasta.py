
import os
import pytest

from jgem import fasta as FA


def test_GenomeFASTAChroms(datadir):
	chromdir = os.path.join(datadir, 'FASTA')
	ga = FA.GenomeFASTAChroms(chromdir, '\n')
	assert set(ga.chromosomes) == set(['chr1','chr2'])
	assert ga.get('chr1',100,110) == 'tcccagatga'
	assert ga.get('chr2',605,610) == 'ATGTC'
	assert ga.ext == '.fa'
	chromdir = os.path.join(datadir, 'FASTA.gz')
	ga = FA.GenomeFASTAChroms(chromdir, '\n')
	assert ga.ext == '.fa.gz'
	assert set(ga.chromosomes) == set(['chr1','chr2'])
	assert ga.get('chr1',100,110) == 'tcccagatga'
	assert ga.get('chr2',605,610) == 'ATGTC'





