"""

.. module:: repeats
    :synopsis:  Repeats (transposon) related stuffs

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""

import csv
import subprocess
import os
import gzip
import glob
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

import pandas as PD
import numpy as N

from jgem import utils as UT
from jgem import fasta as FA



def count_repeats(beddf, genomefastaobj, col='#repbp', returnseq=False, seqcol='seq'):
	"""Looks up genome sequence and counts the number of lower characters.
	(RepeatMaker masked sequence are set to lower characters in UCSC genome)

	Args:
		beddf: Pandas DataFrame with chr,st,ed frame, when calculating repeats bp
		 for genes, unioned bed should be used (use utils.make_union_gene_bed)
		 chr,st,ed columns are required.
		genomefastaobj: an object with get(chr,st,ed) method that returns sequence
		 (use fasta.GenomeFASTAChroms).
		col: column names where counts will be put in
		returnseq (bool): whether to return sequence or not (default False)
		seqcol: column where sequences are put in (default seq)

	Outputs:
		are put into beddf columns with colname col(default #repbp)

	"""
	def _cnt(chrom,st,ed):
		seq = genomefastaobj.get(chrom,st,ed)
		return N.sum([x.islower() for x in seq])

	if returnseq:
		beddf[seqcol] = [genomefastaobj.get(*x) for x in beddf[['chr','st','ed']].values]
		beddf[col] = beddf[seqcol].apply(lambda x: N.sum([y.islower() for y in x]))
	else:
		beddf[col] = [_cnt(*x) for x in beddf[['chr','st','ed']].values]
	return beddf

