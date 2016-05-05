import os, time
import pytest
import pandas as PD
import numpy as N
try:
    from StringIO import StringIO
except:
    from io import StringIO

from jgem import repeats as RP
from jgem import fasta as FA
from jgem import utils as UT

def test_count_repeats(datadir):
	TESTDATA=StringIO("""chr,st,ed,name,sc1,strand
chr1,0,10,a,0,+
chr1,5,20,a,1,-
chr1,25,30,a,1,+
chr1,40,45,b,2,-
chr1,45,50,b,2,+
chr1,49,55,c,2,+
chr2,255,260,d,3,-
chr2,260,270,d,4,+
chr2,370,380,e,4,-
chr2,380,390,e,5,+
	""")
	df = PD.DataFrame.from_csv(TESTDATA, sep=",", index_col=False)
	print(df)
	udf = UT.make_union_gene_bed(df, gidx='name')
	chromdir = os.path.join(datadir, 'FASTA')
	gfc = FA.GenomeFASTAChroms(chromdir)
	RP.count_repeats(udf, gfc, returnseq=True)
	print(udf)
	assert list(udf['#repbp']) == [0,0,0,3,15,20]
	assert udf.iloc[3]['seq'] == 'TCTgag'
	udf = UT.make_union_gene_bed(df, gidx='sc1')
	RP.count_repeats(udf, gfc, returnseq=False)
	print(udf)
	assert list(udf['#repbp']) == [0,0,0,3,5,10,10,10]
	