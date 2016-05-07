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


def test_count_repeats_viz_mp(outdir, testbed):
	TESTDATA=StringIO("""st,ed,name,sc1,chr,strand,_id
0,10,a,0,chr1,+,1
5,20,a,1,chr1,-,2
25,30,a,1,chr1,+,3
40,45,b,2,chr1,-,4
45,50,b,2,chr1,+,5
49,55,c,2,chr1,+,6
255,260,d,3,chr2,-,7
260,270,d,4,chr2,+,8
370,380,e,4,chr2,-,9
380,390,e,5,chr2,+,10
	""")
	TESTDATA2=StringIO("""st,ed,name,sc1,chr,strand
0,5,a1,0,chr1,+
9,20,a2,1,chr1,-
31,35,a3,1,chr1,+
40,45,b1,2,chr1,-
45,47,b2,2,chr1,+
56,70,c1,2,chr1,+
200,210,d1,3,chr2,-
260,280,d2,4,chr2,+
391,400,e1,4,chr2,-
	""")
	df = PD.DataFrame.from_csv(TESTDATA, sep=",", index_col=False)
	rmsk = PD.DataFrame.from_csv(TESTDATA2, sep=",", index_col=False)
	path = os.path.join(outdir, 'rmsktest.bed.gz')
	UT.write_pandas(rmsk[['chr','st','ed','name','sc1','strand']], path, '')
	print(df)
	print(rmsk)
	rslt = RP.count_repeats_viz_mp(df, path, expand=0)
	print(df)
	rslt = RP.count_repeats_viz_mp(df, path, expand=10)
	print(df)
	# assert 0





def test_count_repeats(datadir):
	TESTDATA=StringIO("""st,ed,name,sc1,chr,strand
0,10,a,0,chr1,+
5,20,a,1,chr1,-
25,30,a,1,chr1,+
40,45,b,2,chr1,-
45,50,b,2,chr1,+
49,55,c,2,chr1,+
255,260,d,3,chr2,-
260,270,d,4,chr2,+
370,380,e,4,chr2,-
380,390,e,5,chr2,+
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
	

def test_count_repeats_mp(datadir):
	TESTDATA=StringIO("""st,ed,name,sc1,chr,strand
0,10,a,0,chr1,+
5,20,a,1,chr1,-
25,30,a,1,chr1,+
40,45,b,2,chr1,-
45,50,b,2,chr1,+
49,55,c,2,chr1,+
255,260,d,3,chr2,-
260,270,d,4,chr2,+
370,380,e,4,chr2,-
380,390,e,5,chr2,+
	""")
	df = PD.DataFrame.from_csv(TESTDATA, sep=",", index_col=False)
	print(df)
	udf = UT.make_union_gene_bed(df, gidx='name')
	chromdir = os.path.join(datadir, 'FASTA')
	gfc = FA.GenomeFASTAChroms(chromdir)
	RP.count_repeats_mp(udf, gfc, returnseq=True)
	print(udf)
	assert list(udf['#repbp']) == [0,0,0,3,15,20]
	assert udf.iloc[3]['seq'] == 'TCTgag'
	udf = UT.make_union_gene_bed(df, gidx='sc1')
	udf = RP.count_repeats_mp(udf, gfc, returnseq=False)
	print(udf)
	assert list(udf['#repbp']) == [0,0,0,3,5,10,10,10]
