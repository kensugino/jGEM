import os
import pytest
import pandas as PD
try:
    from StringIO import StringIO
except:
    from io import StringIO

from jgem import utils as UT
from jgem import gtfgffbed as GGB
import jgem.cy.bw as cybw

def test_unionex2bed12():
	TESTDATA=StringIO("""chr,st,ed,name,sc1,strand
chr1,0,10,a,0,+
chr1,15,20,a,1,-
chr1,25,30,a,1,+
chr1,40,45,b,2,-
chr1,47,50,b,2,+
chr1,49,55,c,2,+
chr2,58,60,d,3,-
chr2,65,70,d,4,+
chr2,70,80,e,4,-
chr2,85,90,e,5,+
	""")
	df = PD.DataFrame.from_csv(TESTDATA, sep=",", index_col=False)
	df['_gidx'] = df['name']
	print(df)
	b12 = GGB.unionex2bed12(df)
	print(b12)
	# assert 0


def test_read_gtf_cy(g4gtfpath):
	#gtf = GGB.read_gtf(g4gtfpath)
	recs, cols = cybw.read_gtf_helper(g4gtfpath, GGB.DEFAULT_GTF_PARSE)
	gtf = PD.DataFrame(recs, columns=cols)
	assert len(gtf) == 112665


def test_sjtab2sjbed(sampleinfo, datadir, outdir):
	rec = sampleinfo.iloc[0]
	sjtab = os.path.join(datadir, 'SJ', rec['sjtab'])
	sjbed = os.path.join(outdir, rec['sjbed'])
	aligned = rec['aligned']
	sj = GGB.sjtab2sjbed(sjtab,sjbed,aligned)
	assert os.path.exists(sjbed)
	SJCOLS = ['chr','st','ed','strand2','motif','annotated','ureads','mreads','maxoverhang']
	sji = PD.read_table(sjtab, names=SJCOLS)
	assert len(sj)==len(sji)
	#cols = ['chr','st','ed','name','strand','ucnt','mcnt']
	#sjo = PD.read_table(sjbed, compression='gzip', names=cols)
	sjo = GGB.read_bed(sjbed)
	assert all(sj[GGB.BEDCOLS[:7]] == sjo)
	#assert len(sjo)==len(sji)
	#assert os.path.getsize(sjbed) == 114793

def test_read_gtf(g4gtfpath):
	gtf = GGB.read_gtf(g4gtfpath)
	assert len(gtf) == 112665


# def test_read_gff():
# 	pass

# def test_read_bed():
# 	pass

# def test_write_gtf():
# 	pass

# def test_write_gff():
# 	pass

# def test_write_bed():
# 	pass

# def test_gtf2gff():
# 	pass

# def test_gtf2bed12():
# 	pass

# def test_bed2gtf():
# 	pass

# def chop_chrs_gtf():
# 	pass

	