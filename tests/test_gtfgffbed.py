import os
import pytest
import pandas as PD

from jgem import utils as UT
from jgem import gtfgffbed as GGB


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
	assert len(gtf) == 50942


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

	