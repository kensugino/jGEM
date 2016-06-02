"""

.. module:: bedtools
    :synopsis: interface to bedtools

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
import os
import subprocess
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)
import gzip
from collections import defaultdict
from collections import Counter
import glob
import shutil

import gzip
import csv


import pandas as PD
import numpy as N

from jgem import utils as UT
from jgem import gtfgffbed as GGB
from jgem import fasta as FA
import jgem.cy.bw as cybw

import inspect


### BAM,BED,WIGGLE,BIGWIG ##############################################

# decorators to separate logic
def compressQ(outname, noerr=0):
    """ decorator for checking file compression and error """
    def deco(func):
        argnames, varargs, keywords, defaults = inspect.getargspec(func)
        pos = argnames.index(outname)
        def wrap(*args,**kwargs):
            # check output '.gz'
            if outname in kwargs:
                opath = kwargs[outname]
            else:
                opath = args[pos]
                args = list(args)
            if opath[-3:]=='.gz':
                compress = True
                opath = opath[:-3]
            else:
                compress = False
            UT.makedirs(os.path.dirname(opath))
            if outname in kwargs:
                kwargs[outname] = opath
            else:
                args[pos] = opath
            err = func(*args, **kwargs)
            if err != noerr:
                LOG.warning('bederror:{0}, err={1}'.format(func.__name__, err))
                raise RuntimeError(func.__name__)
            if compress:
                return UT.compress(opath)
            return opath
        return wrap
    return deco

def logerr(noerr=0):
    def deco(func):
        def wrap(*args, **kwargs):
            err = func(*args, **kwargs)
            if err != noerr:
                LOG.warning('bederror:{0}, err={1}'.format(func.__name__, err))
                raise RuntimeError(func.__name__)
            return err
        return wrap
    return deco

@compressQ('bedpath', None)
def bam2bed(bampath, bedpath):
    """Convert BAM to BED7

    BED name field (column 4) contains read id (so that together with map id (col 7) multi-mapper can be identified)
    BED tst field (column 7) contains map id (so that split reads can be identified)

    BED sc1 field (column 5) is from bedtools bamtobed and contains mapping quality
    """

    cmd1 = ['bedtools','bamtobed','-i', bampath, '-split','-bed12']
    awkscript = 'BEGIN{OFS="\t";c=1;}{ n=split($11,a,","); n=split($12,b,","); for(i=1;i<=n;i++){st=$2+b[i]; print $1,st,st+a[i],$4,$5,$6,NR}}'    
    # above keep the original name so that you can go back to fastq
    # awkscript = 'BEGIN{OFS="\t";c=1;}{if(d[$4]){$4=d[$4];}else{d[$4]=c;$4=c;c++;} n=split($11,a,","); n=split($12,b,","); for(i=1;i<=n;i++){st=$2+b[i]; print $1,st,st+a[i],$4,$5,$6,NR}}'
    cmd2 = ['awk',awkscript]
    with open(bedpath, 'wb') as fp:
        p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=fp)
        err = p2.communicate()[1]
    return err



@compressQ('bedpath', None)
def bam2bed12(bampath, bedpath):
    """Convert BAM to BED12

    BED name field (column 4) contains read id (so that multi-mapper can be identified)
    BED tst field (column 7) contains map id 
    BED sc1 field (column 5) is from bedtools bamtobed and contains mapping quality
    """

    cmd1 = ['bedtools','bamtobed','-i', bampath, '-split', '-bed12']
    awkscript = 'BEGIN{OFS="\t";c=1;}{$7=NR; print $0;}'
    #awkscript = 'BEGIN{OFS="\t";c=1;}{if(a[$4]){$4=a[$4];}else{a[$4]=c;$4=c;c++;}; $7=NR; print $0;}'
    cmd2 = ['awk',awkscript]
    with open(bedpath, 'wb') as fp:
        p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(cmd2, stdin=p1.stdout, stdout=fp)
        err = p2.communicate()[1]
    return err

@compressQ('wigpath')
def bed2wig(bedpath, chromsizes, wigpath, scale=None):
    """Runs BEDTOOLS genomecov. Takes BED, makes WIGGLE."""
    if scale is None:
        cmd1 = ['bedtools','genomecov','-bg', '-split', '-i', bedpath, '-g', chromsizes]
    else:
        cmd1 = ['bedtools','genomecov','-bg', '-split', '-i', bedpath, '-g', chromsizes, '-scale', str(scale)]
    with open(wigpath,'wb') as fp:
        p1 = subprocess.Popen(cmd1, stdout=fp)
        err = p1.wait()
    return err

@compressQ('wigpath')
def bam2wig(bampath, chromsizes, wigpath, scale=None):
    """Runs BEDTOOLS genomecov. Takes BAM, makes WIGGLE."""
    if scale is None:
        cmd1 = ['bedtools', 'genomecov', '-split', '-bg', '-ibam', bampath, '-g', chromsizes]
    else:
        cmd1 = ['bedtools', 'genomecov', '-split', '-bg', '-ibam', bampath, '-g', chromsizes, '-scale', str(scale)]
    with open(wigpath,'wb') as fp:
        p1 = subprocess.Popen(cmd1, stdout=fp)
        err = p1.wait()
    return err

@logerr(0)
def wig2bw(wigpath, chromsizes, bwpath):
    """Generate bigwig coverage from WIGGLE.
    Runs Kent's tool wigToBigWig.
    """
    cmd = ['wigToBigWig', wigpath, chromsizes, bwpath]
    UT.makedirs(os.path.dirname(bwpath))
    err = subprocess.call(cmd)
    return err

def bam2bw(bampath, chromsizes, bwpath, scale=None):
    """Generate bigwig coverage from BAM. """
    wigpath = bwpath+'.wig'
    bam2wig(bampath, chromsizes, wigpath, scale)
    wig2bw(wigpath, chromsizes, bwpath)
    os.unlink(wigpath)

def bed2bw(bedpath, chromsizes, bwpath, scale=None):
    """Generate bigwig coverage from BED. """
    wigpath = bwpath+'.wig'
    bed2wig(bedpath, chromsizes, wigpath, scale)
    wig2bw(wigpath, chromsizes, bwpath)
    os.unlink(wigpath)

def make_bw_from_bed0(bedpath, chromsizes, bwpath):
    """ DEPRECATED convert BED to BIGWIG, normalize average coverage to 1 """
    totbp,covbp = get_total_bp_bedfile(bedpath)
    scale = float(covbp)/totbp # 1/avgcov
    bed2bw(bedpath, chromsizes, bwpath, scale)

def make_bw_from_bam0(bampath, chromsizes, bedpath, bwpath):
    """ DEPRECATED convert BAM to BIGWIG, normalize average coverage to 1 """
    bam2bed(bampath, bedpath)
    make_bw_from_bed(bedpath, chromsizes, bwpath)
    
def bed12_bed6(bed):
    """ convert BED12 to BED6 uses cython helper """
    # BED12 ['chr', 'st', 'ed', 'name', 'sc1', 'strand', 'tst', 'ted', 'sc2', '#exons', 'esizes', 'estarts']
    # BED6 ['chr', 'st', 'ed', 'name', 'sc1', 'strand'] flatten exons, collect unique
    # BED12 tid => BED6 name=tid+exon_number
    bed8 = bed[['chr','st','ed','name','sc1','strand','esizes','estarts']]
    # def _gen():
    #     for x in df.values:
    #         esi = str(x[-2])
    #         est = str(x[-1])
    #         if esi[-1]==',':
    #             esi = esi[:-1]
    #             est = est[:-1]
    #         for y,z in zip(esi.split(','), est.split(',')):
    #             x[-2] = y
    #             x[-1] = z
    #             yield x
    # fbed = PD.DataFrame([x for x in _gen()], columns = df.columns)
    fbed = PD.DataFrame(cybw.flatten_bed8(bed8.values), columns=bed8.columns)
    # adjust st, ed
    fbed['st'] = fbed['st'] + fbed['estarts']
    fbed['ed'] = fbed['st'] + fbed['esizes']
    return fbed[['chr','st','ed','name','sc1','strand','esizes']]


@compressQ('bed6path')
def bed12ToBed6(bed12path, bed6path):
    """ uses bedtools bed12ToBed6 """
    cmd = ['bed12ToBed6', '-i', bed12path]
    with open(bed6path, 'wb') as fp:
        p1 = subprocess.Popen(cmd, stdout=fp)
        err = p1.wait()
    return err



### Normalization Scale ##############################################

def save_bed_covstats(bedpath, dstpath, bed12=False, checkuniq=False):
    tdic,cdic = get_total_bp_bedfile(bedpath, bed12, returndics=True, checkuniq=checkuniq)
    df = PD.DataFrame({c: {'totbp':tdic[c], 'covbp':cdic[c]} for c in cdic}).T
    df['acov'] = df['totbp']/df['covbp']
    df = df.sort_values('covbp',ascending=False)
    return UT.write_pandas(df, dstpath, 'ih')

def get_total_bp_bedfile(bedpath, bed12=False, chroms=None, returndics=False, checkuniq=False):
    """ Returns total mapped base pairs (totbp) and covered base pairs (covbp).
    The ratio totbp/covbp gives average coverage. Process without reading entire data
    into the RAM. Process non BED12 file.

    Args:
        bedpath (str): a path to BED file
        bed12 (bool): whether format is BED12 (default False)
        chroms (list): chromosomes to consider, if None (default) use all

    Returns:
        totbp: total base pairs in BED
        covbp: covered base pairs

    See:
        :py:func:`jgem.bigwig.get_totbp_covbp_bw ` (>6x faster if you have bigwig)

    """
    if bed12:
        totbpdic,covbpdic = cybw.get_total_bp_bed12file_helper(bedpath)
    else:
        if checkuniq:
            totbpdic,covbpdic = cybw.get_total_bp_bedfile_helper_check_uniq(bedpath)
        else:
            totbpdic,covbpdic = cybw.get_total_bp_bedfile_helper(bedpath)
    # fix key bytes => str
    tdic = {}
    cdic = {}
    for b in covbpdic.keys():
        u = b.decode()
        tdic[u] = totbpdic[b]
        cdic[u] = covbpdic[b]
    if returndics:
        return tdic, cdic    
    totbp = 0
    covbp = 0
    if chroms is None:
        chroms = cdic.keys()
    for chrom in chroms:
        if chrom not in cdic:
            LOG.warning('{0} not found in the data'.format(chrom))
            continue
        totbp += tdic[chrom]
        covbp += cdic[chrom]
    return totbp, covbp

def get_total_bp(beddf, returndics=False):
    """ Returns total mapped base pairs (totbp) and covered base pairs (covbp).
    The ratio totbp/covbp gives average coverage.

    Args:
        beddf: a BED dataframe (a standard non BED12 format)

    Returns:
        totbp: total base pairs in BED
        covbp: covered base pairs
    """
    # total bp
    beddf['len'] = beddf['ed']-beddf['st']
    # totbp = beddf['len'].sum() # total bp
    # covered bp: calculate using chopped intervals
    # first remove redundancy
    cols = ['st','ed']
    cdic = {}
    tdic = {}
    for chrom in beddf['chr'].unique():
        sub = beddf[beddf['chr']==chrom]
        tdic[chrom] = sub['len'].sum()
        sted = sub.groupby(cols).first().reset_index()
        a = N.array(sted[cols].sort_values(cols).values, dtype=N.int32)
        b = cybw.union_intervals(a)
        cdic[chrom] = N.sum([y-x for x,y in b])
    if returndics:
        return tdic, cdic
    totbp = N.sum(list(tdic.values()))
    covbp = N.sum(list(cdic.values()))
    return totbp, covbp


### INTERSECT ########################################################


def bedtoolintersect(aname, bname, cname, **kwargs):
    return _bedtoolscatcherror('intersect', aname, bname, cname, **kwargs)

def bedtoolmerge(aname, cname, **kwargs):
    return _bedtoolscatcherror2('merge',aname, cname, **kwargs)

def bedtoolcomplement(aname, cname, chromsizes):
    return _runbedtools2('complement',aname,cname,g=chromsizes)

def bedtoolsubtract(aname, bname, cname, **kwargs):
    return _bedtoolscatcherror('subtract', aname, bname, cname, **kwargs)

def _runbedtools2(which, aname, cname, **kwargs):
    cmd = ['bedtools',which, '-i', aname]
    for k,v in kwargs.items():
        if isinstance(v,bool):# in [True,False]: 2016-03-27 fix
            cmd += ['-'+k]
        else:
            cmd += ['-'+k, str(v)]
    with open(cname, "wb") as outfile:
        ret = subprocess.call(cmd, stdout=outfile)
    if ret!=0:
        msg = 'bederror return code:{0}, cmd:{1}'.format(ret, cmd)
        LOG.warning(msg)
        # delete output
        os.unlink(cname)
        raise RuntimeError(msg)
    return ret

def _runbedtools3(which, aname, bname, cname, **kwargs):
    cmd = ['bedtools',which,'-a',aname,'-b',bname]
    for k,v in kwargs.items():
        if isinstance(v,bool):# in [True,False]: 2016-03-27 fix
            cmd += ['-'+k]
        else:
            cmd += ['-'+k, str(v)]
    with open(cname, "wb") as outfile:
        ret = subprocess.call(cmd, stdout=outfile)
    if ret !=0:
        msg = 'bederror return code:{0}, cmd:{1}'.format(ret, cmd)
        LOG.warning(msg)
        # delete output
        os.unlink(cname)
        raise RuntimeError(msg)
    return ret

def _bedtoolscatcherror(which, aname, bname, cname, **kwargs):
    if not os.path.exists(aname):
        raise ValueError('{0} does not exists'.format(aname))
    if not os.path.exists(bname):
        raise ValueError('{0} does not exists'.format(bname))
        
    if cname.endswith('.gz'):
        cname = cname[:-3]
        compress=True
    else:
        compress=False
    try:
        ret = _runbedtools3(which,aname,bname,cname,**kwargs)
    except RuntimeError:
        LOG.warning('bedtool error: repeating on uncompressed a:{0},b:{1},c:{2}'.format(aname,bname,cname))
        aname2 = UT.uncompresscopy(aname)
        bname2 = UT.uncompresscopy(bname)
        ret = _runbedtools3(which,aname2,bname2,cname,**kwargs)
        if aname2 != aname:
            os.unlink(aname2)
        if bname2 != bname:
            os.unlink(bname2)
    if compress:
        return UT.compress(cname)
    return cname

def _bedtoolscatcherror2(which, aname, cname, **kwargs):
    if not os.path.exists(aname):
        raise ValueError('{0} does not exists'.format(aname))        
    if cname.endswith('.gz'):
        cname = cname[:-3]
        compress=True
    else:
        compress=False
    try:
        ret = _runbedtools2(which,aname,cname,**kwargs)
    except RuntimeError:
        LOG.warning('bedtool error: repeating on uncompressed a:{0},c:{1}'.format(aname,cname))
        aname2 = UT.uncompresscopy(aname)
        ret = _runbedtools2(which,aname2,cname,**kwargs)
        if aname2 != aname:
            os.unlink(aname2)
    if compress:
        return UT.compress(cname)
    return cname

def calc_ovlratio(aname, bname, tname, nacol, nbcol, idcol=['chr','st','ed'], returnbcols=False):
    """Calculate overlapped portion of b onto a. 
    Will check existence of result file (tname) and uses it if newer than input files.

    Args:
        aname (str): bed file name 1
        bname (str): bed file name 2
        tname (str): result file name
        nacol (int): number of columns in file 1
        nbcol (int): number of columns in file 2

    Optional:
        idcol (list of str): columns which specify unique entry

    Returns:
        A Pandas DataFrame which contains overlap info
    """
    # requirement: no overlap within b
    # cache?
    if UT.notstale([aname,bname], tname):
        return UT.read_pandas(tname)
    # calculate bedtools intersect
    tmpsuf='.ovlbed.txt'
    cname = aname+tmpsuf
    if nacol==12:
        cname = bedtoolintersect(aname, bname, cname, wao=True, split=True)
    else:
        cname = bedtoolintersect(aname, bname, cname, wao=True)
    # read tmp file
    acols = GGB.BEDCOLS[:nacol]
    bcols = ['b_'+x for x in GGB.BEDCOLS[:nbcol]]
    cols = acols + bcols +['ovl']
    df = UT.read_pandas(cname, names=cols)
    dfg = df.groupby(idcol) #['chr','st','ed'])
    if returnbcols:
        dfa = dfg.first().reset_index()[acols+bcols]
    else:
        dfa = dfg.first().reset_index()[acols]        
    if nacol==12:# sum of exon sizes
        dfa['len'] = [N.sum(map(int, x.split(',')[:-1])) for x in dfa['esizes']]
    else: 
        dfa['len'] = dfa['ed']-dfa['st']
    # since b does not overlap by itself total overlap of an element of a to b is 
    # sum of overlap to individual b
    dfa['ovl'] = dfg['ovl'].sum().values
    dfa['ovlratio'] = dfa['ovl'].astype(float)/dfa['len']
    dfa['notcovbp'] = dfa['len'] - dfa['ovl']
    # clean up
    os.unlink(cname)
    # save
    UT.save_tsv_nidx_whead(dfa, tname)
    return dfa

def fillgap(binfile, gapfile, gap=50):
    if gapfile[-3:]=='.gz':
        gapfile = gapfile[:-3]
    #gapfile = binfile[:-7]+'.gap%d.bed' % gap
    if UT.notstale(binfile, gapfile+'.gz'):
        return gapfile+'.gz'
    gapfile = bedtoolmerge(binfile, gapfile, d=gap)
    return gapfile

def read_ovl(c, acols, bcols=None):
    if bcols is None:
        cols = acols+['b_'+x for x in acols]+['ovl']
    else:
        cols = acols+['b_'+x for x in bcols]+['ovl']
    return UT.read_pandas(c, names=cols)


### MAPBED to WIG ########################################################

# dict read_id => set{map_id}
# multimapper = dup = size(set{map_id})>1
# weight = 1/dup
# 1st pass calculate this map read_id => weight
# for uniq.bw only use weight==1
# for all.bw use all but use weight


def splitbedgz(bedgz, prefix):
    """Split gzipped bed file into separate files according to chromosome. 
    Uses zcat and awk. 

    Args:
        bedgz: path to gzipped bed file
        prefix: output path prefix

    """
    cmd1 = ['zcat', bedgz]
    awkscript = 'BEGIN{{FS="\t"}}{{print > "{0}."$1".bed"}}'.format(prefix)
    #print(awkscript)
    cmd2 = ['awk', awkscript]
    p1 = subprocess.Popen(cmd1, stdout=subprocess.PIPE)
    p2 = subprocess.Popen(cmd2, stdin=p1.stdout)
    err = p2.communicate()[1]
    return err



# SJTABMOTIF = {0:'non-canonical',1:'GT/AG',2:'CT/AC',3:'GC/AG',4:'CT/GC',5:'AT/AC',6:'GT/AT'}
STED2STRAND = dict(
    GTAG='+',
    CTAC='-',
    GCAG='+',
    CTGC='-',
    ATAC='+',
    GTAT='-',
)

def _scan_make_map(paths, dstpath):
    cnt = defaultdict(set)
    #csp = defaultdict(int)
    for path in paths:
        if path[-3:]=='.gz':
            with gzip.open(path) as gz_file:
                with io.BufferedReader(gz_file) as fp:
                    for line in fp:
                        rec = line.strip().split(b'\t')
                        cnt[rec[3]].add(rec[6]) # for each read how many locations?
        else:
            with open(path,'rb') as fp:
                for line in fp: # chr,st,ed,name,sc1,strand,tst
                    rec = line.strip().split(b'\t') # read_id:name(3), map_id:tst(6)
                    cnt[rec[3]].add(rec[6]) # for each read how many locations?
                    # csp[rec[6]] += 1 # count # segments in a read if >1 spliced
    try:# py2
        dup = PD.DataFrame({k:len(v) for k,v in cnt.iteritems() if len(v)>1}, index=['cnt']).T
    except:
        dup = PD.DataFrame({k:len(v) for k,v in cnt.items() if len(v)>1}, index=['cnt']).T
    UT.write_pandas(dup, dstpath,'ih')
    
def pathcode(sse, strand):
    # sse: splice [(st,ed),...]
    if strand in ['+','.']:
        return ','.join(['{0}|{1}'.format(*x) for x in sse])
    return ','.join(['{1}|{0}'.format(*x) for x in sse[::-1]])

def pcode2pos(pcode):
    tmp = [[int(x) for x in y.split('|') for y in pcode.split(',')]]
    if tmp[0][0]<tmp[0][1]: # pos strand
        return tmp
    return [x[::-1] for x in tmp[::-1]]

def process_mapbed(bedpath, dstpre, genome, chromdir, stranded='.', np=10):
    """
    Args:
        bedpath: path to gzipped BED7 file (converted from BAM)
        dstpre: path prefix to destination
        genome: UCSC genome (mm10 etc.)
        chromdir: directory containing chromosome sequence in FASTA
        np: number of CPU to use

    Outputs:
        1. dstpre+'.ex.p.bw'
        2. dstpre+'.ex.n.bw'
        3. dstpre+'.ex.u.bw'
        4. dstpre+'.sj.p.bw'
        5. dstpre+'.sj.n.bw'
        6. dstpre+'.sj.u.bw'
        7. dstpre+'.ex.p.uniq.bw'
        8. dstpre+'.ex.n.uniq.bw'
        9. dstpre+'.ex.u.uniq.bw'
        10. dstpre+'.sj.p.uniq.bw'
        11. dstpre+'.sj.n.uniq.bw'
        12. dstpre+'.sj.u.uniq.bw'
        13. dstpre+'.sjpath.bed' BED12 (sc1:ucnt, sc2:jcnt=ucnt+mcnt)
    """
    chroms = UT.chroms(genome)
    chromdf = UT.chromdf(genome)
    chromsizes = UT.chromsizes(genome)

    # split into chroms
    UT.makedirs(dstpre)
    splitbedgz(bedpath, dstpre) # ~30sec
    duppath = dstpre+'.dupitems.txt.gz'
    chroms = [c for c in chroms if os.path.exists(dstpre+'.{0}.bed'.format(c))]
    files = [dstpre+'.{0}.bed'.format(c) for c in chroms]
    _scan_make_map(files, duppath)

    files0 = [dstpre+'.{0}.bed'.format(c) for c  in chromdf['chr'].values] # to be deleted
    args = [(dstpre, x, genome, chromdir, stranded) for x in chroms]
    # spread to CPUs
    rslts = UT.process_mp(_process_mapbed_chr, args, np=np, doreduce=False)
    # concatenate chr files
    files1 = []
    dstpath = dstpre+'.sjpath.bed'
    LOG.info('making {0}...'.format(dstpath))
    with open(dstpath, 'wb') as dst:
        for c in chroms:
            srcpath = dstpre+'.{0}.sjpath.bed'.format(c)
            files1.append(srcpath)
            with open(srcpath, 'rb') as src:
                shutil.copyfileobj(src, dst)
    dstpath = UT.compress(dstpath)

    for kind in ['.ex','.sj']:
        for strand in ['.p','.n','.u']:
            for suf in ['','.uniq']:
                pre = dstpre+kind+suf+strand
                wigpath = pre+'.wig'
                bwpath = pre+'.bw'
                with open(wigpath, 'wb') as dst:
                    for c in chroms:
                        srcpath = pre+'.{0}.wig'.format(c)
                        files1.append(srcpath)
                        if os.path.exists(srcpath):
                            with open(srcpath,'rb') as src:
                                shutil.copyfileobj(src, dst)
                LOG.info('making {0}...'.format(bwpath))
                if os.path.getsize(wigpath)>0:
                    wig2bw(wigpath, chromsizes, bwpath)
                files1.append(wigpath)

    # clean up temp files
    for x in files0+files1:
        if os.path.exists(x):
            LOG.debug('deleting {0}...'.format(x))
            os.unlink(x)

STRANDMAP = {('+','+'):'.p',
             ('+','-'):'.n',
             ('+','.'):'.u',
             ('-','+'):'.n',
             ('-','-'):'.p',
             ('-','.'):'.u',
             ('.','+'):'.u',
             ('.','-'):'.u',
             ('.','.'):'.u'}

def _process_mapbed_chr(dstpre, chrom, genome, chromdir, stranded):
    # 1st pass: calc dupdic
    bedpath = dstpre+'.{0}.bed'.format(chrom)
    dupids = UT.read_pandas(dstpre+'.dupitems.txt.gz', index_col=[0]).index
    # 2nd pass make wiggles
    gfc = FA.GenomeFASTAChroms(chromdir)
    chromsize = UT.df2dict(UT.chromdf(genome), 'chr', 'size')[chrom]
    
    # mqth MAPQ threshold there are ~6% <10
    # generator which makes an array
    fp = open(bedpath,'rb')

    wigs = {}
    wigpaths = {}
    for kind in ['.ex','.sj']:
        wigs[kind] = {}
        wigpaths[kind] = {}
        for strand in ['.p','.n','.u']:
            wigs[kind][strand] = {}
            wigpaths[kind][strand] = {}
            for suf in ['','.uniq']:
                wigpath = dstpre+kind+suf+strand+'.{0}.wig'.format(chrom)
                if os.path.exists(wigpath):
                    os.unlink(wigpath)
                wigpaths[kind][strand][suf] = wigpath
                wigs[kind][strand][suf] = N.zeros(chromsize, dtype=float)

    sjs = [] # path: (chr, st, ed, pcode, ucnt, strand, acnt)
    # pcode = a(apos)d(dpos) = a(ed)d(st) if strand=='+' else a(st)d(ed)
    # ucnt = unique read counts
    # acnt = multi-read adjusted all counts (=ucnt+Sum(mcnt(i)/dup(i)))
    # delete previous
    sjbed12 = dstpre+'.{0}.sjpath.bed'.format(chrom)
    if os.path.exists(sjbed12):
        os.unlink(sjbed12)

    def _write_arrays():
        for kind in ['.ex','.sj']:
            for strand in ['.p','.n','.u']:
                for suf in ['','.uniq']:
                    cybw.array2wiggle_chr64(wigs[kind][strand][suf], chrom,  wigpaths[kind][strand][suf], 'w')
        
    def _write_sj(sjs):
        # sjs = [(chr,st,ed,pathcode(name),ureads(sc1),strand,tst,ted,areads(sc2),cse),...]
        sjdf = PD.DataFrame(sjs, columns=GGB.BEDCOLS[:9]+['cse'])
        sjdfgr = sjdf.groupby('name')
        sj = sjdfgr.first()
        sj['sc1'] = sjdfgr['sc1'].sum().astype(int) # ucnt
        sj['sc2'] = sjdfgr['sc2'].sum().astype(int) # jcnt=ucnt+mcnt
        sj['st'] = sjdfgr['st'].min()
        sj['ed'] = sjdfgr['ed'].max()
        sj['#exons'] = sj['cse'].apply(len)+1
        sj['ests'] = [[0]+[z[1]-st for z in cse] for st,cse in sj[['st','cse']].values]
        sj['eeds'] = [[z[0]-st for z in cse]+[ed-st] for st,ed,cse in sj[['st','ed','cse']].values]
        esizes = [[u-v for u,v in zip(x,y)] for x,y in sj[['eeds','ests']].values]
        sj['estarts'] = ['{0},'.format(','.join([str(y) for y in x])) for x in sj['ests']]
        sj['esizes'] = ['{0},'.format(','.join([str(y) for y in x])) for x in esizes]
        sj = sj.reset_index()
        with open(sjbed12, 'w') as f:
            sj[GGB.BEDCOLS].to_csv(f, index=False, header=False, sep='\t', quoting=csv.QUOTE_NONE)
            
    def _append_sj(cse, css, csj, chrom,ureads,areads):
        if (len(cse)>0): # spits out splice rec
            # chr,st,ed,pathcode,ureads,strand,tst,ted,areads
            tst = cse[0][0]
            ted = cse[-1][1]
            if len(css)>0:
                strand = Counter(css).most_common()[0][0]
            else:
                strand = '.'
            name = pathcode(cse, strand)
            st = int(csj[0][1]) # first segment start
            ed = int(csj[-1][2]) # last segment end
            sjs.append((chrom,st,ed,name,ureads,strand,tst,ted,areads,cse))   
    
    def _add_to_ex_arrays(st,ed,dup,strand):
        kind='.ex'
        strand = STRANDMAP[(strand,stranded)]
        dic = wigs[kind][strand]
        dic[''][st:ed] += 1
        if not dup:
            dic['.uniq'][st:ed] += 1

    def _add_to_sj_arrays(sst,sed,dup,strand):
        kind='.sj'
        s = {'+':'.p','-':'.n','.':'.u'}[strand]
        dic = wigs[kind][s]
        # add to the arrays
        dic[''][sst:sed] += 1
        if not dup:
            dic['.uniq'][sst:sed] += 1
            ureads,areads = 1,1
        else:
            ureads,areads = 0,1
        return ureads,areads
        
    csj = [] # current collection of spliced reads
    css = [] # current strands
    cse = [] # current (sst,sed)
    csn = 0 # current segment number
    ureads,areads = 1,1 # uniq, total reads it's either 1,1 or 0,1
    pmid = None # previous map id common to spliced segments
    for line in fp:
        rec = line.strip().split(b'\t')
        # 7 column bed: chr(0), st(1), ed(2), name(3), mapq(4), strand(5), mapid(6)
        cchr = rec[0].decode()
        st,ed = int(rec[1]),int(rec[2])
        dup = rec[3] in dupids #dic[rec[3]]
        estrand = rec[5]
        _add_to_ex_arrays(st,ed,dup,estrand)
        # process splice
        if pmid != rec[6]: # new map 
            _append_sj(cse, css, csj, chrom, ureads, areads)
            csj,css,cse,csn = [rec],[],[],0 # reset running params
        else: # add segments
            csj.append(rec)            
            prec = csj[-2] # previous rec
            sst = int(prec[2]) # ed of previous segment
            sed = int(rec[1]) # st of current segment
            cse.append((sst,sed))
            # find strand
            sted = gfc.get(chrom,sst,sst+2)+gfc.get(chrom,sed-2,sed)
            strand = STED2STRAND.get(sted,'.')
            if strand != '.':
                css.append(strand)
            ureads,areads = _add_to_sj_arrays(sst,sed,dup,strand)
        pmid = rec[6]

    _append_sj(cse, css, csj, chrom, ureads, areads)

    _write_arrays()
    _write_sj(sjs)


def sj02wig(sjchr, chrom, chromsize, pathtmpl):
    a = {'+':N.zeros(chromsize, dtype=N.float64),
         '-':N.zeros(chromsize, dtype=N.float64),
         '.':N.zeros(chromsize, dtype=N.float64)}
    for st,ed,v,strand in sjchr[['st','ed','jcnt','strand']].values:
        a[strand][st-1:ed] += v
    for strand in a:
        path = pathtmpl.format(strand)
        cybw.array2wiggle_chr64(a[strand], chrom, path)
    

STRANDMAP0 = {'+':'.p','-':'.n','.':'.u'}

def sj02bw(sj0, pathpre, genome, np=12):
    chroms = UT.chroms(genome)
    chromdf = UT.chromdf(genome).sort_values('size',ascending=False)
    chroms = [x for x in chromdf['chr'] if x in chroms]
    chromdic = UT.df2dict(chromdf, 'chr', 'size')
    if 'jcnt' not in sj0:
        sj0['jcnt'] = sj0['ucnt']+sj0['mcnt']
    files = []
    args = []
    for c in chroms:
        f = '{0}.{1}.{{0}}.wig'.format(pathpre,c)
        args.append((sj0[sj0['chr']==c], c, chromdic[c], f))
        files.append(f)
    rslts = UT.process_mp(sj02wig, args, np=np, doreduce=False)
    rmfiles = []
    for strand in ['+','-','.']:
        s = STRANDMAP0[strand]
        wig = pathpre+'.sj{0}.wig'.format(s)
        bwpath = pathpre+'.sj{0}.bw'.format(s)
        with open(wig, 'w') as dst:
            for tmpl in files:
                f = tmpl.format(strand)
                with open(f,'r') as src:
                    shutil.copyfileobj(src, dst)
                rmfiles.append(f)
        rmfiles.append(wig)
        wig2bw(wig, UT.chromsizes(genome), bwpath)
    for f in rmfiles:
        os.unlink(f)
    os.unlink(wig)
    