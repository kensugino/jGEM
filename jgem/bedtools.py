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

import pandas as PD
import numpy as N

from jgem import utils as UT
from jgem import gtfgffbed as GGB
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
    awkscript = 'BEGIN{OFS="\t";c=1;}{if(d[$4]){$4=d[$4];}else{d[$4]=c;$4=c;c++;} n=split($11,a,","); n=split($12,b,","); for(i=1;i<=n;i++){st=$2+b[i]; print $1,st,st+a[i],$4,$5,$6,NR}}'
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
    awkscript = 'BEGIN{OFS="\t";c=1;}{if(a[$4]){$4=a[$4];}else{a[$4]=c;$4=c;c++;}; $7=NR; print $0;}'
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

def make_bw_from_bed(bedpath, chromsizes, bwpath):
    """ convert BED to BIGWIG, normalize average coverage to 1 """
    totbp,covbp = get_total_bp_bedfile(bedpath)
    scale = float(covbp)/totbp # 1/avgcov
    bed2bw(bedpath, chromsizes, bwpath, scale)

def make_bw_from_bam(bampath, chromsizes, bedpath, bwpath):
    """ convert BAM to BIGWIG, normalize average coverage to 1 """
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
    return UT.compress(cname)

def _bedtoolscatcherror2(which, aname, cname, **kwargs):
    if not os.path.exists(aname):
        raise ValueError('{0} does not exists'.format(aname))        
    if cname.endswith('.gz'):
        cname = cname[:-3]
    try:
        ret = _runbedtools2(which,aname,cname,**kwargs)
    except RuntimeError:
        LOG.warning('bedtool error: repeating on uncompressed a:{0},c:{1}'.format(aname,cname))
        aname2 = UT.uncompresscopy(aname)
        ret = _runbedtools2(which,aname2,cname,**kwargs)
        if aname2 != aname:
            os.unlink(aname2)
    return UT.compress(cname)

def calc_ovlratio(aname, bname, tname, nacol, nbcol, idcol=['chr','st','ed']):
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



