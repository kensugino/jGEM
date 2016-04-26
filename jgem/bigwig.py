"""

.. module:: bigwig
    :synopsis: BIGWIG file related stuffs

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""

# system
import os
import subprocess
import time
import multiprocessing
from operator import iadd
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)
import shutil

# 3rd party
import pandas as PD
import numpy as N

# for reading BIGWIG file
#ngslib is faster than bx but fails when parallerized
#import bx
#from bx.bbi.bigwig_file import BigWigFile
#from ngslib import wWigIO
#bx-python is not compatible with Python3 => use modified version 

from jgem import utils as UT
from jgem import bedtools as BT
import jgem.cy.bw  as cybw #import array2wiggle_chr # Cython version
from jgem.cy.bw import array2wiggle_chr
from jgem.bxbbi.bigwig_file import BigWigFile

MAXSIZE = int(300e6)  # 300Mbp bigger than chr1,chrX

# BAM to genomecov BIGWIG ##################################################
def cnt_bam(fpath):
    """
    Uses samtools to count reads in BAM file. 

    Args:
        fpath (str): path to BAM file

    Returns:
        int. The number of aligned reads.

    Raises:
        RuntimeError
    """
    cache = fpath+'.flagstat'
    bai = fpath+'.bai'
    if not os.path.exists(fpath):
        raise RuntimeError('file {0} does not exist'.format(fpath))
    if not os.path.exists(bai):
        cmd = ['samtools', 'index', fpath]
        subprocess.call(cmd)
    if os.path.exists(cache) and os.path.getmtime(cache)>os.path.getmtime(fpath):
        out = open(cache,'r').read()
    else:
        cmd = ['samtools', 'flagstat', fpath]
        out = subprocess.check_output(cmd)
        open(cache,'w').write(out)
    firstline = out.split('\n')[0].split()
    return int(firstline[0])+int(firstline[2])


def wig2bw(wigpath, chromsizes, bwpath):
    pass

def bam2bw(fpath, chromsizes, bpath, aligned=None):
    """
    Generate normalized coverage from BAM

    Args:
        fpath (str): path to BAM
        chromsizes (str): path to chromsizes file 
        bpath (str): path to BIGWIG
        aligned (int): number of aligned reads, if None uses samtools to find it from BAM

    Requires Bedtools (genomeCoverageBed) and Kent Tool (wigToBigWig)

    """
    # countreads
    if aligned is None:
        aligned =   cnt_bam(fpath)
    scale = 1000000./float(aligned)
    # convert_to_wig
    tpath = bpath +'.wig'
    tfobj = open(tpath,'w')
    cmd1 = ['genomeCoverageBed', '-split', '-bg', '-ibam', fpath, '-g', chromsizes, '-scale', str(scale)]
    p1 = subprocess.Popen(cmd1, stdout=tfobj)
    p1.wait()
    tfobj.close()

    # convet_wig_to_bigwig
    cmd2 = ['wigToBigWig', tpath, chromsizes, bpath]
    p2 = subprocess.call(cmd2)

    # remove_temporary_file
    os.remove(tpath)


# bw2bed based on #########################################################
# https://github.com/CGATOxford/cgat/blob/master/scripts/wig2bed.py
# [TODO] report bug to AndreasHeger line 80 "elif => if"

def block_iter(infile, chrom, chunk=int(10e6)):
    "BigWig file iterator"
    with open(infile, mode='rb') as fobj:
        bw = BigWigFile(fobj)
        for x in range(0,MAXSIZE,chunk): # 10Mbp chunk
            iterator = bw.get(chrom, x, x+chunk)
            if iterator is None:
                raise StopIteration
            for v in iterator:
                yield v

# def block_iter_ngs(infile, chrom, chunk=int(10e6)):
#     # ngslib is faster than bx but fails when parallerized
#     wWigIO.open(infile)
#     for x in range(0,MAXSIZE,chunk): # 10Mbp chunk
#         iterator = wWigIO.getIntervals(infile, chrom, x, x+chunk)
#         if not iterator:
#             raise StopIteration
#         for v in iterator:
#             yield v
#     wWigIO.close(infile)

def apply_threshold(infile, threshold, chroms):
    """Only returns intervals exceeding the threshold

    Args:
        infile: path to bigwig
        threshold: positive float
        chroms: list of chromosome names

    Yields:
        intervals (chr,st,ed)
    """
    for chrom in chroms:
        last_start, last_end = -1, 0
        #for start, end, value in block_iter_ngs(infile, chrom):
        for start, end, value in block_iter(infile, chrom): 
            # start=>end has value
            d = start - last_end 
            # d: distance from last_end
            # d> 0 => last_end==>start value is 0 assume threshold>=0
            if (d > 0 or value <= threshold):
                if last_start >= 0: # if there's un-yielded interval then yield
                    yield chrom, last_start, last_end
                last_start = -1 # reset (no hit)
            #elif last_start < 0 and value > threshold: 
            # this part is a bug from the original code
            # original code will skip the next interval
            if last_start < 0 and value > threshold:
                # last_start <0 => no current interval
                # and value above threshod ==> start new interval
                last_start = start
            last_end = end
        if last_start >= 0:
            yield chrom, last_start, end

def bw2bed(bwfile, bedfile, chroms, th):
    """Transform BigWig genomeCov to binary BED by thresholding. 
    Makes result file (bwfile[:-3]+'.binary%g.bed'.format(th))

    Args:
        bwfile: path to BigWig file
        chroms: list of chromosome names
        th: coverage threshold

    Returns:
        path to generated BED file
    """
    #bedfile = '{0}.binary{1:g}.bed'.format(bwfile[:-3], th)
    if UT.notstale(bwfile, bedfile+'.gz'):
        return bedfile+'.gz'
    # make sure bwfile exists
    if not ( os.path.exists(bwfile) ):
        raise RuntimeError('BigWig file {0} does not exist.'.format(bwfile))
    processor = apply_threshold(bwfile,th, chroms)
    out = open(bedfile,'w')
    out.write(''.join(['%s\t%i\t%i\n' % x for x in processor]))
    out.write('\n')
    out.close()
    return UT.compress(bedfile)

### Merge BigWigs ##########################################################

def get_bigwig_as_array(bwfile, chrom, st, ed):
    """Get BIGWIG coverage values as array of size (ed-st). Array start corresponds to st.
    0-based

    Args:
        bwfile: path to BIGWIG
        chrom (str): chromosome name
        st (int): start position
        ed (int): end position

    Returns:
        Numpy array of size (ed-st)
    """
    # with open(bwfile, mode='rb') as fobj:
    #     bw = BigWigFile(fobj)
    #     it = bw.get(chrom,st,ed)
    #     a = N.zeros(ed-st)
    #     for s,e,v in it:
    #         a[s-st:e-st] += v
    # return a
    with open(bwfile, mode='rb') as fobj:
        bw = BigWigFile(fobj)
        a = bw.get_as_array(chrom,st,ed)
        a[N.isnan(a)]=0.
    return a

def merge_bigwigs_chr(bwfiles, chrom, chromsize, dstpath, scale):
    # merge4-allsample.bw chr1 89026991 intervals ~50%
    # better to just use dense array than sparse array
    a = N.zeros(chromsize)
    for fpath in bwfiles:
        with open(fpath,mode='rb') as fobj:
            bw = BigWigFile(fobj)
            it = bw.get(chrom, 0, chromsize)
            if it is not None:
                for s,e,v in it:
                    a[s:e] += v
    #a = a/float(len(bwfiles))
    if scale is not None:
        a = a*scale
    a = N.array(a, dtype=N.float32)
    cybw.array2wiggle_chr(a, chrom, dstpath)
    return (chrom, dstpath)

def merge_bigwigs_mp(bwfiles, genome, dstpath, scale=None, np=7):
    chroms = UT.chroms(genome)
    chromfile = UT.chromsizes(genome)
    chromsizes = UT.df2dict(UT.chromdf(genome), 'chr', 'size')
    args = [(bwfiles, c, chromsizes[c], dstpath+'.{0}.wig'.format(c), scale) for c in chroms]

    rslts = UT.process_mp(merge_bigwigs_chr, args, np, doreduce=False)

    dic = dict(rslts)
    LOG.debug('concatenating chromosomes...')
    wigpath = dstpath+'.wig'
    with open(wigpath, 'w') as dst:
        for c in chroms:
            with open(dic[c],'r') as src:
                shutil.copyfileobj(src, dst)

    LOG.debug('converting wiggle to bigwig')
    BT.wig2bw(wigpath, chromfile, dstpath)

    # clean up 
    for c in chroms:
        f = dstpath+'.{0}.wig'.format(c)
        if os.path.exists(f):
            os.unlink(f)
    if os.path.exists(wigpath):
        os.unlink(wigpath)
    
# def array2wiggle_chr(a, chrom, dstpath):
    # possibly Cythonify
    # def _gen():
    #     i = 0
    #     # skip initial 0
    #     while(a[i]==0):
    #         i+=1
    #     st = i
    #     c = a[st] # initial non-zero
    #     i+=1
    #     while(i<len(a)):
    #         # advance until change
    #         while(i<len(a) and a[i]==c):
    #             i+=1
    #         if c!=0:
    #             # 0-based => 1-based
    #             yield '{0}\t{1}\t{2}\t{3}\n'.format(chrom,st,i,c)
    #         if i<len(a):
    #             st = i
    #             c = a[st]
    # with open(dstpath,'w') as fobj:
    #     if N.sum(a)>0: # some elements are not zero
    #         cnt = 0
    #         txt = []
    #         for line in _gen():
    #             txt.append(line)
    #             cnt += 1
    #             if cnt == 100000:
    #                 fobj.write(''.join(txt))
    #                 cnt = 0
    #                 txt = []
    #         fobj.write(txt)
    # return dstpath

