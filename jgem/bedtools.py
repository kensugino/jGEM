"""Copyright (c) 2015-2016 Ken Sugino

.. module:: bedtools
    :synopsis: interface to bedtools

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
import os
import subprocess
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

#import pandas as PD
import numpy as N

from jgem import utils as UT
from jgem import gtfgffbed as GGB

def bedtoolintersect(aname, bname, cname, **kwargs):
    return _bedtoolscatcherror('intersect', aname, bname, cname, **kwargs)

def bedtoolmerge(aname, cname, **kwargs):
    return _runbedtools2('merge',aname, cname, **kwargs)

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
    with open(cname, "w") as outfile:
        ret = subprocess.call(cmd, stdout=outfile)
        if ret!=0:
            print(ret, cmd)
    return UT.compress(cname)

def _runbedtools3(which, aname, bname, cname, **kwargs):
    cmd = ['bedtools',which,'-a',aname,'-b',bname]
    for k,v in kwargs.items():
        if isinstance(v,bool):# in [True,False]: 2016-03-27 fix
            cmd += ['-'+k]
        else:
            cmd += ['-'+k, str(v)]
    with open(cname, "w") as outfile:
        ret = subprocess.call(cmd, stdout=outfile)
    if ret !=0:
        print('bederror', ret, cmd)
    return ret

def _bedtoolscatcherror(which, aname, bname, cname, **kwargs):
    if not os.path.exists(aname):
        raise ValueError('{0} does not exists'.format(aname))
    if not os.path.exists(bname):
        raise ValueError('{0} does not exists'.format(bname))
        
    if cname.endswith('.gz'):
        cname = cname[:-3]
    ret = _runbedtools3(which,aname,bname,cname,**kwargs)
    if ret !=0: # try uncompressed 
        aname = UT.uncompresscopy(aname)
        bname = UT.uncompresscopy(bname)
        ret = _runbedtools3(which,aname,bname,cname,**kwargs)
        if ret !=0:
            raise RuntimeError('bedtools error')
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



