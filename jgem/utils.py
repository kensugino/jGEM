"""Copyright (c) 2015-2016 Ken Sugino

.. module:: utils
    :synopsis: a collection of utility functions

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
import subprocess
import gzip
import os
import errno
import math
import sys
import csv
import time
import collections
import functools
from functools import partial
import shutil
import uuid

import pandas as PD
import numpy as N
import scipy.optimize as SO

import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

#### DATE ##################################################################
def today():
    return time.strftime('%Y-%m-%d')

### MISC  ##################################################################
def isstring(obj):
    try:
        return isinstance(obj, basestring) #  python 2
    except NameError:
        return isinstance(obj, str) # python 3

### PANDAS UTIL     ########################################################

def calc_locus(df,chrfld='chr',stfld='st',edfld='ed'):
    """Create locus (chr*:start-end). 

    Args:
        df : Pandas.DataFrame 
        chrfld (str): field (column name) for chromosome 
        stfld (str): field (column name) for start position
        edfld (str): field (column name) for end position

    Returns:
        Pandas.Series
    """
    return df[chrfld]+':'+df[stfld].astype(str)+'-'+df[edfld].astype(str)

def calc_locus_strand(df,chrfld='chr', stfld='st',edfld='ed',strandfld='strand'):
    """Create locus with strand (chr*:start-end:+ or -). 

    Args:
        df : Pandas.DataFrame 
        chrfld (str): field (column name) for chromosome 
        stfld (str): field (column name) for start position
        edfld (str): field (column name) for end position
        strandfld (str); Field (column name) for strand

    Returns:
        Pandas.Series
    """
    return df[chrfld]+':'+df[stfld].astype(str)+'-'+df[edfld].astype(str)+':'+df[strandfld]


def flattendf(df, col):
    """Flattens a column whose elements are lists. One row gets expanded into multiple rows,
    each row containing a value from the list. All other columns in the expanded rows are the
    same as the original row. 

    Args:
        df: DataFrame
        col: target column which contains lists
    """
    i = list(df.columns).index(col)
    def _gen():
        for x in izipcols(df, df.columns):
            for y in x[i]:
                x1 = list(x)
                x1[i] = y
                yield x1
    df1 = PD.DataFrame([x for x in _gen()], columns = df.columns)
    return df1

def series2dict(series):
    "Make dict from Pandas.Series"
    return dict(zip(series.index.values, series.values))

def df2dict(df,f1,f2):
    """Make dict from Pandas.DataFrame

    Args:
        df: DataFrame
        f1: column for key
        f2: column for value

    Returns:
        dict
    """
    def _getval(fld):
        if fld=='index':
            return df.index.values
        return df[fld].values
    return dict(zip(_getval(f1), _getval(f2)))

def izipcols(df, cols, index=False):
    """Return an iterator to go through rows of Pandas.DataFrame
    (Much faster than DataFrame.rows() which involves instantiation of Series objects)

    Args:
        df: DataFrame
        cols: list of column names
        index: if True, includue index at the beginning (default False)
    """
    if index:
        l = [df.index]+[df[x] for x in cols]
    else:
        l = [df[x] for x in cols]
    #return izip(*l) # python2
    return zip(*l)


def remove_ercc(srcfile, dstfile):
    """Removes ERCC portion of the data

    Args:
        srcfile: source (GTF, GFF or BED)
        dstfile: destination

    """
    if srcfile[:-3]=='.gz':
        lines = gzip.open(srcfile).readlines()
    else:
        lines = open(srcfile).readlines()
    if dstfile[:-3]=='.gz':
        t = gzip.open(dstfile,'w')
    else:
        t = open(dstfile,'w')
    for line in lines:
        if not line.startswith('ERCC-'):
            t.write(line)
    t.close()



#### directory/file ########################################################
def makedirs(path):
    """Make all the directories along the path, calls os.makedirs but 
    intercept irrelevant exceptions
    """
    try:
        os.makedirs(path)
    except OSError as err:
        if err.errno != errno.EEXIST or not os.path.isdir(path):
            raise

def notstale(opath,dpath,override=False):
    """Check dpath exists and newer than opath.

    Args:
        opath (str or list): dependency file(s)
        dpath (str): cache file
        override (bool): if True always return False

    Returns:
        bool: whether dpath exists and newer than all of opath
    """
    if override:
        return False

    if not isinstance(opath, list):
        opath = [opath]
    for o in opath:
        if not os.path.exists(o):
            #raise RuntimeError('file {0} does not exist'.format(o))
            return False # stale
    notstale = os.path.exists(dpath)
    for o in opath:
        if not notstale:
            return notstale
        notstale = os.path.getmtime(dpath)>=os.path.getmtime(o)
    return notstale


#### read/write PANDAS #####################################################

def save_tsv_nidx_nhead(df, path, gzip=True, **kwargs):
    """Write tab separated file. No index, no header. 

    Args:
        df: Pandas DataFrame
        path (str): destination path 
        gzip (bool): whether to gzip compress or not (default True)
        kwargs: keyward arguments to be passed to Pandas.to_csv 

    Returns:
        saved file path
    """
    kwargs['index'] = False
    kwargs['header'] = False
    return save_tsv(df,path,gzip=gzip,**kwargs)
    
def save_tsv_nidx_whead(df, path, gzip=True, **kwargs):
    """Write tab separated file. No index, with header. 

    Args:
        df: Pandas DataFrame
        path (str): destination path 
        gzip (bool): whether to gzip compress or not (default True)
        kwargs: keyward arguments to be passed to Pandas.to_csv 

    Returns:
        saved file path
    """
    kwargs['index'] = False
    kwargs['header'] = True
    return save_tsv(df,path,gzip=gzip,**kwargs)

def save_tsv(df, path, gzip=True, **kwargs):
    """Write tab separated file. (Write index and header)

    Args:
        df: Pandas DataFrame
        path (str): destination path 
        gzip (bool): whether to gzip compress or not (default True)
        kwargs: keyward arguments to be passed to Pandas.to_csv 

    Returns:
        saved file path
    """
    if path[-3:]=='.gz':
        path = path[:-3]
    makedirs(os.path.dirname(path))
    df.to_csv(path, sep='\t', quoting=csv.QUOTE_NONE, **kwargs)
    if gzip:
        return compress(path)
    return path

def write_pandas(df, path, fm='h', **kwargs):
    """Write gzip compressed tab separated file.

    Args:
        df: Pandas DataFrame
        path (str): destination path 
        fm (str): 'i':write index, 'h':write header, 'ih':write both index and header
        kwargs: keyward arguments to be passed to Pandas.to_csv 

    Returns:
        saved file path
    """
    # fm = 'i','h','ih'
    if 'i' in fm:
        kwargs['index'] = True
    else:
        kwargs['index'] = False
    if 'h' in fm:
        kwargs['header'] = True
    else:
        kwargs['header'] = False
    return save_tsv(df, path, gzip=True, **kwargs)


def read_pandas(path,**kwargs):
    """Read tab separated file into Pandas DataFrame. 
    Automatically recognize gzip compressed or uncompressed files.

    Args:
        path (str): path to the file
        kwargs: keyward arguments to pass to Pandas.read_table
    """
    if not isstring(path):
        return path

    if path[-3:] == '.gz':
        if os.path.exists(path):
            return PD.read_table(path, compression='gzip', **kwargs)
        if os.path.exists(path[:-3]):
            return PD.read_table(path[:-3], **kwargs)
    else:
        if os.path.exists(path):
            return PD.read_table(path, **kwargs)
        if os.path.exists(path+'.gz'):
            return PD.read_table(path+'.gz', compression='gzip', **kwargs)
    raise RuntimeError('file {0} do not exists'.format(path))


#### GZIP ###############################################################
def compress(fname):
    "Uses gzip to compress file"
    if fname[-3:]=='.gz':
        return fname
    if os.path.exists(fname+'.gz'):
        os.unlink(fname+'.gz')
    ret = subprocess.call(['gzip', fname])
    if ret!=0:
        LOG.warning('Error compressing file {0}\n'.format(fname))
        return fname
    return fname+'.gz'
    
def uncompress(fname):
    "Unzip .gz file"
    if(fname[-3:]!='.gz'):
        return fname
    ret = subprocess.call(['gunzip', fname])
    if ret!=0:
        LOG.warning('Error uncompressing file {0}'.format(fname))
        return fname
    return fname[:-3]

def uncompresscopy(fname):
    "Uses unique file name to uncompress to avoid corruption in case of parallel computation."
    if(fname[-3:]!='.gz'):
        return fname
    # copy
    fname2 = fname[:-3]
    if os.path.exists(fname2):
        return fname2
    tmp = fname2+str(uuid.uuid4())
    with open(tmp,'w') as fobj:
        ret = subprocess.call(['gunzip','-c', fname], stdout=fobj)
        if ret!=0:
            LOG.warning('Error uncompressing file {0}'.format(fname))
            return fname
        os.rename(tmp, fname2)
    return fname2


#### EX,SJ attributes ###############################################

def exid(ex):
    "calculate exon id (chr:st-ed:strand)"
    if isinstance(ex,list): # ['chr','st','ed','name','sc1','strand']
            return '{0}:{1}-{2}:{5}'.format(*ex)
    return '{chr}:{st}-{ed}:{strand}'.format(**ex)
    
def sjid(sj):
    "calculate junction id (chr:st^ed:strand)"
    if isinstance(sj,list): # ['chr','st','ed','name','sc1','strand']
            return '{0}:{1}^{2}:{5}'.format(*sj)
    return '{chr}:{st}^{ed}:{strand}'.format(**sj)

## donor/acceptor, st/ed pos ids
def set_ids(df):
    "set '_id' field"
    df.sort_values(['chr','st','ed'], inplace=True)
    # 2016-02-27 df = df.sort_values(['chr','st','ed'])
    # is a separate copy of DataFrame and doesn't change original
    # fixed
    df['_id'] = N.arange(len(df))
    df.index = N.arange(len(df))
    
def set_ptyp(e):
    """[Deprecated]
    Set ptyp (exon type, internal 3prime, 5prime, single, etc.) according to exon name.

    """
    #LOG.warning('deprecated use set_info or set_category instead')

    il = e['name'].str.startswith('[')
    ir = e['name'].str.endswith(']')
    srp = e['strand']=='+'
    srn = e['strand']=='-'

    i3 = ((srp)&(il)&(~ir))|((srn)&(~il)&(ir))
    i5 = ((srp)&(~il)&(ir))|((srn)&(il)&(~ir))
    ie = (~srp)&(~srn)&(((~il)&(ir))|((il)&(~ir)))
    ii = il&ir
    ise = (~il)&(~ir)

    e.loc[ii,'ptyp'] = 'i' # internal
    e.loc[i3,'ptyp'] = '3' # 3prime
    e.loc[i5,'ptyp'] = '5' # 5prime
    e.loc[ie,'ptyp'] = 'e' # edge exon
    e.loc[ise,'ptyp'] = 's' # single exon

def _set_pos(bed, fld, stfld, edfld):
    "helper for set_ad_info"
    pos = bed[stfld].astype(int).astype(str)
    neg = bed[edfld].astype(int).astype(str)
    ppos = bed['chr']+':'+pos+':+'
    npos = bed['chr']+':'+neg+':-'
    opos = bed['chr']+':'+pos+':.'
    bed.loc[bed['strand']=='+', fld] = ppos
    bed.loc[bed['strand']=='-', fld] = npos
    bed.loc[bed['strand']=='.', fld] = opos

def set_ad_info(sj, me):
    """assign donor/acceptor id

    Args:
        sj: splice junction DataFrame
        me: (multi-)exon DataFrame

    SideEffects:
        add d_pos (donor position), a_pos (acceptor position), d_id, a_id, 
        d_degree, a_degree to sj and me

    """
    sj['st-1'] = sj['st']-1

    _set_pos(sj, 'd_pos', 'st-1', 'ed')
    _set_pos(sj, 'a_pos', 'ed', 'st-1')
    _set_pos(me, 'd_pos', 'ed', 'st')
    _set_pos(me, 'a_pos', 'st', 'ed')

    nullid = 0
    # donor id (did)
    dids = sj['d_pos'].unique()
    #dids = dids[dids!='']
    dpos2did = dict([(x,i+1) for i,x in enumerate(dids)])
    sj['d_id'] = [dpos2did[x] for x in sj['d_pos']]
    #me['d_id'] = [dpos2did.get(x,-1) for x in me['d_pos']]
    me['d_id'] = [dpos2did.get(x,nullid) for x in me['d_pos']]
    # acceptor id (aid)
    aids = sj['a_pos'].unique()
    #aids = aids[aids!='']
    apos2aid = dict([(x,i+1) for i,x in enumerate(aids)])
    sj['a_id'] = [apos2aid[x] for x in sj['a_pos']]
    #me['a_id'] = [apos2aid.get(x,-1) for x in me['a_pos']]
    me['a_id'] = [apos2aid.get(x,nullid) for x in me['a_pos']]
    
    # donor degree, acceptor degree (how many edges)
    ddeg = sj.groupby('d_pos').size()
    adeg = sj.groupby('a_pos').size()
    dpos2deg = dict(zip(ddeg.index.values, ddeg.values))
    apos2deg = dict(zip(adeg.index.values, adeg.values))
    sj['d_degree'] = [dpos2deg[x] for x in sj['d_pos']]
    sj['a_degree'] = [apos2deg[x] for x in sj['a_pos']]
    me['d_degree'] = [dpos2deg.get(x,0) for x in me['d_pos']]
    me['a_degree'] = [apos2deg.get(x,0) for x in me['a_pos']]    
    
def set_pos_info(sj, me):
    # assign donor/acceptor id
    sj['st-1'] = sj['st']-1
    
    def _set_sted_pos(bed, stfld, edfld):
        bed['st_pos'] = bed['chr']+':'+bed[stfld].astype(int).astype(str)
        bed['ed_pos'] = bed['chr']+':'+bed[edfld].astype(int).astype(str)

    _set_sted_pos(sj, 'st-1', 'ed')
    _set_sted_pos(me, 'ed', 'st' )

    # st_pos id (stid)
    dids = sj['st_pos'].unique()
    dpos2did = dict([(x,i+1) for i,x in enumerate(dids)])
    sj['st_id'] = [dpos2did[x] for x in sj['st_pos']]
    me['st_id'] = [dpos2did.get(x,-1) for x in me['st_pos']]
    # acceptor id (aid)
    aids = sj['ed_pos'].unique()
    apos2aid = dict([(x,i+1) for i,x in enumerate(aids)])
    sj['ed_id'] = [apos2aid[x] for x in sj['ed_pos']]
    me['ed_id'] = [apos2aid.get(x,-1) for x in me['ed_pos']]

def find_nullidx(ex):
    amin = ex['a_id'].min()
    if amin<-1: # merged
        nullidx = 0
    elif amin==-1: # old
        nullidx = -1
    else:
        nullidx = 0
    return nullidx

def set_exon_category(sj,ex):
    # 5', 3', i, s
    if ('a_pos' not in sj.columns) or ('a_pos' not in ex.columns):
        set_ad_info(sj, ex) 
    nullidx = find_nullidx(ex)
    aidx = ex['a_id']!=nullidx
    didx = ex['d_id']!=nullidx
    idxs = (~aidx)&(~didx)
    idx3 = aidx&(~didx)
    idx5 = (~aidx)&didx
    idxi = aidx&didx
    ex.loc[idxs, 'cat'] = 's'
    ex.loc[idx5, 'cat'] = '5'
    ex.loc[idx3, 'cat'] = '3'
    ex.loc[idxi, 'cat'] = 'i'
    
def set_info(sj,me):
    if '_id' not in sj.columns:
        set_ids(sj)
    if '_id' not in me.columns:
        set_ids(me)
    if ('a_id' not in sj.columns) or ('d_id' not in sj.columns) or \
       ('a_id' not in me.columns) or ('d_id' not in me.columns):
        set_ad_info(sj, me) 
    

## ME/SE
def mese(ex):
    #LOG.warning('deprecated, use me_se instead')
    idx = ex['cat']=='s'
    se = ex[idx]
    me = ex[~idx]
    return me, se
    
def me_se(ex):
    "Separate exons into multi-exons and single-exons"
    nullidx = find_nullidx(ex)
    seidx = (ex['a_id']==nullidx)&(ex['d_id']==nullidx)
    se = ex[seidx]
    me = ex[~seidx]
    return me, se



#### BINNED AVG ########################################################
def calc_binned(xs, ys, num=500, yth=0, returnsorted=False, returnminmax=False, method='mean'):
    idx = range(len(xs))
    tmp = sorted(zip(xs, idx), reverse=True)
    n = int(math.ceil(len(tmp)/float(num))) # number of bins
    bidx = N.array([x*num for x in range(n)]+[len(tmp)]) # start (& end) indices for bins
    h0 = bidx[1:]-bidx[:-1] # number of items in each bin
    sted = list(zip(bidx[:-1],bidx[1:])) # starts and ends
    func = getattr(N, method)
    avgx = N.array([func([x[0] for x in tmp[st:ed]]) for st,ed in sted])
    avgy = N.array([func([ys[x[1]] for x in tmp[st:ed] if x[1]>=yth]) for st,ed in sted])
    idxst = N.array([st for st,ed in sted])
    if returnsorted:
        sx = [x[0] for x in tmp]
        sy = [ys[x[1]] for x in tmp]
        return avgx,avgy,idxst[::-1],sx[::-1],sy[::-1]
    if returnminmax:
        minx = N.array([N.min([x[0] for x in tmp[st:ed]]) for st,ed in sted])
        maxx = N.array([N.max([x[0] for x in tmp[st:ed]]) for st,ed in sted])
        cnt = N.array([len(tmp[st:ed]) for st,ed in sted])
        return avgx,avgy,minx,maxx, cnt
    return avgx, avgy

def calc_binned_percentile(xs, ys, percent=99, minbinw=0.25, minnum=500):
    idx = range(len(xs))
    tmp = sorted(zip(xs, idx), reverse=True)
    def _genbin():
        ci = 0
        yield 0
        for i in range(1,len(tmp)):
            if (i-ci > minnum) & (tmp[ci][0]-tmp[i][0] > minbinw):
                yield i
                ci = i
    bidx = [x for x in _genbin()]
    sted = list(zip(bidx[:-1],bidx[1:])) # starts and ends
    px = N.array([N.mean([x[0] for x in tmp[st:ed]]) for st,ed in sted])
    py = N.array([N.percentile([ys[x[1]] for x in tmp[st:ed]], percent) for st,ed in sted])
    return px,py

#### Sigmoid fitting #########################################################

def sigmoid(x,x0,k):
    return 1./(1.+N.exp(-k*(x-x0)))

def invsig(y,x0,k):
    return x0 - N.log((1.-y)/y)/k

def fit_sigmoid(x,y,xlim=(0,7),yth=0.99):
    popt,pcov= SO.curve_fit(sigmoid, x, y)
    x2 = N.linspace(*xlim)
    y2 = sigmoid(x2, *popt)
    xth = invsig(yth, *popt)
    return x2,y2,xth


#### Genome info #############################################################

def chromsizes(genome):
    "returns the path to chrom.sizes file"
    return os.path.join(os.path.dirname(__file__), 'data', '{0}.chrom.sizes'.format(genome))

def chromdf(genome):
    "returns chromosome size dataframe"
    return read_pandas(chromsizes(genome), names=['chr','size'])

def chroms(genome):
    "returns chromosome names"
    df = read_pandas(chromsizes(genome), names=['chr','size'])
    return list(df['chr'])


### Atomic Intervals ##################################################

def chopintervals(exons, fname=None, sort=True, idcol='_id'):
    """Separate into intervals where overlaps are constant over each interval.

    Args:
        exons: DataFrame containing exons
        fname (str): (optional) save results into this file
        sort (bool): whether to sort the DataFrame
        idcol (str): column name to use for id

    Returns:
        A DataFrame containing chopped intervals
    """
    # Assumes st<ed (not even st==ed)
    if sort:
        exons = exons.sort_values(['chr','st','ed'])
    if idcol=='_id' and '_id' not in exons.columns:
        exons['_id'] = N.arange(len(exons))
        exons.index = N.arange(len(exons))
    # process by chrom
    def _todf2(df,tgt,kval,idfld):
        tmp = PD.DataFrame(df[tgt].values, columns=['pos'])
        tmp['kind'] = kval
        tmp['id'] = df[idfld].astype(str).values
        return tmp
    def _gen():
        for chrom in exons['chr'].unique():
            ec = exons[exons['chr']==chrom]
            # generate array with positions (st,ed) and kind (st,ed) and id (_id)
            posarr = PD.concat([_todf2(ec,'st',1,idcol),
                                _todf2(ec,'ed',-1,idcol)],
                               ignore_index=True).sort_values('pos')
            posarr['cumsum'] = N.cumsum(posarr['kind'])
            p1,cs1,id1 = posarr.iloc[0][['pos','cumsum','id']]
            idset = set([id1])
            for p2,cs2,id2 in posarr[['pos','cumsum','id']].values[1:]:
                name = ','.join(idset)
                if cs2>cs1: # start
                    idset.add(id2)
                else:
                    idset.remove(id2)
                if (cs1>0) and (p2>p1):
                    yield (chrom,p1,p2,name)
                p1,cs1,id1 = p2,cs2,id2
    ci = PD.DataFrame([x for x in _gen()], columns=['chr','st','ed','name'])
    ci['id'] = N.arange(len(ci))
    if fname:
        save_tsv_nidx_nhead(ci, fname)
    return ci
