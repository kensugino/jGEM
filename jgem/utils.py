"""

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
from functools import partial, reduce
from operator import iadd
from itertools import repeat, groupby
try:
    from itertools import izip_longest as zip_longest
except:
    from itertools import zip_longest
import json
import shutil
import uuid
import multiprocessing

import pandas as PD
import numpy as N
import scipy.optimize as SO

from jgem import taskqueue as TQ

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

def make_dict(df,f1,f2):
    """Make dict but maps to sets
    """
    dic = {}
    for k,v in df[[f1,f2]].values:
        dic.setdefault(k,set()).add(v)
    dic = {k:list(dic[k]) for k in dic}
    return dic

def subdict(dic, keys):
    return {k:dic[k] for k in keys}

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
    UT.makedirs(os.path.dirname(dstfile))
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
    if path=='':
        path='./'
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

def transpose_csv(src, dst, linesep=b'\n', colsep=b'\t'):
    if src[-3:]=='.gz':
        sp = gzip.open(src,'r') # byte
    else:
        sp = open(src,'rb')
    rows = [x.strip().split(colsep) for x in sp]
    sp.close()
    
    n = len(rows[0])
    cols = [colsep.join([x[i] for x in rows]) for i in range(n)]
    if dst[-3:]=='.gz':
        dp = gzip.open(dst,'w')
    else:
        dp = open(dst,'wb')
    dp.write(linesep.join(cols))
    dp.close()

def save_json(obj, fname):
    """ save obj (usually dict) in JSON format """
    makedirs(os.path.dirname(fname))
    with open(fname,'w') as fp:
        json.dump(obj, fp)

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
    df = check_int_nan(df, ['chr','st','ed'], ['st','ed'])
    df.to_csv(path, sep='\t', quoting=csv.QUOTE_NONE, **kwargs)
    if gzip:
        return compress(path)
    return path

def check_integer(df, cols):
    if any([df.dtypes.get(x,int)!=int for x in cols]):
        LOG.warning('col not integer: copy and converting')
        df = df.copy()
        for x in cols:
            df[x] = df[x].astype(int)
    return df

def check_int_nan(d, when=['chr','st','ed'], icols=['st','ed']):
    if not all([x in d.columns for x in when]):
        return d # only check when there are all of these
    idx = d[when[0]].isnull()
    for c in when[1:]:
        idx = idx |d[c].isnull()
    nannum = N.sum(idx)
    notint = any([d.dtypes[c] != int for c in icols])
    if nannum>0 or notint:
        d = d[~idx].copy()
        if nannum>0:
            LOG.warning('{0} NaN in chr/st/ed, discarding'.format(N.sum(idx)))
        if notint:
            LOG.warning('st,ed not integer')
        d['st'] = d['st'].astype(int)
        d['ed'] = d['ed'].astype(int)
    return d


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
    if path[-3:]=='.gz':
        gzip=True
    else:
        gzip=False
    return save_tsv(df, path, gzip=gzip, **kwargs)


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

def make_empty_df(colnames):
    """ make an empty Pandas dataframe with colnames """
    return PD.DataFrame(N.zeros((0,len(colnames))),columns=colnames) 


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
    ext = '.'+fname[:-3].split('.')[-1]
    tmp = fname2+'.'+str(uuid.uuid4())+ext
    with open(tmp,'wb') as fobj:
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
    bed.loc[bed['strand'].isin(['+','.+']), fld] = ppos
    bed.loc[bed['strand'].isin(['-','.-']), fld] = npos
    bed.loc[bed['strand'].isin(['.',]), fld] = opos

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
    if ('a_id' not in sj.columns) or ('a_id' not in ex.columns):
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

def calc_tlen(ex, ci, gidx='_gidx'):
    """ calculate gene transcript length, here just union of exons """
    ci['len'] = ci['ed']-ci['st']
    ci['eid'] = ci['name'].astype(str).apply(lambda x: [int(y) for y in x.split(',')])
    if ('id' not in ci.columns):
        ci['id'] = ci['sc1']
    ccf = flattendf(ci[['eid','chr','st','ed','len','id','name']], 'eid')
    e2g = df2dict(ex, '_id',gidx)
    ccf[gidx] = [e2g[x] for x in ccf['eid']]
    # ignore duplicated ci, just take unique set within the gene
    ccf1 = ccf.groupby([gidx,'id']).first().reset_index() 
    gbed = ccf1.groupby(gidx)['len'].sum().to_frame('tlen')
    gbed['ltlen'] = N.log10(gbed['tlen'])
    return gbed # column: tlen, ltlen, index:_gidx

def set_glen_tlen(ex,ci,gidx='_gidx'):
    if 'len' not in ex.columns:
        ex['len'] = ex['ed']-ex['st']        
    tlen = calc_tlen(ex, ci, gidx)
    g2tlen = df2dict(tlen, 'index', 'tlen')
    ex['tlen'] = [g2tlen[x] for x in ex[gidx]]
    gr = ex.groupby(gidx)
    g2gst = series2dict(gr['st'].min())
    g2ged = series2dict(gr['ed'].max())
    ex['gst'] = [g2gst[x] for x in ex[gidx]]
    ex['ged'] = [g2ged[x] for x in ex[gidx]]
    ex['glen'] = ex['ged']-ex['gst']
    # glen = gr['ed'].max() - gr['st'].min()
    # g2glen = series2dict(glen)
    # ex['glen'] = [g2glen[x] for x in ex['_gidx']]


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
    try:
        popt,pcov= SO.curve_fit(sigmoid, x, y)
    except RuntimeError:
        x2 = N.linspace(*xlim)
        y2 = N.array([N.nan]*len(x2))
        xth = N.nan
        return x2,y2,xth
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
    path = os.path.join(os.path.dirname(__file__), 'data', '{0}.chromosomes'.format(genome))
    df = read_pandas(path,names=['chr'])
    return list(df['chr'])


### Atomic Intervals ##################################################

def chopintervals(exons, fname=None, sort=True, idcol='_id', poscol='_pid'):
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
    if idcol not in exons.columns:
        exons[idcol] = N.arange(len(exons))
        exons.index = N.arange(len(exons))
    if poscol not in exons.columns:
        exg = exons.groupby(['chr','st','ed'])[[idcol]].first()
        exg['_pid'] = N.arange(len(exg)) # position id
        cse2pid = dict(zip(exg.index,exg['_pid']))
        exons['_pid'] = [cse2pid[tuple(x)] for x in exons[['chr','st','ed']].values]

    exu = exons.groupby('_pid').first().reset_index()
    # process by chrom
    def _todf2(df,tgt,kval,idfld):
        tmp = PD.DataFrame(df[tgt].values, columns=['pos'])
        tmp['kind'] = kval
        tmp['id'] = df[idfld].astype(str).values
        return tmp
    def _gen():
        for chrom in exu['chr'].unique():
            ec = exu[exu['chr']==chrom]
            # generate array with positions (st,ed) and kind (st,ed) and id (_id)
            posarr = PD.concat([_todf2(ec,'st',1,poscol),
                                _todf2(ec,'ed',-1,poscol)],
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
        write_pandas(ci, fname, '')
    return ci


def read_ci(cipath):
    """ ci has txt ext but don't have header """
    # should have been assigned bed ext 
    header = ['chr','st','ed','name','id']
    return read_pandas(cipath, names=header)


def union_contiguous(beddf, returndf=True, pos_cols=['chr','st','ed']):
    """ Union contiguous records into one. 
    Uses chr,st,ed (not strand) columns.
    For other fields, values from the first record is used. 
    """
    # pos_cols = ['chr','st','ed']
    beddf = beddf.sort_values(pos_cols)
    cols = list(beddf.columns)
    pos_idx = [cols.index(x) for x in pos_cols]
    def _gen():
        rec0 = beddf.iloc[0].values
        chrpos,stpos,edpos = pos_idx
        for rec1 in beddf[cols].values[1:]:
            if rec1[chrpos] != rec0[chrpos]: # chrom1 != chrom0 switch chrom
                yield rec0
                rec0 = rec1
            elif rec0[edpos] < rec1[stpos]: # ed1<st0 new interval
                yield rec0
                rec0 = rec1
            else: # overlapping/contiguous
                #rec0[edpos] = rec1[edpos] # update the end <== bug if next record
                # is include in the previous, end will be wrong
                rec0[edpos] = max(rec0[edpos], rec1[edpos])
        yield rec0
    recs = [x for x in _gen()]
    if not returndf:
        return recs
    df = PD.DataFrame(recs, columns=cols)
    return df

def union_contiguous_intervals(arr):
    """ Union contiguous intervals.
    arr = list of [(st,ed),...]
    """
    arr = sorted(arr)
    def _gen():
        st0,ed0 = arr[0]
        for st1,ed1 in arr[1:]:
            if ed0 < st1: # ed1<st0 new interval
                yield [st0,ed0]
                st0,ed0 = st1,ed1
            else: # overlapping/contiguous
                ed0 = max(ed0, ed1)
        yield [st0,ed0]
    return [x for x in _gen()]

def make_unionex(ex, gidx='_gidx'):
    """Makes gene bed df where overlapping exons belonging to a genes
    are concatenated.

    (Takes ~2min for gencodeVM4)
    """
    def _gen():
        for k, v in ex.groupby(gidx):
            for x in union_contiguous(v,returndf=False):
                yield x
    recs = [x for x in _gen()]
    return PD.DataFrame(recs, columns=ex.columns)


def name2gidx(s):
    if s[:2]=='JN':
        return -int(s[3:])
    if s[:2]=='JP':
        return int(s[3:])
    if s[:2]=='JS':
        return int(s[2:])

def tname2gidx(tname):
    gname = tname.split('.')[0]
    return name2gidx(gname)

#### multiprocessing ##################################################

def mp_worker(args):
    func, arg = args
    return func(*arg)

def process_mp(func, args, np, doreduce=True):
    rslts = []
    if np==1:
        for i, arg in enumerate(args):
            rslts.append(func(*arg))
            LOG.debug(' processing: {0}/{1}...'.format(i+1,len(args)))
    else:
        try:
            p = multiprocessing.Pool(np)
            a = zip(repeat(func), args)
            rslts = p.map(mp_worker, a)
        finally:
            LOG.debug('closing pool')
            # p.join()
            p.close()
    if doreduce:
        rslts = reduce(iadd, rslts, [])
    return rslts    


def process_mp2(func, args, np, doreduce=True):
    rslts = []
    if np==1:
        for i, arg in enumerate(args):
            rslts.append(func(*arg))
            LOG.debug(' processing: {0}/{1}...'.format(i+1,len(args)))
    else:
        status = {}
        server = TQ.Server(np=np)
        with server:
            for i,a in enumerate(args):
                tname = 'func.{0}'.format(i)
                task = TQ.Task(tname, func, a)
                server.add_task(task)
            while server.check_error():
                try:
                    name, rslt = server.get_result(timeout=5)
                except TQ.Empty:
                    name, rslt = None, None
                if name is not None:
                    subid = name.split('.')[-1]
                    status[subid] = rslt
                    if len(status)==len(args):
                        break
        rslts = status.values()
    if doreduce:
        rslts = reduce(iadd, rslts, [])
    return rslts    




def grouper(iterable, n, fill=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper('ABCDEFG', 3, 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return zip_longest(fillvalue=fillvalue, *args)



