"""Copyright (c) 2015-2016 Ken Sugino

.. module:: calccov
    :synopsis: calculate coverages, gene length using nnls or 

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
# system imports
import multiprocessing
import os
import time
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)
from itertools import repeat
from functools import reduce

# 3rd party imports
import pandas as PD
import numpy as N
from scipy.optimize import nnls

# library imports
from jgem import utils as UT
from jgem import bigwig as BW
from jgem import gtfgffbed as GGB


### All in one func  ##################################################

def calc_cov_ovl_mp(srcname, bwname, dstname, np=1, covciname=None, 
    ciname=None, colname='cov', override=False):
    """Calculate coverage (from BigWig) over intervals (from srcname). 
    A column (default 'cov') which contains coverages is added to source dataframe 
    and the source is overwritten. 

    Args:
        srcname: path to exons tsv
        bwname: path to bigwig
        dstname: path for result
        np: number of processors
        covciname: path to covci (coverage for chopped interval dataframe)
        ciname: path to ci (chopped interval dataframe)
        colname: name for column which contain calculated coverages

    Returns:
        source dataframe with column (cov) added

    SideEffects:
        source tsv is overwritten with new column added

    """
    if UT.isstring(srcname):
        exons = UT.read_pandas(srcname)
    else:
        exons = srcname
    # cache
    if covciname is None:
        assert(UT.isstring(srcname))
        covciname = srcname[:-7]+'.covci.txt.gz'
    if ciname is None:
        assert(UT.isstring(srcname))
        ciname = srcname[:-7]+'.ci.txt.gz'

    if override or (not os.path.exists(covciname)):
        LOG.debug('calculating covci...')
        _sttime = time.time()
        if override or not (os.path.exists(ciname)):
            ci = UT.chopintervals(exons, ciname)
        else:
            ci = UT.read_pandas(ciname, names=['chr','st','ed','name','id'])
        covci = calc_cov_mp(ci,bwname,covciname,np)
        LOG.debug(' time: {0:.3f}s'.format(time.time()-_sttime))
    else:
        LOG.debug('loading cached covci...')
        covci = UT.read_pandas(covciname)

    # covci: chopped interval's cov => reverse
    # ci => exon id ====> revers exon => ci indices
    # exon cov = sum(cicov*cilen)/totlen
    LOG.debug('calculating exon cov...')
    if 'id' not in covci.columns:
        covci['id'] = covci['sc1']
        
    _sttime = time.time()
    e2c = {}
    for i,name in covci[['id','name']].values:
        for eid in name.split(','):
            e2c.setdefault(int(eid),[]).append(i)
    covci['len'] = covci['ed']-covci['st']
    covci['val'] = covci['cov']*covci['len']            
    def _gen():
        for eid in exons['_id']:
            for cid in e2c[eid]:
                yield (cid,eid)
    tmp = PD.DataFrame(list(set([x for x in _gen()])),columns=['cid','eid'])
    c2len = dict(covci[['id','len']].values)
    c2val = dict(covci[['id','val']].values)
    tmp['val'] = [c2val[x] for x in tmp['cid']]
    tmp['len'] = [c2len[x] for x in tmp['cid']]
    tmpg = tmp.groupby('eid')[['val','len']].sum().reset_index()
    tmpg['cov'] = tmpg['val']/tmpg['len']
    e2cov = dict(tmpg[['eid','cov']].values)
    exons[colname] = [e2cov[x] for x in exons['_id']]
    
    UT.save_tsv_nidx_whead(exons, dstname)
    return exons

### ecov calc with nnls ################################################
def mp_worker(args):
    func, arg = args
    return func(*arg)

def calc_ecov_chrom(covci,blocksize):
    # covci: required cols
    # name (eid str concat ',')
    # name1 ([eids,...])
    # id: cid
    # cov: coverage for cid
    covci = covci.sort_values('id')
    if 'name1' not in covci.columns:
        covci['name1'] = covci['name'].astype(str).apply(lambda x: [int(y) for y in x.split(',')])
    if 'eidmax' not in covci.columns:
        covci['eidmax'] = covci['name1'].apply(lambda x:max(x))
    if 'eidmin' not in covci.columns:
        covci['eidmin'] = covci['name1'].apply(lambda x:min(x))
    e2c = {} # dict
    # find chunk boundary with approximate blocksize
    def _find_chunk(st, bs):
        # first take stcid => stcid+blocksize
        # reduce to self contained chunk => find isolated cid<=>eid (len(name1)==1) instance        
        # if final size is too small increase blocksize and try again
        if st+bs+1>=len(covci):
            return len(covci)
        for x in range(st+bs,int(st+bs/2),-1):
            if len(covci['name1'].iloc[x])==1:# found isolated exon
                return x+1
        # not found => increase blocksize
        if bs>500: # blocksize too big => NNLS take too much time
            LOG.warning('blocksize exceeded 500')
            return st+bs
        return _find_chunk(st, 2*bs)
    # calc ecov for a chunk
    def _calc_ecov(st,ed):
        # construct matrix
        c = covci.iloc[st:ed]
        emax = reduce(max,c['eidmax'].values, 0)
        emin = reduce(min,c['eidmin'].values, c['eidmax'].max())
        nc = len(c)
        ne = emax+1-emin
        mat = N.zeros((nc,ne))
        for i,n1 in enumerate(c['name1'].values):# fill in rows
            N.put(mat[i], N.array(n1)-emin, 1)
        # solve by NNLS (nonnegative least square)
        ecov,err = nnls(mat, c['cov'].values)
        e2c.update(dict(zip(range(emin,emax+1),ecov)))# e -> cov
    def catcherr(stcid, blocksize):
        bs = blocksize
        done = False
        while((bs>0)&(not done)):
            edcid = _find_chunk(stcid,bs)
            try:
                _calc_ecov(stcid,edcid)
                done = True
            except RuntimeError:
                bs = bs/2
        if bs==0: # unsuccessfull  set to zero
            e2c[stcid] = 0.
            LOG.warning('nnls did not converge for {0}, setting to 0'.format(stcid))
            return stcid+1
        return edcid
    # iterate
    stcid=0
    while(stcid<len(covci)):
        #edcid = _find_chunk(stcid,blocksize)
        #_calc_ecov(stcid,edcid)
        #stcid = edcid
        stcid = catcherr(stcid, blocksize)
    return e2c

def calc_ecov_mp(covci,fname,np,blocksize=100):
    """
    WARNING: this assumes _id is assinged according to sorted (chr,st,ed)
    """
    LOG.debug('calc_ecov...')
    chroms = sorted(covci['chr'].unique())
    if 'name1' not in covci.columns:
        covci['name1'] = covci['name'].astype(str).apply(lambda x: [int(y) for y in x.split(',')])
    if 'eidmax' not in covci.columns:
        covci['eidmax'] = covci['name1'].apply(lambda x:max(x))
    if 'eidmin' not in covci.columns:
        covci['eidmin'] = covci['name1'].apply(lambda x:min(x))
    args = [(covci[covci['chr']==c].copy(), blocksize) for c in chroms]
    e2cs = {}
    if np==1:
        # for c,bwname,chrom,d in data:
        for arg in args:
            e2cs.update(calc_ecov_chrom(*arg))
    else:
        try:
            p = multiprocessing.Pool(np)
            rslts = p.map(mp_worker, zip(repeat(calc_ecov_chrom), args))
        finally:
            LOG.debug('closing pool')
            p.close()
        for x in rslts:
            e2cs.update(x)
    LOG.debug('writing rslts...')
    if fname is None:
        return e2cs
    ccf = UT.flattendf(covci, 'name1')
    ccfg = ccf.groupby('name1')
    e2chr = dict(UT.izipcols(ccfg['chr'].first().reset_index(), ['name1','chr']))
    e2st = dict(UT.izipcols(ccfg['st'].min().reset_index(), ['name1','st']))
    e2ed = dict(UT.izipcols(ccfg['ed'].max().reset_index(), ['name1','ed']))
    df = PD.DataFrame(e2cs,index=['ecov']).T
    df.index.name='eid'
    df = df.reset_index()
    df['chr'] = [e2chr[x] for x in df['eid']]
    df['st'] = [e2st[x] for x in df['eid']]
    df['ed'] = [e2ed[x] for x in df['eid']]
    UT.save_tsv_nidx_whead(df[['eid','chr','st','ed','ecov']], fname)
    return df

### gcov, gmax calc low level ##########################################

def worker_cov(c,bwname,chrom, path):
    it = BW.block_iter(bwname, chrom)
    recs = calc_cov_chrom(c,it)
    return recs

def worker_max(c,bwname,chrom, path):
    it = BW.block_iter(bwname, chrom)
    recs = calc_max_chrom(c,it)
    return recs
        
def calc_cov_mp(bed, bwname, fname, np, which='cov'):
    if which=='cov':
        worker=worker_cov
    elif which=='max':
        worker=worker_max

    if UT.isstring(bed):
        bed = GGB.read_bed(bed)
    #cols = list(bed.columns)+['cov']
    cols = list(bed.columns)+[which]
    chroms = bed['chr'].unique()
    #LOG.debug(chroms)
    cdir = os.path.dirname(__file__)
    data = [(bed[bed['chr']==c].copy(), bwname, c, cdir) for c in chroms]
    recs = []
    if np==1:
        # for c,bwname,chrom,d in data:
        for arg in data:
            LOG.debug('cov calculation: processing {0}...'.format(arg[-2]))
            recs += worker(*arg)
    else:
        LOG.debug('{1} calculation: np={0}'.format(np,which))
        try:
            p = multiprocessing.Pool(np)
            a = zip(repeat(worker), data)
            rslts = p.map(mp_worker, a)
            for v in rslts:
                recs += v
            LOG.debug('done {1} calculation: np={0}'.format(np,which))
        finally:
            LOG.debug('closing pool')
            p.close()
            #p.join()
        #recs = reduce(iadd, rslts)
    LOG.debug('writing rslts...')
    df = PD.DataFrame(recs,columns=cols) 
    UT.save_tsv_nidx_whead(df, fname)
    return df

def calc_max_chrom(c, b_iter):
    nomore = False
    try:
        start,end,value = next(b_iter)
    except StopIteration:
        nomore = True
    prev = None
    rslts = [None]*len(c)
    for i, row in enumerate(c.values): #c.iterrows():
        row = list(row)
        st,ed = row[1],row[2]
        # advance until intersect
        if not nomore:
            while(end<st):
                if prev is None:
                    try:
                        start,end,value=next(b_iter)
                    except StopIteration:
                        nomore = True
                        break
                else:
                    start,end,value = curr
                    prev = None  
        if nomore:
            row += [0.]
            rslts[i] = row
        else:
            # once intersect collect values
            v = 0.
            prev = (start,end,value) # remember one before
            while(start<ed):
                st0 = max(start,st-1)
                ed0 = min(end,ed)
                #v += value*(ed0 - st0)
                v = max(v, value)
                prev = (start,end,value) # remember one before
                try:
                    start,end,value=next(b_iter)
                except StopIteration:
                    nomore=True
                    break
            #row += [v/float(ed-st+1)]
            row += [v]
            rslts[i] = row
            curr = start,end,value
            start,end,value = prev
    return rslts

def calc_cov_chrom(c, b_iter):
    nomore = False
    try:
        start,end,value = next(b_iter)
    except StopIteration:
        nomore = True
    prev = None
    rslts = [None]*len(c)
    for i, row in enumerate(c.values): #c.iterrows():
        row = list(row)
        st,ed = row[1],row[2]
        # advance until intersect
        if not nomore:
            while(end<=st): # equal boundary important, otherwise if end==st it will skip one interval
                if prev is None:
                    try:
                        start,end,value=next(b_iter)
                    except StopIteration:
                        nomore = True
                        break
                else:
                    start,end,value = curr
                    prev = None
                    
        if nomore:
            row += [0.]
            rslts[i] = row
        else:
            # once intersect collect values
            v = 0.
            prev = (start,end,value) # remember one before
            while(start<ed):
                st0 = max(start,st)
                ed0 = min(end,ed)
                v += value*(ed0 - st0)
                prev = (start,end,value) # remember one before
                try:
                    start,end,value=next(b_iter)
                except StopIteration:
                    nomore=True
                    break
            row += [v/float(ed-st)]
            rslts[i] = row
            curr = start,end,value
            start,end,value = prev

    return rslts

### high level        ##################################################

def calc_sjcnt(sjpath1, sjpath2, dstprefix, override=False, np=4):
    pass

def calc_ecov(expath, cipath, bwpath, dstprefix, override=False, np=4):
    """Calculate exon coverages.

    Args:
        expath: merged ex
        cipath: chopped interval for ex
        bwpath: bigwig file (sample)
        dstprefix: prefix for outputs

    Outputs:
        1. dstprefix+'.covci.txt.gz': coverage for ci
        2. dstprefix+'.ecov.txt.gz' : DataFrame(cols: eid, chr, st, ed, ecov)

    """
    covcipath = dstprefix+'.covci.txt.gz'
    ecovpath = dstprefix+'.ecov.txt.gz'

    if UT.notstale([expath, cipath], covcipath, override):
        cc = UT.read_pandas(covcipath)
    else:
        if UT.notstale(expath, cipath, override):
            ci = UT.read_pandas(cipath, names=['chr','st','ed','name','id'])
        else:
            ex = UT.read_pandas(expath)
            ci = UT.chopintervals(ex, cipath, idcol='_id')
        cc = calc_cov_mp(ci, bwpath, covcipath, np=np)

    # if override or (not os.path.exists(covcipath)):
    #     # calc covci
    #     if not os.path.exists(cipath):
    #         ex = UT.read_pandas(expath)
    #         ci = UT.chopintervals(ex, cipath, idcol='_id')
    #     else:
    #         ci = UT.read_pandas(cipath, names=['chr','st','ed','name','id'])
    #     cc = calc_cov_mp(ci, bwpath, covcipath, np=np)
    # else:
    #     cc = UT.read_pandas(covcipath)

    if 'id' not in cc.columns:
        cc['id'] = cc['sc1']
    if 'eid' not in cc.columns:
        cc['eid'] = cc['name'].astype(str).apply(lambda x: [int(y) for y in x.split(',')])
        cc['name1'] = cc['eid']
    ccf = UT.flattendf(cc[['chr','st','ed','eid']], 'eid')
    ccfg = ccf.groupby('eid')
    df = ccfg[['chr']].first()
    df['st'] = ccfg['st'].min()
    df['ed'] = ccfg['ed'].max()
    df.reset_index(inplace=True)
    e2cs = calc_ecov_mp(cc, None, np)
    df['ecov'] = [e2cs[x] for x in df['eid']]
    UT.save_tsv_nidx_whead(df, ecovpath)
    return df

def calc_gcov(expath, cipath, bwpath, dstprefix, override=False, np=4):
    """Calculate gene coverages.

    Args:
        expath: merged ex
        cipath: chopped interval for ex
        bwpath: bigwig file (sample)
        dstprefix: prefix for outputs

    Outputs:
        1. dstprefix+'.covci.txt.gz'
        2. dstprefix+'.gcov.txt.gz' : DataFrame(col:_gidx,len,val,gcov,len2,gcov2,cids)
            len2: calculate length from ci with cov > 0
            (normal length = use entire ci's belonging to the gene)
            gcov2 = val/len2
            cids: cid with cov > for the gene ','.joined
    """
    ex = UT.read_pandas(expath)
    covcipath = dstprefix+'.covci.txt.gz'
    gcovpath = dstprefix+'.gcov.txt.gz'

    if UT.notstale([expath, cipath], covcipath, override):
        cc = UT.read_pandas(covcipath)
    else:
        if UT.notstale(expath, cipath, override):
            ci = UT.read_pandas(cipath, names=['chr','st','ed','name','id'])
        else:
            ci = UT.chopintervals(ex, cipath, idcol='_id')
        cc = calc_cov_mp(ci, bwpath, covcipath, np=np)

    # if override or (not os.path.exists(covcipath)):
    #     # calc covci
    #     if not os.path.exists(cipath):
    #         ci = UT.chopintervals(ex, cipath, idcol='_id')
    #     else:
    #         ci = UT.read_pandas(cipath, names=['chr','st','ed','name','id'])
    #     cc = calc_cov_mp(ci, bwpath, covcipath, np=np)
    # else:
    #     cc = UT.read_pandas(covcipath)

    if 'id' not in cc.columns:
        cc['id'] = cc['sc1']
    if 'eid' not in cc.columns:
        cc['eid'] = cc['name'].astype(str).apply(lambda x: [int(y) for y in x.split(',')])
    cc['len'] = cc['ed']-cc['st']
    cc['val'] = cc['cov']*cc['len']
    ccf = UT.flattendf(cc[['id','eid','len','val','st','ed']], 'eid')
    e2g = dict(UT.izipcols(ex, ['_id','_gidx']))
    ccf['_gidx'] = [e2g[x] for x in ccf['eid']]
    # for normal gcov: take unique combination of (gid, id) (id=cid)
    # for gocv2 : first select ccf with val>0
    ccf2 = ccf[ccf['val']>0].groupby(['_gidx','id']).first().reset_index()
    ccf2g = ccf2.groupby('_gidx')
    df2 = ccf2g[['len','val']].sum()
    df2['gcov2'] = df2['val']/df2['len']
    df2['cids'] = ccf2g['id'].apply(lambda x: ','.join([str(y) for y in x]))
    df2['gst2'] = ccf2g['st'].min()
    df2['ged2'] = ccf2g['ed'].max()
    df2['glen2'] = df2['ged2']-df2['gst2']

    df2 = df2.reset_index()

    ccf1 = ccf.groupby(['_gidx','id']).first().reset_index()
    ccf1g = ccf1.groupby('_gidx')
    df = ccf1g[['len','val']].sum()
    df['gcov'] = df['val']/df['len']
    df['st'] = ccf1g['st'].min()
    df['ed'] = ccf1g['ed'].max()
    df['glen'] = df['ed'] - df['st']
    df = df.reset_index()
    g2chr = dict(UT.izipcols(ex, ['_gidx','chr']))
    df['chr'] = [g2chr[x] for x in df['_gidx']]

    def _set_df2prop(src,tgt, default):
        dic = dict(UT.izipcols(df2, ['_gidx',src]))
        df[tgt] = [dic.get(x,default) for x in df['_gidx']]
    _set_df2prop('gcov2','gcov2', 0)
    _set_df2prop('len','len2', 0)
    _set_df2prop('cids','cids', '')
    _set_df2prop('gst2','st2', -1)
    _set_df2prop('ged2','ed2', -1)
    _set_df2prop('glen2','glen2', 0)

    cols = ['_gidx','chr','st','ed','len','val','gcov','glen','len2','gcov2','cids','st2','ed2','glen2']
    df = df[cols]
    UT.save_tsv_nidx_whead(df, gcovpath)
    return df

# just use trimed ex to calculate gcov using calc_gcov
# def calc_gcov1000(expath, cipath, bwpath, dstprefix, override=False, np=4):
    # """Calculate gene coverage but only use 1000bp from 3' end.

    # Args:
    #     expath: merged ex
    #     cipath: chopped interval for ex
    #     bwpath: bigwig file (sample)
    #     dstprefix: prefix for outputs
    #     1. dstprefix+'.covci.txt.gz'
    #     2. dstprefix+'.gcov.txt.gz' : DataFrame(col:_gidx,len,val,gcov,cids)
    #         cids: cid with cov > for the gene ','.joined
    # """
    # ex = UT.read_pandas(expath)
    # covcipath = dstprefix+'.covci.txt.gz'
    # gcovpath = dstprefix+'.gcov1000.txt.gz'
    # if override or (not os.path.exists(covcipath)):
    #     # calc covci
    #     if not os.path.exists(cipath):
    #         ci = UT.chopintervals(ex, cipath, idcol='_id')
    #     else:
    #         ci = UT.read_pandas(cipath, names=['chr','st','ed','name','id'])
    #     cc = calc_cov_mp(ci, bwpath, covcipath, np=np)
    # else:
    #     cc = UT.read_pandas(covcipath)

    # if 'id' not in cc.columns:
    #     cc['id'] = cc['sc1']
    # if 'eid' not in cc.columns:
    #     cc['eid'] = cc['name'].astype(str).apply(lambda x: map(int, x.split(',')))
    # cc['len'] = cc['ed']-cc['st']
    # cc['val'] = cc['cov']*cc['len']
    # ccf = UT.flattendf(cc[['id','eid','len','val','st','ed']], 'eid')
    # e2g = dict(UT.izipcols(ex, ['_id','_gidx']))
    # ccf['_gidx'] = [e2g[x] for x in ccf['eid']]
    # # for normal gcov: take unique combination of (gid, id) (id=cid)
    # ccf1 = ccf.groupby(['_gidx','id']).first().reset_index()

    # # MODIFICATION to only use <cumsum 2000bp 
    # g2strand = dict(UT.izipcols(ex.groupby('_gidx').first().reset_index(), ['_gidx','strand']))
    # ccf1['strand'] = [g2strand[x] for x in ccf1['_gidx']]
    # ccf1 = ccf1.sort_values(['_gidx','ed','st'],ascending=False)
    # ccf1['cumsum+'] = ccf1.groupby('_gidx')['len'].cumsum()
    # ccf1 = ccf1.sort_values(['_gidx','st','ed'],ascending=True)
    # ccf1['cumsum-'] = ccf1.groupby('_gidx')['len'].cumsum()
    # ccf1 = ccf1.sort_values(['_gidx','ed','st'],ascending=False)
    
    # size1p = ccf1.groupby('_gidx')['cumsum+'].first().to_frame('v').reset_index()
    # g2s1p = dict(UT.izipcols(size1p, ['_gidx','v']))

    # ccf1 = ccf1.sort_values(['_gidx','st','ed'],ascending=True)
    # size1n = ccf1.groupby('_gidx')['cumsum-'].first().to_frame('v').reset_index()
    # g2s1n = dict(UT.izipcols(size1n, ['_gidx','v']))

    # ccf1['size1+'] = [g2s1p[x] for x in ccf1['_gidx']]
    # ccf1['size1-'] = [g2s1n[x] for x in ccf1['_gidx']]

    # idx =(((ccf1['strand']=='+')&((ccf1['cumsum+']<2000)|(ccf1['size1+']>=2000)))|
    #       ((ccf1['strand']!='+')&((ccf1['cumsum-']<2000)|(ccf1['size1-']>=2000))))
    # ccf2 = ccf1[idx]
    # ccf1g = ccf2.groupby('_gidx')


    # df = ccf1g[['len','val']].sum()
    # df['gcov'] = df['val']/df['len']
    # df['st'] = ccf1g['st'].min()
    # df['ed'] = ccf1g['ed'].max()
    # df['glen'] = df['ed'] - df['st']
    # df['cids'] = ccf1g['id'].apply(lambda x: ','.join(map(str, x)))

    # df = df.reset_index()
    # g2chr = dict(UT.izipcols(ex, ['_gidx','chr']))
    # df['chr'] = [g2chr[x] for x in df['_gidx']]
    
    # cols = ['_gidx','chr','st','ed','len','val','gcov','glen','cids']
    # df = df[cols]
    # UT.save_tsv_nidx_whead(df, gcovpath)
    # return df

# old method
# def calc_gene_cov(ex, cc, gidfld='_gidx'):
    # """
    # ex: exon df
    # cc: covci df
    # return dfg ['val','len','cov']
    # """
    # if 'id' not in cc.columns:
    #     cc['id'] = cc['sc1']
    # if 'cid' not in ex.columns:
    #     e2c = {}
    #     for i,name in cc[['id','name']].values:
    #         for eid in map(int, name.split(',')):
    #             e2c.setdefault(eid,[]).append(i)
    #     ex['cid'] = [e2c[x] for x in ex['_id']]
    # # flatten
    # def _gen():
    #     for g,cids in UT.izipcols(ex, [gidfld,'cid']):
    #         for c in cids:
    #             yield (c,g)
    # df = PD.DataFrame(list(set([x for x in _gen()])),columns=['cid',gidfld])
    # cc['len'] = cc['ed'] - cc['st']
    # cc['val'] = cc['cov']*cc['len']
    # c2len = dict(cc[['id','len']].values)
    # c2val = dict(cc[['id','val']].values)
    # df['len'] = [c2len[x] for x in df['cid']]
    # df['val'] = [c2val[x] for x in df['cid']]
    # dfg = df.groupby(gidfld)[['val','len']].sum()
    # dfg['cov'] = dfg['val']/dfg['len']
    # return dfg
 
def calc_glen(ex, cipath):
    ci = GGB.read_bed(cipath) # 5 col bed, name:eids, sc1:cid
    ci['len'] = ci['ed']-ci['st']
    ci['cid'] = ci['sc1']
    c2l = dict(UT.izipcols(ci, ['cid','len']))
    if 'cid' not in ex.columns:
        e2c = {}
        for i,name in ci[['cid','name']].values:
            for eid in name.split(','):
                e2c.setdefault(int(eid),[]).append(i)
        ex['cid'] = [e2c[x] for x in ex['_id']]
    def _gen():
        for g,cids in UT.izipcols(ex, ['_gidx','cid']):
            for c in cids:
                yield (c,g)
    df = PD.DataFrame(list(set([x for x in _gen()])),columns=['cid','_gidx'])
    df['len'] = [c2l[x] for x in df['cid']]
    glen = df.groupby('_gidx')['len'].sum()
    return dict(zip(glen.index, glen.values))
    

