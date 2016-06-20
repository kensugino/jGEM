"""

.. module:: gtfgffbed
    :synopsis: GTF/GFF/BED related functions

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""

import csv
import subprocess
import os
import gzip
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)
from itertools import repeat, groupby

import pandas as PD
import numpy as N

from jgem import utils as UT
import jgem.cy.bw as cybw


GTFCOLS = [u'chr',u'src',u'typ',u'st',u'ed',u'sc1',u'strand',u'sc2',u'extra']
GFFCOLS = [u'chr',u'src',u'typ',u'st',u'ed',u'sc1',u'strand',u'sc2',u'attr']
BEDCOLS = [u'chr', u'st', u'ed', u'name', u'sc1', u'strand', u'tst', u'ted', u'sc2', u'#exons', u'esizes', u'estarts']
KGCOLS = [u'name',u'chr',u'strand',u'st',u'ed',u'cst',u'ced',u'excnts',u'exstarts',u'exends',u'proteinID',u'alignID']
RGCOLS = [u'bin',u'name',u'chr',u'strand',u'st',u'ed',u'cst',u'ced',u'excnts',u'exstarts',u'exends',
            u'score',u'name2',u'cdsStartStat',u'cdsEndStat', u'exonFrames']

SJCOLS = [u'chr', u'st', u'ed', u'name', u'ucnt', u'strand', u'mcnt']

SJTABCOLS = [u'chr',u'st',u'ed',u'strand2',u'motif',u'annotated',u'ureads',u'mreads',u'maxoverhang']
SJTABMOTIF = {0:'non-canonical',1:'GT/AG',2:'CT/AC',3:'GC/AG',4:'CT/GC',5:'AT/AC',6:'GT/AT'}
SJTABSTRAND = {1:'+',2:'-',0:'.'}

DEFAULT_GTF_PARSE = [u'gene_id',u'transcript_id',u'exon_number',u'gene_name',u'cov',u'FPKM']

# SJ.out.tab to SJBED ###################################################################

def sjtab2sjbed(sjtab, sjbed, scale=None):
    """Generate splice junction input file from STAR SJ.out.tab

    Args:
        sjtab (str): path to SJ.out.tab file
        sjbed (str): path to output bed
        scale (float): a number to multiply (default None, no change)

    Returns:
        Pandas dataframe

    """
    sj = PD.read_table(sjtab, names=SJTABCOLS)
    sj['name'] = ['%s-k%d-u%d-m%d-o%d' % (SJTABMOTIF[x], a, u, m, o) for x,a, u,m,o in \
                  sj[['motif','annotated','ureads','mreads','maxoverhang']].values]
    sj['strand'] = [SJTABSTRAND[x] for x in sj['strand2']]
    #scale = 1e6/float(aligned)
    if scale is None:
        sj['ucnt'] = sj['ureads']
        sj['mcnt'] = sj['mreads']
    else:
        sj['ucnt'] = sj['ureads']*scale
        sj['mcnt'] = sj['mreads']*scale
    #sj['jcnt'] = [x or y for x, y in sj[['ucnt','mcnt']].values]
    #cols = ['chr','st','ed','name','strand','ucnt','mcnt']#,'jcnt']
    #UT.write_pandas(sj[cols], sjbed, '')
    sj['sc1'] = sj['ucnt'] #sj['ureads']*scale
    sj['tst'] = sj['mcnt'] #sj['mreads']*scale
    #cols = ['chr','st','ed','name','sc1','strand','tst'] 
    #UT.write_pandas(sj[cols], sjbed, '')
    write_bed(sj, sjbed, ncols=7)
    return sj

def read_sj(path, parsename=False):
    # read BED (input) or TXT (output) with consistent column names
    if path[-7:]=='.bed.gz' or path[-4:]=='.bed':
        df = read_bed(path).rename(columns={'sc1':'ucnt','tst':'mcnt'})
        if parsename:
            # name is encoded as above 'motif-k0[k1]-u(reads)-m(reads)-o(maxoverhang)'
            # motif(0), known(1), u(2), m(3), o(4)
            tmp = df['name'].str.split('-')
            df['motif'] = tmp.str[0]
            df['annotated'] = tmp.str[1].str[1]
            df['maxoverhang'] = tmp.str[4].str[1:].astype(int)
    else:
        df = UT.read_pandas(path) # header should be there
    return df


# READ/WRITE       ######################################################################    
def get_gff_attr_col(gff, aname):
    "helper function for read_gff"
    return [dict([y.split('=') for y in line.split(';') if '=' in y]).get(aname,'') for line in gff['attr']]

def read_gff(gffname, onlytypes=[], parseattrs=[]):
    """ Read in whole GFF, parse id & parent

    Args:
        gffname: path to GFF file
        onlytypes: only keep these types. If [] or None, then keep all (default). 
        parseattrs: extra attributes to parse

    Returns:
        Pandas.DataFrame containing GFF data
    
    """
    if not UT.isstring(gffname):
        return gffname
        
    if gffname.endswith('.gz'):
        gff = PD.read_table(gffname, names=GFFCOLS, comment='#', compression='gzip')
    else:
        gff = PD.read_table(gffname, names=GFFCOLS, comment='#')

    for c in ['ID','Parent']+parseattrs:
        gff[c] = get_gff_attr_col(gff, c)

    # set gid, tid, eid by default
    gff['gid'] = ''
    gff['tid'] = ''
    gff['eid'] = ''
    # genes
    gidx = gff['typ']=='gene'
    gff.ix[gidx, 'gid'] = gff['ID']
    # transcripts
    tidx = gff['typ']=='transcript'
    gff.ix[tidx, 'gid'] = gff['Parent']
    gff.ix[tidx, 'tid'] = gff['ID']
    # exons
    tid2gid = dict(gff[tidx][['tid','gid']].values)
    eidx = gff['typ']=='exon'
    gff.ix[eidx, 'tid'] = gff['Parent']
    gff.ix[eidx, 'gid'] = [tid2gid[x] for x in gff[eidx]['tid'].values]
    if N.sum(gff[eidx]['ID']=='')==0: # ID already set
        gff.ix[eidx, 'eid'] = gff['ID']
    else: # sometimes exons don't have ID => create
        edf = gff[eidx]
        en = edf.groupby('tid')['gid'].transform(lambda x: N.arange(1,len(x)+1))
        eids = edf['tid']+':'+ en.astype(str)
        gff.ix[eidx, 'eid'] = eids
        attr = 'ID='+eids+';Parent='+edf['Parent']
        gff.ix[eidx, 'attr'] = attr

    if onlytypes:
        gff = gff[gff['typ'].isin(onlytypes)]

    return gff

# this function is awfully slow => replace with vectorized version
def get_gtf_attr_col(gtf, aname):
    "helper function for read_gtf"
    def _attr(line):
        #if isinstance(line, basestring):
        #if isinstance(line, str):
        if UT.isstring(line):
            dic = dict([(x[0],x[1][1:-1]) for x in [y.split() for y in line.split(';')] if len(x)>1])
            return dic.get(aname,'')
        return ''
    return [_attr(line) for line in gtf['extra']]

# ~35 sec to read Gencode.vM4
def read_gtf(gtfname, onlytypes=[], parseattrs=DEFAULT_GTF_PARSE, rename={}):
    """ Read in whole GTF, parse gene_id, transcript_id from column 9

    Args:
        gtfname: path to GTF file
        onlytypes: only keep these types. If [] or None, then keep all (default).
        parseattrs: which column attributes to parse.

    Returns:
        Pandas DataFrame containing GTF data

    """
    recs,cols = cybw.read_gtf_helper(gtfname, parseattrs, '#')
    if len(recs)==0 or len(recs[0])!=len(cols):
        return UT.make_empty_df(cols)
    df = PD.DataFrame(recs, columns=cols)
    df['st'] = df['st'].astype(int)
    df['ed'] = df['ed'].astype(int)
    df.replace('', N.nan, inplace=True)
    if onlytypes:
        df = df[df['typ'].isin(onlytypes)].copy()
    if rename:
        df.rename(columns=rename, inplace=True)    
    return df

# ~3.5 min to read Gencode.vM4
def read_gtf1(gtfname, onlytypes=[], parseattrs=DEFAULT_GTF_PARSE, rename={}):
    """ Read in whole GTF, parse gene_id, transcript_id from column 9

    Args:
        gtfname: path to GTF file
        onlytypes: only keep these types. If [] or None, then keep all (default).
        parseattrs: which column attributes to parse.

    Returns:
        Pandas DataFrame containing GTF data

    """
    if not UT.isstring(gtfname):
        return gtfname
        
    if gtfname.endswith('.gz'):
        gtf = PD.read_table(gtfname, names=GTFCOLS, compression='gzip', comment='#')
    else:
        gtf = PD.read_table(gtfname, names=GTFCOLS, comment='#')
    if onlytypes:
        gtf = gtf[gtf['typ'].isin(onlytypes)].copy()
    LOG.debug( "extracting ids...")
    # field1 "field1 value"; field2 "field2 value"; ...
    # 14.283 sec (using get_gtf_attr_col) ==> 5.830 sec (using below)
    tmp = gtf['extra'].str.split(';',expand=True) # each column: field "field val"
    cols = gtf[['chr']].copy()
    for c in tmp.columns:
        kv = tmp[c].str.split(expand=True) # key col(0) and value col(1)
        if kv.shape[1]<2: # split unsuccessful
            LOG.debug('column {0} did not split'.format(c))
            LOG.debug('col:{0} first line:{1}, kv.shape:{2}'.format(c, tmp[c].iloc[0], kv.shape))
            continue
        #LOG.debug((kv[0].unique(),kv.shape, kv.columns))
        LOG.debug(kv[0].unique())
        if len(kv[0].unique())==1:
            k = kv.iloc[0][0]
            if UT.isstring(k) and kv.shape[1]==2:
                cols[k] = kv[1].str.replace('"','') # strip "
        else: # multiple  fields => make cols for each
            for k in kv[0].unique():
                if UT.isstring(k):
                    idx = kv[0]==k
                    cols.loc[idx, k] = kv[1][idx].str.replace('"','')

    for c in parseattrs:
        if c in cols:
            gtf[c] = cols[c] # get_gtf_attr_col(gtf, c)
        else:
            LOG.warning('column {0} not found'.format(c))
    if rename:
        gtf.rename(columns=rename, inplace=True)
    return gtf

# ~3.5 min to read Gencode.vM4 (original version)
def read_gtf0(gtfname, onlytypes=[], parseattrs=DEFAULT_GTF_PARSE, rename={}):
    """ Read in whole GTF, parse gene_id, transcript_id from column 9

    Args:
        gtfname: path to GTF file
        onlytypes: only keep these types. If [] or None, then keep all (default).
        parseattrs: which column attributes to parse.

    Returns:
        Pandas DataFrame containing GTF data

    """
    if not UT.isstring(gtfname):
        return gtfname
        
    if gtfname.endswith('.gz'):
        gtf = PD.read_table(gtfname, names=GTFCOLS, compression='gzip', comment='#')
    else:
        gtf = PD.read_table(gtfname, names=GTFCOLS, comment='#')
    if onlytypes:
        gtf = gtf[gtf['typ'].isin(onlytypes)].copy()
    LOG.debug( "extracting ids...")
    # field1 "field1 value"; field2 "field2 value"; ...
    for c in parseattrs:
        gtf[c] = get_gtf_attr_col(gtf, c)
    if rename:
        gtf.rename(columns=rename, inplace=True)
    return gtf

def read_bed(fpath):
    """Read BED file

    Args:
        fpath: path to BED file (no header)

    Returns:
        Pandas DataFrame containing BED data

    """
    if not UT.isstring(fpath):
        return fpath

    if fpath.endswith('.gz'):
        d = PD.read_table(fpath, header=None, compression='gzip')
        d.columns = BEDCOLS[:len(d.columns)]
    else:
        d = PD.read_table(fpath, header=None)
        d.columns = BEDCOLS[:len(d.columns)]
    # if calcextra:
    #     d['locus'] = d['chr'].astype(str) + ':'+ d['st'].astype(str)+'-'+ d['ed'].astype(str)
    #     d['min.exon.size'] = d['esizes'].apply(lambda x: N.min(list(map(int, x[:-1].split(',')))))
    #     d['max.exon.size'] = d['esizes'].apply(lambda x: N.max(list(map(int, x[:-1].split(',')))))
    #     d['length'] = d['ed']-d['st']
    # make sure st,ed are integers
    # discard entries with NaN in (chr,st,ed)

    idx = (d['chr'].isnull())|(d['st'].isnull())|(d['ed'].isnull())
    if N.sum(idx)>0:
        LOG.warning('{1} NaN in chr/st/ed in file {0}, discarding'.format(fpath, N.sum(idx)))
        d = d[~idx].copy()
    d['st'] = d['st'].astype(int)
    d['ed'] = d['ed'].astype(int)
    return d

def write_gff(df, fname):
    """Write GFF file.

    Args:
        df: Pandas.DataFrame containing GFF data
        fname: path (if ends with .gz, gzipped)

    Returns:
        actual path written
    """
    return write_ggb(df, fname, GFFCOLS)
    
def write_gtf(df, fname):
    """Write GTF file.

    Args:
        df: Pandas.DataFrame containing GTF data
        fname: path (if ends with .gz, gzipped)

    Returns:
        actual path written
    """
    return write_ggb(df, fname, GTFCOLS)
    
def write_bed(df, fname, ncols=None, mode='w'):
    """Write BED file.

    Args:
        df: Pandas.DataFrame containing BED data
        fname: path (if ends with .gz, gzipped)
        ncols: number of bed columns (default 12)

    Returns:
        actual path written
    """ 
    if ncols is None:
        for i, x in enumerate(BEDCOLS):
            if x not in df.columns:
                ncols = i
                break
        else:
            ncols = 12
    return write_ggb(df, fname, BEDCOLS[:ncols], mode=mode)
    
def write_ggb(df, fname, cols, mode='w'):    
    # df.loc[:,'st'] = df['st'].astype(int)
    # df.loc[:,'ed'] = df['ed'].astype(int)
    if fname[-3:]=='.gz':
        compress=True
        fname = fname[:-3]
    else:
        compress=False
    if (df.dtypes['st'] != int) or (df.dtypes['ed'] != int):
        LOG.warning('st,ed not integer: copy and converting')
        df = df.copy()
        df['st'] = df['st'].astype(int)
        df['ed'] = df['ed'].astype(int)
    UT.makedirs(os.path.dirname(fname))
    with open(fname, mode) as f:
        df[cols].to_csv(f, index=False, header=False, sep='\t', quoting=csv.QUOTE_NONE)
    if compress:
        return UT.compress(fname)
    return fname

# UCSC knownGene/refGene ##############################################################

def read_ucsc_knownGene(path):
    return UT.read_pandas(path, names=KGCOLS)

def read_ucsc_refGene(path):
    return UT.read_pandas(path, names=RGCOLS)


# CONVERSION     ######################################################################

def gtf2gff(gtfname,gffname, memt=True):
    """Convert GTF to GFF.

    Args:
        gtfname: path to GTF file
        gffname: path for converted GFF file
        memt: only select multiexon, multitranscript

    Returns:
        Pandas.DataFrame containing converted GFF data
    """
    eids = read_gtf(gtfname, 
                    onlytypes=['exon'], 
                    parseattrs=['gene_id','transcript_id','exon_number','gene_name'],
                    rename={'gene_id':'gid','transcript_id':'tid','gene_name':'gname','exon_number':'e#'})
    if N.sum(eids['e#']=='')>0: # recalculate exon_number
        eids['e#'] = eids.groupby('tid')['gid'].transform(lambda x: N.arange(1,len(x)+1))
    else:
        eids['e#'] = eids['e#'].astype(int)
    eids['ID'] = eids['tid']+':'+eids['e#'].astype(str)
    eids['attr'] = 'ID='+eids['ID']+';Parent='+eids['tid']

    # groupby tid and get transcript records
    LOG.debug( "calculating transcripts...")
    gtid = eids.groupby('tid')
    tids = gtid.first().copy() # in general get first record
    tids['typ'] = 'transcript'  # fix typ
    tids['st'] = gtid['st'].min() # fix st
    tids['ed'] = gtid['ed'].max() # fix ed
    tids['#exons'] = gtid.size()
    if memt:
        tids = tids[tids['#exons']>1]
    tids = tids.reset_index()
    tids['e#'] = 0
    tids['attr'] = 'ID='+ tids['tid']+';Parent='+tids['gid']+\
                   ';num_exons='+tids['#exons'].astype(str)+\
                   ';gene_name='+tids['gname']

    # groupby gid and get gene records
    LOG.debug( "calculating genes...")
    ggid = tids.groupby('gid')
    gids = ggid.first().copy()
    gids['typ'] = 'gene'
    gids['st'] = ggid['st'].min()
    gids['ed'] = ggid['ed'].max()
    gids['#trans'] = ggid.size()
    gids = gids.reset_index()
    if memt:
        gids = gids[gids['#trans']>1] # multi transcript
        tids = tids[tids['gid'].isin(gids['gid'].values)]
        eids = eids[eids['tid'].isin(tids['tid'].values)]
    gids['tid'] = ''
    gids['e#'] = -1
    gids['attr'] = 'ID='+gids['gid']+';num_trans='+gids['#trans'].astype(str)

    LOG.debug( "merging exons, transcripts, genes...")
    gte = PD.concat([gids,tids,eids],ignore_index=True)
    # sort by gid,tid,st,ed
    gte = gte.sort_values(['chr','gid','tid','e#'])
    # write out
    LOG.debug( "writing GFF...")
    write_gff(gte, gffname)
    return gte

def gtf2bed12(fpath, compress=True):
    """Convert GTF to BED. Uses gtfToGenePred, genePredToBed (UCSC Kent Tools)

    Args:
        gtfname: path to GTF file
        compress: whether to gzip (default True)

    Returns:
        Pandas.DataFrame containing converted BED12 data
    """
    if fpath.endswith('.gz'):
        base = fpath[:-7]
        cmd = ['gunzip',fpath]
        LOG.debug( "expanding compressed ...{0}".format(base))
        subprocess.call(cmd)
    else:
        base = fpath[:-4]
    cmd = ['gtfToGenePred','-genePredExt','-ignoreGroupsWithoutExons',base+'.gtf',base+'.gp']
    LOG.debug( "converting to GenPred...{0}".format(base))
    ret = subprocess.call(cmd)
    if ret != 0:
        LOG.debug("error converting to GenPred...code {0}".format(ret))
        raise Exception
    cmd = ['genePredToBed', base+'.gp', base+'.bed']
    LOG.debug( "converting to Bed12...", base)
    ret = subprocess.call(cmd)
    if ret != 0:
        LOG.debug("error converting to GenPred...code {0}".format(ret))
        raise Exception
    os.unlink(base+'.gp')
    # gzip
    LOG.debug("gzipping ...{0}.bed".format(base))
    bdpath = base+'.bed'
    if compress:
        bdpath = UT.compress(bdpath)
    if fpath.endswith('.gz'):
        LOG.debug( "gzipping ...{0}".format(fpath[:-3]))
        p = subprocess.call(['gzip',fpath[:-3]])
        LOG.debug( "subprocess result: {0} ".format(p))
    return bdpath

def bed2gtf(fpath, compress=True):
    """Convert BED to GTF. Uses bedToGenePred, genePredToGtf (UCSC Kent Tools)

    Args:
        gtfname: path to BED file
        compress: whether to gzip (default True)

    Returns:
        Pandas.DataFrame containing converted GTF data
    """
    if fpath.endswith('.gz'):
        base = fpath[:-7]
        cmd = ['gunzip',fpath]
        LOG.debug( "expanding compressed ...", base)
        subprocess.call(cmd)
    else:
        base = fpath[:-4]
    gppath = base+'.genePred'
    bdpath = base+'.gtf'
    cmd = ['bedToGenePred',base+'.bed', gppath]
    LOG.debug( "converting to GenPred...", base)
    subprocess.call(cmd)
    cmd = ['genePredToGtf','-source=.','file', gppath, bdpath]
    LOG.debug( "converting to GTF...", base)
    subprocess.call(cmd)
    os.unlink(gppath)
    # gzip
    LOG.debug( "gzipping ...", bdpath)
    if compress:
        UT.compress(bdpath)
        bdpath=bdpath+'.gz'
    if fpath.endswith('.gz'):
        LOG.debug( "gzipping ...", fpath[:-3])
        subprocess.call(['gzip',fpath[:-3]])
    return bdpath

def bed12Tobed6(bed12):
    def _gen():
        for rec in bed12.values:
            esizes = [int(x) for x in rec[-2].split(',') if x!='']
            estarts = [int(x) for x in rec[-1].split(',') if x!='']
            for y,z in zip(estarts, esizes):
                st = rec[1]+y
                ed = st+z
                yield (rec[0],st,ed,rec[3],rec[4],rec[5])
    bed6 = PD.DataFrame([x for x in _gen()], names=['chr','st','ed','name','sc1','strand'])
    return bed6

def unionex2bed12(uex, gidx='_gidx', name='name', sc1='sc1', sc2='sc2'):
    """
    Args:
        uex: unionex
        gidx: colname for unique gene id
        name: colname for BED name field
        sc1: colname for BED sc1 field
        sc2: colname for BED sc2 field
    """
    cols0 = uex.columns
    if sc1 not in cols0:
        uex[sc1] = 0
    if sc2 not in cols0:
        uex[sc2] = 0
    if name not in cols0:
        uex[name] = uex[gidx]
    bed = uex[[gidx,'chr','st','ed',name,sc1,'strand',sc2]].sort_values([gidx,'chr','st','ed'])
    def _gen():
        for gid, exs in groupby(bed.values, lambda x: x[0]):
            exs = list(exs)
            nex = len(exs)
            st0 = exs[0][2]
            ed0 = exs[-1][3]
            esiz = ','.join([str(x[3]-x[2]) for x in exs])+','
            ests = ','.join([str(x[2]-st0) for x in exs])+','
            x = exs[0]
            yield (x[1],st0,ed0,x[4],x[5],x[6],st0,ed0,x[7],nex,esiz,ests)
    #bcols=['chr','st','ed','name','sc1','strand','tst','ted','sc1','#exons','esizes','estarts']
    bcols = BEDCOLS
    df = PD.DataFrame([x for x in _gen()], columns=bcols).sort_values(['chr','st','ed'])
    return df




# UTILS         ######################################################################

def chop_chrs_gtf(gtfname, chrs, outdir=None):
    """Separate chromosomes into different files.

    Args:
        gtfname: path to GTF
        chrs: list of chromosome names
        outdir: output directory, if None (default), then use same directory as input
        
    """
    #chrs = ['chr%d' % (x+1,) for x in range(19)] +['chrX','chrY']
    if outdir is None:
        outdir = os.path.dirname(gtfname)
    base = os.path.basename(gtfname)[:-4]
    outnames = [os.path.join(outdir, base+'-%s.gtf' % x) for x in chrs]
    if all([UT.notstale(gtfname, x) for x in outnames]):
        # all files already exist and newer than gtfname
        return outnames
    gtf = read_gtf(gtfname, parseattrs=[]) # don't parse attrs
    for c,fname in zip(chrs,outnames):
        LOG.debug( "writing %s to %s..." % (c, fname))
        sub = gtf[gtf['chr']==c]
        write_gtf(sub, fname, compress=False)
    return outnames    





