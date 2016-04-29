"""

.. module:: bw.pyx
    :synopsis: cython speedup of wiggle/bigwig related stuffs

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""

import numpy as N
cimport numpy as N
import gzip
import codecs
import os
import errno

def makedirs2(path):
    """Make all the directories along the path, calls os.makedirs but 
    intercept irrelevant exceptions
    """
    dirpath = os.path.dirname(path)
    try:
        os.makedirs(dirpath)
    except OSError as err:
        if err.errno != errno.EEXIST or not os.path.isdir(dirpath):
            raise


ctypedef N.float32_t F32_t
ctypedef N.float64_t F64_t
ctypedef N.int32_t I32_t
ctypedef N.ndarray ADTYPE_t
cimport cython

@cython.boundscheck(False) # turns off bounds-checking for entire function
@cython.wraparound(False) # turns off negative indexing checking
cpdef list union_intervals(ADTYPE_t[I32_t,ndim=2] a):
    """ input: 2d int32 array of intervals, assume sorted
    returns union intervals
    """
    cdef list rows = []
    cdef int cst,ced,st,ed
    cdef Py_ssize_t n = len(a)
    cst,ced = a[0] # first
    for i in range(1,n):
        st,ed = a[i]
        if st>ced: # new interval
            rows.append((cst,ced))
            cst,ced = st,ed
        else: # update end
            ced = ed
    rows.append((cst,ced))
    return rows


@cython.boundscheck(False) # turns off bounds-checking for entire function
@cython.wraparound(False) # turns off negative indexing checking
cpdef get_total_bp_bed12file_helper( bed12path):
    """ Returns total mapped base pairs (totbp) and covered base pairs (covbp).
    The ratio totbp/covbp gives average coverage. Process without reading entire data
    into the RAM.

    Args:
        bedpath: a path to BED12

    Returns:
        totbp: total base pairs in BED
        covbp: covered base pairs
    """
    cdef dict chrdic = {}
    cdef dict covbpdic = {}
    cdef dict totbpdic = {}
    cdef int st0,siz,st1,ed1
    cdef bytes line,chrom,st,ed,name,sc1,strand,tst,ted,sc2,nexons,esizes,estarts,y,z
    # using bytes rather than decoding is much faster
    cdef ADTYPE_t[I32_t,ndim=2] a
    cdef list b
    cdef set v

    if bed12path[len(bed12path)-3:len(bed12path)]=='.gz':
        fp = gzip.open(bed12path,'r')
    else:
        fp = open(bed12path,'r')
    for line in fp:#.readlines():
        chrom,st,ed,name,sc1,strand,tst,ted,sc2,nexons,esizes,estarts = line.strip().split(b'\t')
        if chrom not in chrdic:
            chrdic[chrom] = set()
            totbpdic[chrom] = 0
        st0 = int(st)
        if esizes[len(esizes)-1]==b',':
            esizes = esizes[:len(esizes)-1]
            estarts = estarts[:len(esizes)-1]
        for y,z in zip(esizes.split(b','),estarts.split(b',')):
            siz = int(y)
            st1 = st0+int(z)
            ed1 = st1+siz
            totbpdic[chrom] += siz
            chrdic[chrom].add((st1,ed1))
    fp.close()
    # union intervals
    for chrom, v in chrdic.items():
        a = N.array(sorted(list(v)), dtype=N.int32) 
        # numpy a.sort(axis=1) will sort two columns independently! Be careful.
        b = union_intervals(a) # union intervals
        covbpdic[chrom] = N.sum([ed1-st1 for st1,ed1 in b])
    return totbpdic, covbpdic

@cython.boundscheck(False) # turns off bounds-checking for entire function
@cython.wraparound(False) # turns off negative indexing checking
cpdef get_total_bp_bedfile_helper( bedpath):
    """ Returns total mapped base pairs (totbp) and covered base pairs (covbp).
    The ratio totbp/covbp gives average coverage. Process without reading entire data
    into the RAM. Process non BED12 file.

    Args:
        bedpath: a path to non BED12 (BED3,...)

    Returns:
        totbp: total base pairs in BED
        covbp: covered base pairs
    """
    cdef dict chrdic = {}
    cdef dict covbpdic = {}
    cdef dict totbpdic = {}
    # cdef set chrsted = set()
    cdef int st0,ed0
    #cdef str line,chrom,st,ed,y,z
    cdef bytes line,chrom,st,ed,y,z
    cdef set v
    cdef list b
    cdef ADTYPE_t[I32_t,ndim=2] a

    if bedpath[len(bedpath)-3:len(bedpath)]=='.gz':
        fp = gzip.open(bedpath,'r')
        #reader = codecs.getreader("utf-8")
        #fp = reader( gfp )
    else:
        fp = open(bedpath,'r')
        #gfp = fp
    for line in fp:#.readlines():
        chrom,st,ed = line.split(b'\t')[:3]
        if chrom not in chrdic:
            chrdic[chrom] = set()
            totbpdic[chrom] = 0
        st0,ed0 = int(st),int(ed)
        chrdic[chrom].add((st0,ed0))
        totbpdic[chrom] += ed0-st0
    #gfp.close()
    fp.close()
    # union intervals
    for chrom, v in chrdic.items():
        a = N.array(sorted(list(v)), dtype=N.int32)
        # numpy a.sort(axis=1) will sort two columns independently! Be careful.
        b = union_intervals(a) # union intervals
        covbpdic[chrom] = N.sum([ed0-st0 for st0,ed0 in b])
    return totbpdic, covbpdic

@cython.boundscheck(False) # turns off bounds-checking for entire function
@cython.wraparound(False) # turns off negative indexing checking
cpdef get_total_bp_bedfile_helper_check_uniq( bedpath):
    """ Returns total mapped base pairs (totbp) and covered base pairs (covbp).
    The ratio totbp/covbp gives average coverage. Process without reading entire data
    into the RAM. Process non BED12 file.

    Check the 4th column (name) as read id. Omits duplicates (but take the first one).

    Args:
        bedpath: a path to non BED12 (BED3,...)

    Returns:
        totbp: total base pairs in BED
        covbp: covered base pairs
    """
    cdef dict chrdic = {}
    cdef dict covbpdic = {}
    cdef dict totbpdic = {}
    # cdef set chrsted = set()
    cdef int st0,ed0,cid0,pid0,rid0, prid0
    #cdef str line,chrom,st,ed,y,z
    cdef bytes line,chrom,st,ed,y,z,cid,sc1,strand,rid
    cdef set v
    cdef list b
    cdef ADTYPE_t[I32_t,ndim=2] a

    if bedpath[len(bedpath)-3:len(bedpath)]=='.gz':
        fp = gzip.open(bedpath,'r')
        #reader = codecs.getreader("utf-8")
        #fp = reader( gfp )
    else:
        fp = open(bedpath,'r')
        #gfp = fp
    pid0 = 0 # previous read id 
    prid0 = 0 # previous map id
    for line in fp:#.readlines():
        chrom,st,ed,cid,sc1,strand,rid = line.split(b'\t')[:7]
        cid0,rid0 = int(cid),int(rid)
        if (cid0>pid0) or (rid0==prid0):
            if chrom not in chrdic:
                chrdic[chrom] = set()
                totbpdic[chrom] = 0
            st0,ed0 = int(st),int(ed)
            chrdic[chrom].add((st0,ed0))
            totbpdic[chrom] += ed0-st0
        pid0 = cid0
        prid0 = rid0
    #gfp.close()
    fp.close()
    # union intervals
    for chrom, v in chrdic.items():
        a = N.array(sorted(list(v)), dtype=N.int32)
        # numpy a.sort(axis=1) will sort two columns independently! Be careful.
        b = union_intervals(a) # union intervals
        covbpdic[chrom] = N.sum([ed0-st0 for st0,ed0 in b])
    return totbpdic, covbpdic

@cython.boundscheck(False) # turns off bounds-checking for entire function
@cython.wraparound(False) # turns off negative indexing checking
cpdef list flatten_bed8(ADTYPE_t bed8):
    cdef list rslt, x1
    cdef int i
    cdef str siz,sta,y,z
    cdef Py_ssize_t n = len(bed8)
    rslt = []
    for i in range(n):
        x = bed8[i]
        siz = str(x[6])
        sta = str(x[7])
        if siz[len(siz)-1]==',':
            siz = siz[:len(siz)-1]
            sta = sta[:len(siz)-1]
        for y,z in zip(siz.split(','),sta.split(',')):
            x1 = list(x)
            x1[6] = int(y)
            x1[7] = int(z)
            rslt.append(x1)
    return rslt


@cython.boundscheck(False) # turns off bounds-checking for entire function
@cython.wraparound(False) # turns off negative indexing checking
cpdef str array2wiggle_chr(N.ndarray[F32_t] a,  chrom,  dstpath):
    #cdef unicode txt
    cdef int i,j,st
    cdef F32_t c

    makedirs2(dstpath)
    fobj = open(dstpath,'w')
    i = 0
    cdef Py_ssize_t n = len(a)
    for i in range(n):
        if a[i]!=0:
            break
    st = i
    c = a[st] # initial non-zero
    for i in range(st+1,n):
        if a[i]==c: # same
            continue
        else:
            if c!=0:
                fobj.write(u'{0}\t{1}\t{2}\t{3}\n'.format(chrom,st,i,c))
            st = i
            c = a[st]
    if c!=0:
        fobj.write(u'{0}\t{1}\t{2}\t{3}\n'.format(chrom,st,i,c))
    fobj.close()
    return dstpath

@cython.boundscheck(False) # turns off bounds-checking for entire function
@cython.wraparound(False) # turns off negative indexing checking
cpdef str array2wiggle_chr64(N.ndarray[F64_t] a,  chrom,  dstpath):
    #cdef unicode txt
    cdef int i,j,st
    cdef F64_t c

    makedirs2(dstpath)
    fobj = open(dstpath,'w')
    i = 0
    cdef Py_ssize_t n = len(a)
    for i in range(n):
        if a[i]!=0:
            break
    st = i
    c = a[st] # initial non-zero
    for i in range(st+1,n):
        if a[i]==c: # same
            continue
        else:
            if c!=0:
                fobj.write(u'{0}\t{1}\t{2}\t{3}\n'.format(chrom,st,i,c))
            st = i
            c = a[st]
    if c!=0:
        fobj.write(u'{0}\t{1}\t{2}\t{3}\n'.format(chrom,st,i,c))
    fobj.close()
    return dstpath


@cython.boundscheck(False) # turns off bounds-checking for entire function
@cython.wraparound(False) # turns off negative indexing checking
cpdef read_gtf_helper( gtfpath, list parseattrs,  comment='#'):
    """ Returns total mapped base pairs (totbp) and covered base pairs (covbp).
    The ratio totbp/covbp gives average coverage. Process without reading entire data
    into the RAM.

    Args:
        bedpath: path to GTF
        parseattrs: attributes (in extra fields) to parse
        comment: comment string at the beginning of the lines

    Returns:
        row: list of parsed rows
        cols: column names

    """
    cdef str line,chrom,src,typ,st,ed,sc1,strand,sc2,extra,x,y
    cdef list l,r
    cdef list GTFCOLS = ['chr','src','typ','st','ed','sc1','strand','sc2','extra']
    cdef list recs = []
    cdef int st0, ed0
    cdef dict dic
    cdef list cols0 = GTFCOLS+parseattrs
    cdef list e = ['']*len(parseattrs)

    if gtfpath[len(gtfpath)-3:len(gtfpath)]=='.gz':
        gfp = gzip.open(gtfpath,'r')
        reader = codecs.getreader("utf-8")
        fp = reader( gfp )
    else:
        fp = open(gtfpath,'r')
        gfp = fp

    for line in fp: # skip initial comments
        if len(line)>0 and line[0]!=comment:
            break
    # process the first line
    r = line.strip().split('\t')
    if len(r)==9:
        chrom,src,typ,st,ed,sc1,strand,sc2,extra = r
        dic = dict([(l[0],l[1][1:-1]) for l in [y.split() for y in extra.split(';')] if len(l)>1])
        recs.append(r+[dic.get(x,'') for x in parseattrs])
    else:
        recs.append(r)        
    for line in fp: # process the rest
        r = line.strip().split('\t')
        if len(r)==9:
            chrom,src,typ,st,ed,sc1,strand,sc2,extra = r
            dic = dict([(l[0],l[1][1:-1]) for l in [y.split() for y in extra.split(';')] if len(l)>1])
            recs.append(r+[dic.get(x,'') for x in parseattrs])
        else:
            recs.append(r)
    gfp.close()

    return recs, cols0

