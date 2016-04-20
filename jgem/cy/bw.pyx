import numpy as N
cimport numpy as N

DTYPE = N.float32
ctypedef N.float32_t DTYPE_t
cimport cython

@cython.boundscheck(False) # turns off bounds-checking for entire function
@cython.wraparound(False) # turns off negative indexing checking
cpdef str array2wiggle_chr(N.ndarray[DTYPE_t] a, str chrom, str dstpath):
    #cdef unicode txt
    cdef int i,j,st
    cdef DTYPE_t c

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
