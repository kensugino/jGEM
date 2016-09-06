"""

.. module:: as2helper.pyx
    :synopsis: cython speedup of assembler2 related stuffs

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""

import numpy as N
cimport numpy as N
import gzip
import codecs
import os
import errno
from collections import defaultdict

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
cpdef I32_t find_maxgap(N.ndarray[float] arr, covfactor):
    cdef float emax = arr.max()
    cdef float th = emax*covfactor
    cdef float emin = arr.min()
    cdef int cmax = 1
    cdef int cst,ced,i

    if (emin>th):
        return 0
    idx = N.nonzero(arr<=th)[0]
    if len(idx)==0:
        return 0
    cst = idx[0]
    ced = cst
    for i in idx[1:]:
        if i==ced+1: # continuous
            ced = i
        else:
            cmax = max(cmax, ced-cst+1)
            cst = ced = i
    cmax = max(cmax, ced-cst+1)
    return cmax  

