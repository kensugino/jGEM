"""

.. module:: fasta
    :synopsis: FASTA related functions, includes retrieving seqs from genome FASTA.

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""

import csv
import subprocess
import os
import gzip
import glob
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

import pandas as PD
import numpy as N


#### FASTA/PANDAS ######################################################
def fasta2panda(fname):
    if fname.endswith('.gz'):
        fa = gzip.open(fname).read()
    else:
        fa = open(fname).read()
    def _parse(x):
        lines = x.split('\n')
        tid = lines[0].split()[0]
        seq = ''.join(lines[1:])
        return tid, seq
    recs = [_parse(x) for x in fa.split('>') if x.strip()]
    fadf = PD.DataFrame(recs, columns=['tid','seq'])
    return fadf

#### GENOME seq  ######################################################

class GenomeFASTAfai(object):
    """ Random access to genome fast with .fai index """
    # .fai format: contig_name, contig_size, initial position in file in bytes, #bp in a line, #bytes in a line
    pass


class GenomeFASTAChroms(object):
    """ Access to genome fast with single chromosome single file format """


    def __init__(self, chromdir, linesep='\n'):
        self.chromdir = chromdir
        self.linesep = linesep
        files = glob.glob(os.path.join(chromdir, '*.fa'))
        self.gz = False
        self.ext = '.fa'
        if len(files)==0:
            files = glob.glob(os.path.join(chromdir, '*.fasta'))
            self.ext = '.fasta'
        if len(files)==0:
            files = glob.glob(os.path.join(chromdir, '*.fa.gz'))
            self.gz = True
            self.ext = '.fa.gz'
        if len(files)==0:
            files = glob.glob(os.path.join(chromdir, '*.fasta.gz'))
            self.gz = True
            self.ext = '.fasta.gz'
        self.chroms = {} # cache whole chrom seq
        self.chromosomes = [os.path.basename(x)[:-len(self.ext)] for x in files]

    def _path(self, chrom):
        return os.path.join(self.chromdir, chrom+self.ext)

    def _read(self, chrom):
        if self.gz:
            with gzip.open(self._path(chrom)) as fp: # bytes (py3) or string (py2)
                txt = fp.read().decode()
        else:
            with open(self._path(chrom),'r') as fp:
                txt = fp.read()
        lines = txt.split(self.linesep) # could be '\r' or '\r\n' or '\n'
        tid = lines[0].split()[0]
        seq = ''.join(lines[1:])
        # assert(tid==chrom)
        self.chroms[chrom] = seq

    def get(self, chrom, st, ed):
        """ zero based index returns chromseq[st:ed] """
        if chrom not in self.chroms:
            self._read(chrom)
        return self.chroms[chrom][st:ed]





