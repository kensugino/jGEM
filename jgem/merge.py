"""Copyright (c) 2015-2016 Ken Sugino

.. module:: merge
    :synopsis: module for merging multiple assemblies

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""

import subprocess
import os
import gzip
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

import pandas as PD
import numpy as N

from jgem import utils as UT
from jgem import gtfgffbed as GGB
from jgem import filenames as FN
from jgem import assembler as AS


class MergeNames(FN.FileNamesBase):
	"""Filelname manager for merging process.

    Attributes:
        sampleinfo: sample info dataframe (with columns: name, sjexbase, bwfile, sjfile)
        genome: UCSC genome name
        code: merge identifier
        outdir: output directory

    All outputs and temporary files are prefixed by **outdir/code**

	"""

    def __init__(self, sampleinfo, genome, code, outdir):
        self.si = sampleinfo
        self.genome = genome
        self.chromsizes = UT.chromsizes(genome)
        self.chroms = UT.chroms(genome)
        self.code = code
        self.outdir = outdir

        prefix = os.path.join(outdir, code)
        super(MergeNames, self).__init__(prefix)

    def 

