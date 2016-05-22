"""

.. module:: assembler
    :synopsis: assemble genes from RNASeq data (normalized genome coverage (bigwig) and junctions)

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
# system imports
import gzip
import os
from functools import reduce
from operator import iadd, iand
from collections import Counter
from itertools import repeat
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

# 3rd party imports
import pandas as PD
import numpy as N
import matplotlib.pyplot as P

# library imports
from jgem import utils as UT

