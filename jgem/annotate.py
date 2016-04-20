"""
.. module:: compare
    :synopsis: compare to a reference and annotate known genes

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
# system imports
import gzip
import os
import subprocess
from collections import Counter
from operator import iadd

# 3rd party libraries
import pandas as PD
import numpy as N
import matplotlib.pyplot as P

# library imports
from jgem import utils as UT
from jgem import gtfgffbed as GGB
from jgem import bedtools as BT
from jgem import bigwig as BW


