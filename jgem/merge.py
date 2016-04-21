"""Copyright (c) 2015-2016 Ken Sugino

.. module:: merge
    :synopsis: module for merging multiple assemblies

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""

import csv
import subprocess
import os
import gzip
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

import pandas as PD
import numpy as N

from jgem import utils as UT
