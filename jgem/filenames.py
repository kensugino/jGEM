"""

.. module:: filenlames
    :synopsis: deals with (mostly temporary) files 

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""

import os
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)
import glob

# 3rd party imports
import numpy as N

# library imports
from jgem import gtfgffbed as GGB
from jgem import utils as UT
from jgem import convert as CV

class FileNamesBase(object):
    """Base calss for handling tons of (temporary) filenames.
    * All filenames have a common prefix. 
    * Each filename belong to a category.

    """
    def __init__(self, prefix):
        """
        Args:
            prefix (str): path prefix

        """
        self._prefix = prefix
        self._fnames = {}  #set()

    def fname(self, suffix, category='temp'):
        """Generate a filename (path).

        Args:
            suffix: (str)
            category: (str) default 'temp'

        Returns:
            (prefix).(suffix)

        """
        fn = '{0}.{1}'.format(self._prefix, suffix)
        self._fnames.setdefault(category, set()).add(fn)
        return fn

    def txtname(self, suffix, category='temp'):
        return self.fname(suffix+'.txt.gz', category)

    def bedname(self, suffix, category='temp'):
        return self.fname(suffix+'.bed.gz', category)

    def bedname2(self, suf, th, cat='temp'):
        return self.bedname('{0}{1:g}'.format(suf, th), cat)

    def delete(self, delete=['temp'], protect=[]):
        """Delete temporary files.

        Args:
            delete (list): categories to delete default ['temp']
            protect (list): output categories (files in these cateogries 
             are protpected)

        """
        outputs = []
        for c in protect:
            if c in self._fnames:
                outputs += self._fnames[c]
        if len(delete)==0:
            delete = [x for x in self._fnames.keys() if x not in outputs]
        for c in delete:
            if c in self._fnames:
                for fpath in self._fnames[c]:
                    if os.path.exists(fpath) and fpath not in outputs:
                        LOG.debug('deleting {0}...'.format(fpath))
                        os.unlink(fpath)                

    def delete_prefixed(self, suffix=''):
        pat = fn = self._prefix+'.{0}*'.format(suffix)
        flist = glob.glob(pat)
        for f in flist:
            LOG.debug('deleting {0}...'.format(f))
            try:
                os.unlink(f)
            except:
                pass

    def write_bed(self, df, suffix, ncols, category='temp', **kw):
        if ncols is None:
            ncols = len(df.columns)
        fname = self.bedname(suffix, category)
        return GGB.write_bed(df, fname, ncols=ncols, **kw)

    def write_txt(self, df, suffix, fm='h', category='temp', **kw):
        fname = self.txtname(suffix, category)
        return UT.write_pandas(df, fname, fm=fm, **kw)

    def read_bed(self, suffix, category='read'):
        return GGB.read_bed(self.bedname(suffix, category))

    def read_txt(self, suffix, category='read'):
        return UT.read_pandas(self.txtname(suffix, category))


class FileNames(FileNamesBase):
    """Filenames for the Assembler.
    
    Holds
    - sample name (sname)
    - bigwig coverage file path (bwfile)
    - junction file path (sjfile)
    - output directory (outdir)
    - reference (refgtf) if using it for finding SE (single exon) coverage threshold

    All filenames have **outdir+sname** as prefix.

    """
    def __init__(self,sname,bwfile,sjfile,outdir,refgtf='.gtf'):
        """Handles filenames.

        Args:
            sname (str): sample name
            bwfile (str): path to normalized coverage bigwig file
            sjfile (str): path to junction file
            outdir (str): path to output directory
            refgtf (str): (optional) path to reference GTF

        """
        self.sname = sname
        self.bwfile = bwfile
        self.sjfile = sjfile
        self.outdir = outdir
        self.refgtf = CV.GTF2SJEX(refgtf)
        
        # FileNamesBase init
        prefix = os.path.join(outdir, sname)
        super(FileNames, self).__init__(prefix)

    def refname(self,suf):
        return self.refgtf.fname(suf+'.txt.gz')

    def refsjex(self):
        return self.refgtf.sjex()

    def ex_out(self):
        return self.txtname('ex', category='output')

    def sj_out(self):
        return self.txtname('sj', category='output')

    def genes_out(self):
        return self.bedname('genes', category='output')

