""" 
    Expression Dataset for analysis of matrix (RNASeq/microarray) data with annotations

"""

import pandas as PD
import numpy as N
from matplotlib import pylab as P
from collections import OrderedDict
from ast import literal_eval
# from ..plot.matrix import matshow_clustered

class ExpressionSet(object):

    def __init__(self, eData, gData=None, sData=None):
        """ 
        eData: expression data (gene x samples) header: MultiIndex (samplename, group)
        fData: gene annotation (gene x gene annotations)
        pData: sample annotation (sample x sample annotations)
        """
        self.eData = eData
        self.gData = gData
        self.sData = sData

    def read(self, eFile, gFile=None, sFile=None):
        pass

    def write(self, eFile, gFile=None, sFile=None):
        self.eData.to_csv(eFile, tupleize_cols=False, sep="\t")
        if gFile is not None:
            self.gData.to_csv(gFile, tupleize_cols=False, sep="\t")
        if sFile is not None:
            self.sData.to_csv(sFile, tupleize_cols=False, sep="\t")

    def find(self, field, pat):
        pass


def read_bioinfo3_data(fname):
    """ read bioinfo3.table.dataset type of data """
    fobj = open(fname)
    groups = OrderedDict()
    cnt = 0
    for line in fobj:
        cnt += 1
        if line[:2]=='#%':
            if line.startswith('#%groups:'):
                gname, members = line[len('#%groups:'):].split('=')
                gname = gname.strip()
                members = members.strip().split(',')
                groups[gname] = members
                datafields = line.strip().split('=')[1].strip().split(',')
            elif line.startswith('#%fields'):
                fields = line.strip().split('=')[1].strip().split(',')
        elif not line.strip():
            continue #  empty line
        else:
            break
    df = PD.read_table(fname, skiprows=cnt-1)
    f2g = {}
    for g,m in groups.items():
        for f in m: 
            f2g[f] = g
            
    df.columns = PD.MultiIndex.from_tuples([(x, f2g.get(x,'')) for x in df.columns], names=['samplename','group'])
    e = ExpressionSet(df)
    return e

def read_multiindex_data(fname, tupleize=True, index_names = ['samplename','group']):
    """ read dataset table with MultiIndex in the header """
    if not tupleize:
        df = PD.read_table(fname, header=range(len(index_names)), index_col=[0], tupleize_cols=False)
        e = ExpressionSet(df)
        return e
    df = PD.read_table(fname, index_col=0)
    df.columns = PD.MultiIndex.from_tuples(df.columns.map(literal_eval).tolist(), names=index_names)
    e = ExpressionSet(df)
    return e

def read_grouped_table(fname, groupfn=lambda x: '_'.join(x.split('_')[:-1])):
    """ Read dataset whose group is encoded in the colname. Column 0 is index. """
    df = PD.read_table(fname)
    f2g = {x:groupfn(x) for x in df.columns}
    df.columns = PD.MultiIndex.from_tuples([(x, f2g[x]) for x in df.columns], names=['samplename','group'])
    e = ExpressionSet(df)
    return e



def concatenate(dic):
    """ dic: dict of DataFrames
        merge all using index and outer join
    """
    keys = list(dic)
    d = dic[keys[0]].merge(dic[keys[1]], left_index=True, right_index=True, how='outer', suffixes=('.'+keys[0],'.'+keys[1]))
    for k in keys[2:]:
        d = d.merge(dic[k], left_index=True, right_index=True, how='outer', suffixes=('','.'+k))
    return d

def calc_mergesortkey(dic, pos_neg_flds):
    conc = concatenate(dic)
    selected = ~N.isnan(conc[pos_neg_flds])
    pos = conc[pos_neg_flds]>0
    neg = conc[pos_neg_flds]<=0
    num_pos = pos.sum(axis=1)
    num_neg = neg.sum(axis=1)
    pos_neg_mix = -1*(num_neg==0) + 1*(num_pos==0)  # pos(-1),  mix(0), neg(1)
    #num_hit = num_pos - num_neg
    num_hit = num_pos + num_neg
    n = len(pos_neg_flds)
    #position = (N.arange(1,n+1)*pos + N.arange(-1,-n-1,-1)*neg).sum(axis=1)
    position = (N.arange(1,n+1)*pos + N.arange(-n,0)*neg).sum(axis=1)
    strength = (conc[pos_neg_flds]*pos).sum(axis=1) + (conc[pos_neg_flds]*neg).sum(axis=1)
    #msk = PD.Series(list(zip(pos_neg_mix, num_hit, position, strength)), index=conc.index)
    #msk.sort()
    conc['mergesortkey'] = list(zip(pos_neg_mix, num_hit, position, strength))
    conc.sort('mergesortkey', inplace=True)
    return conc



