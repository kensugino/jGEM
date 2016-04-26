"""

.. module:: trimex
    :synopsis: makes trimed genes (to avoid length bias in coverage variance)

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""
import multiprocessing
import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

# 3rd party imports
import pandas as PD
import numpy as N
# library imports
from jgem import utils as UT
from jgem import bigwig as BW


# Trim from 3' end    
def trim_ex_worker(args):
    ex,length,gidfld=args
    # for each gene collect exons from 3' end up to length
    cols = list(ex.columns.values)
    idxst = cols.index('st')
    idxed = cols.index('ed')
    idxlen = cols.index('len')
    def _gen():
        for _gidx, gr in ex.groupby(gidfld):
            total = gr['len'].sum()
            if total<length: # too short to trim
                for v in gr.values:
                    yield tuple(v)
            else:
                # look from 3' end
                strand = gr['strand'].iloc[0]
                if strand == '+':
                    gr = gr.sort_values(['ed','st'],ascending=False)
                else:
                    gr = gr.sort_values(['st','ed'],ascending=True)
                grv = gr.values
                # if the first one already larger than length?
                first = grv[0]
                clen = first[idxlen]
                if clen>length: # trim this one and go next
                    if strand == '+':
                        first[idxed] = first[idxst] + length
                        yield tuple(first)
                    else:
                        first[idxst] = first[idxed] - length
                        yield tuple(first)
                # start collecting until it reaches the length
                else:
                    yield first
                    for i in range(1,len(gr)):
                        cur = grv[i]
                        clen += cur[idxlen]
                        if clen<=length: # not enough yet
                            yield tuple(cur)
                            if clen==length:
                                break
                        else: # exceeded => trim
                            if strand=='+':
                                cur[idxst] = cur[idxst] + (clen-length)
                            else:
                                cur[idxed] = cur[idxed] - (clen-length)
                            yield tuple(cur)
                            break
    #nex = PD.DataFrame([x for x in _gen()], columns = cols)
    return [x for x in _gen()]

def trim_ex(expath, dstpath, dstcipath, length=1000, gidfld='_gidx', np=7):
    """Generate trimmed version of genes for calculating coverage to avoid length bias. 

    Args:
        expath (str): path exon tsv
        dstpath (str): path to trimmed exon
        dstcipath (str): path to ci (chopped interval) 
        length (pos int): length to trim from 3' end in base pair (default 1000 bp)
        gidfld (str): column name for gene id (default _gidx)
        np (pos int): number of CPU to use

    Generates:
        Two files (dstpath, dstcipath).

    Returns:
        a dataframe containing trimmed exons
    """
    #ex = UT.read_pandas(MD.paths[code]['ex'])
    #dstpath = MD.trimmedex[code][length]['ex']
    #dstcipath = MD.trimmedex[code][length]['ci']
    ex = UT.read_pandas(expath)
    if 'len' not in ex.columns:
        ex['len'] = ex['ed'] - ex['st']
    if np==1:
        recs = trim_ex_worker((ex, length, gidfld))
    else:
        chroms = sorted(ex['chr'].unique())
        data = [(ex[ex['chr']==c], length, gidfld) for c in chroms]
        recs = []
        try:
            p = multiprocessing.Pool(np)
            for v in p.map(trim_ex_worker, data):
                recs += v
            #recs = reduce(iadd, p.map(trim_ex_worker, *zip(*data)))
        finally:
            p.close()
            # p.join()
    cols = list(ex.columns.values)
    nex = PD.DataFrame(recs, columns = cols)
    nex['len'] = nex['ed'] - nex['st']
    # edge case
    nex.loc[nex['st']==nex['ed'],'ed'] = nex['st'] + 1
    UT.save_tsv_nidx_whead(nex, dstpath)
    UT.chopintervals(nex, dstcipath)

    return nex
