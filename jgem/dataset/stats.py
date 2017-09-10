"""
 stats on datasets

 functions acts on pandas.DataFrame (first argument)

"""
import numpy as N
import pandas as PD

from scipy.stats import chi2
from .p_adjust import r_p_adjust

def select_by_lfc_and_q(df, lfc_th=N.log2(1.5), q_th=0.05, 
	lfc_prefix='Coef.', pval_prefix='cp.', conditions=None, method="BH"):
	conditions = _get_conditions(df, lfc_prefix, conditions)
	sel = {c: df[[lfc_prefix+c,pval_prefix+c]][df[lfc_prefix+c].abs()>=lfc_th] for c in conditions}
	for c in conditions:
		tmp = sel[c]
		pcol = pval_prefix+c
		acol = 'adj.'+pcol
		tmp[acol] = r_p_adjust(tmp[pcol].values, method)
		sel[c] = tmp[tmp[acol]<q_th]
	return sel

def adjust_p(df, pval_prefix='cp.', conditions=None, method="BH"):
	conditions = _get_conditions(df, pval_prefix, conditions)
	pvalcols = [pval_prefix+c for c in conditions]
	cols = ['adj.'+c for c in pvalcols]
	adjpvals = {'adj.'+c:r_p_adjust(df[c].values, method) for c in pvalcols}
	tmp = PD.DataFrame(adjpvals, index=df.index, columns=cols)
	return df.merge(tmp, left_index=True, right_index=True)

def removedup(df, removedup_col='symbol', pval_prefix='p.value.', lfc_prefix='Coef.', conditions=None, twotailed=False):
	conditions = _get_conditions(df, pval_prefix, conditions)
	cp = agg_combined_p(df, removedup_col, pval_prefix, conditions, twotailed)
	lfc = agg_op(df,'absmax', removedup_col, lfc_prefix, conditions)
	# find id to use for plotting (the one with minp)
	pcols = [pval_prefix+c for c in conditions]
	fcols = ['id.'+c for c in conditions]
	minpid = df.groupby(removedup_col).apply(lambda x:PD.Series(x.index[x[pcols].values.argmin(axis=0)], index=fcols))
	ids = agg_index(df, removedup_col)
	siz = agg_size(df, removedup_col)
	d = PD.DataFrame({'ids':ids, 'size':siz})
	d = d.merge(cp, left_index=True, right_index=True)
	d = d.merge(lfc, left_index=True, right_index=True)
	d = d.merge(minpid, left_index=True, right_index=True)
	return d

def removedup0(df, removedup_col='symbol', pval_prefix='p.value.', lfc_prefix='Coef.', conditions=None, twotailed=False):
	conditions = _get_conditions(df, pval_prefix, conditions)
	cp = agg_combined_p(df, removedup_col, pval_prefix, conditions, twotailed)
	minp = agg_op(df,'min', removedup_col, pval_prefix, conditions)
	#avglfc = agg_op(df,'mean', removedup_col, lfc_prefix, conditions)
	avglfc = agg_op(df,'absmax', removedup_col, lfc_prefix, conditions)
	pcols = [pval_prefix+c for c in conditions]
	fcols = [lfc_prefix+c for c in conditions]
	rnf = range(len(fcols))
	minplfc = df.groupby(removedup_col).apply(lambda x:PD.Series(x[fcols].values[x[pcols].values.argmin(axis=0),rnf], index=fcols))
	cp_minp = cp.values > minp.values
	cp1val = cp.values*(~cp_minp) + minp.values*cp_minp
	cp1 = PD.DataFrame(cp1val, index=cp.index, columns=cp.columns)
	lfcval = avglfc.values*(~cp_minp) + minplfc.values*cp_minp
	lfc1 = PD.DataFrame(lfcval, index=avglfc.index, columns=avglfc.columns)
	ids = agg_index(df, removedup_col)
	siz = agg_size(df, removedup_col)
	d = PD.DataFrame({'ids':ids, 'size':siz})
	d = d.merge(cp1, left_index=True, right_index=True)
	d = d.merge(lfc1, left_index=True, right_index=True)
	return d


def _get_conditions(df, prefix, conditions):
	if conditions is not None:
		return conditions
	lp = len(prefix)
	return [x[lp:] for x in df.columns if x.startswith(prefix)]

def agg_combined_p(df, groupby=None, prefix='p.value.', conditions=None, twotailed=True):
	""" Calculate Fisher's combined P value for the column specified by **pvalcol** for each group, 
		defined by **groupby** column. If groupby is None, use index for grouping. 

		Returns single col (pandas.Series) which contains combined p values. 
	"""
	conditions = _get_conditions(df, prefix, conditions)
	pvalcols = [prefix + c for c in conditions]
	cpcols = ['cp.'+c for c in conditions]
	pchisq = chi2.cdf
	# first calculate -2*log(p) for each group
	if twotailed:
		logpval = -2*N.log(df[pvalcols]/2.)
	else:
		logpval = -2*N.log(df[pvalcols])
	if groupby is None:
		gb = logpval.groupby(level=0)
	else:
		if type(groupby)!=list:
			logpval[groupby] = df[groupby]
		else:
			for col in groupby:
				logpval[col] = df[col]
		gb = logpval.groupby(groupby)
	all_sump = gb.sum()
	all_degf = gb.size()
	cp = 1. - pchisq(all_sump, 2*all_degf[:,N.newaxis])
	return PD.DataFrame(cp, index=all_sump.index, columns=cpcols)

def agg_op(df, op='median', groupby=None, prefix='Coef.', conditions=None):
	conditions = _get_conditions(df, prefix, conditions)
	cols = [prefix+c for c in conditions]
	if groupby is None:
		gb = df.groupby(level=0)
	else:
		gb = df.groupby(groupby)
	if op=='absmax':
		df_mean = gb[cols].agg(lambda x: x[N.abs(x).argmax()])
	else:
		df_mean = getattr(gb[cols], op)()
	#df_mean.rename(columns = {prefix+c:'mean.'+c for c in conditions}, inplace=True)
	return df_mean


def agg_index(df, groupby):
	return df.groupby(groupby).apply(lambda x: ','.join(x.index))

def agg_size(df, groupby):
	return df.groupby(groupby).size()

 
def detect_noise(df, level='group', cv_th=0.9, mean_th=3):

	avg = df.groupby(level=level,axis=1).mean()
	std = df.groupby(level=level,axis=1).std()
	cv = std/avg
	cv[N.isnan(cv)] = N.inf
	accept = (cv<=cv_th)*(avg>mean_th)
	idx =  accept.apply(lambda x: any(x), axis=1)
	return idx
