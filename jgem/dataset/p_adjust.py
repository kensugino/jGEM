from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import FloatVector
#from numba import autojit, jit

stats = importr('stats')

def r_p_adjust(pvals, method):
	return N.array(stats.p_adjust(FloatVector(pvals), method = method))

import numpy as N
#import pandas as PD

METHODS = ["holm","hochberg","hommel","bonferroni","BH","BY","fdr","none" ]

#@autojit
def pmin(*p):
	return N.minimum(*p)

def pmax(*p):
	return N.maximum(*p)

def order(p, decreasing=False):
	tmp = [(x,i) for i,x in enumerate(p)]
	tmp.sort(reverse=decreasing)
	return [x[1] for x in tmp]

def cummax(p):
	return N.maximum.accumulate(p)

def cummin(p):
	return N.minimum.accumulate(p)

def seq_len(x):
	return N.arange(1,x+1)

def seq(a,b):
	if a<b:
		return N.arange(a,b+1)
	return N.arange(a,b-1,-1)

#@autojit
def p_adjust(pvals, method = "BH", n=None):
	"""
	from R p.adjust

	p: array of pvalues
	method: one of "holm","hochberg","hommel","bonferroni","BH","BY","fdr","none" 
	n: number of comparison

	"""
	assert method in METHODS
	if method == "fdr":
		method = "BH"
	#p0 = N.array(pvals)
	p0 = pvals.copy()
	nna = N.nonzero(~N.isnan(p0))
	p = p0[nna]
	if n is None:
		n = len(p)
	lp = len(p)

	if n <= 1: # no multiple comparison
		return p0
	if n == 2 and method == "hommel":
		method = "hochberg"

	#@autojit
	def bonferroni():
		return pmin(1., n * p)

	def holm():
		i = seq_len(lp)
		o = order(p)
		ro = order(o)
		return pmin(1., cummax((n - i + 1) * p[o]))[ro]

	def hommel():
		if n > lp:
			p = N.concatenate(p, N.ones(n - lp))
		i = seq_len(n)
		o = order(p)
		p = p[o]
		ro = order(o)
		q = N.ones(n) * N.min(float(n) * p/i)
		pa = N.ones(n) * N.min(float(n) * p/i)
		for j in range((n - 1), 1, -1):
			ij = seq_len(n - j + 1)
			i2 = seq(n - j + 2, n)
			q1 = N.min(j * p[i2]/N.arange(2,j+1))
			q[ij] = pmin(j * p[ij], q1)
			q[i2] = q[n - j + 1]
			pa = pmax(pa, q)
		i = ro[:lp] if lp < n else ro
		return pmax(pa, p)[i]

	def hochberg():
		i = seq(lp,1)
		o = order(p, decreasing = True)
		ro = order(o)
		return pmin(1., cummin((n - i + 1) * p[o]))[ro]

	def BH():
		i = seq(lp,1)
		o = order(p, decreasing = True)
		# tmp = [(x,i) for i,x in enumerate(p)]
		# tmp.sort(reverse=True)
		# o = [x[1] for x in tmp]
		ro = order(o)
		# tmp = [(x,i) for i,x in enumerate(o)]
		# tmp.sort(reverse=False)
		# ro = [x[1] for x in tmp]
		return pmin(1., cummin(float(n)/i * p[o]))[ro]

	def BY():
		i = seq(lp,1)
		o = order(p, decreasing = True)
		ro = order(o)
		q = N.sum(1./seq(1,n))
		return pmin(1., cummin(q * n/i * p[o]))[ro]
	
	def none():
		return p

	p0[nna] = locals()[method]()
	return p0








