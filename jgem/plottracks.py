"""Basic parts for plotting bigwig (coverage etc.), genes, ideograms. 
"""
import os
import re
try:
	from itertools import izip
except:
	izip = zip
from itertools import chain

import numpy as N
import pandas as PD

import matplotlib.pyplot as PP
from matplotlib.collections import BrokenBarHCollection
import matplotlib.patches as MP

#from ngslib import wWigIO
#import bx
#from bx.bbi.bigwig_file import BigWigFile
from jgem.bxbbi.bigwig_file import BigWigFile

import jgem.plotutils as PU
import jgem.plotgenes as PG

FONTSIZE=6

####################### util for drawing Ideograms
# modified based on https://www.biostars.org/p/147364/#147637
IDEOCOLS = ['chrom', 'start', 'end', 'name', 'gieStain']

class Ideograms(object):
	
	color_lookup = {
		'gneg': (1., 1., 1.),
		'gpos25': (.6, .6, .6),
		'gpos33': (.5,.5,.5),
		'gpos50': (.4, .4, .4),
		'gpos66': (.3,.3,.3),
		'gpos75': (.2, .2, .2),
		'gpos100': (0., 0., 0.),
		'acen': (.8, .4, .4),
		'gvar': (.8, .8, .8),
		'stalk': (.9, .9, .9),
	}

	chrom_list = ['chr%s' % i for i in (list(range(1, 20)) + ['X', 'Y'])]
	
	def __init__(self, fpath, chrom_list=None):
		if chrom_list is not None:
			self.chrom_list = chrom_list
		self.ideo = PD.read_table(fpath,skiprows=1,names=IDEOCOLS)
		self.ideo = self.ideo[[x in self.chrom_list for x in self.ideo.chrom]]
		self.ideo['width'] = self.ideo.end - self.ideo.start
		self.ideo['colors'] = [self.color_lookup[x] for x in self.ideo.gieStain]

	def make_collection(self, df, chrom, ypos=0.25, height=0.5, **kwargs):
		if 'width' not in df.columns:
			df['width'] = df.end - df.start
		df0 = df[df['chrom']==chrom]
		xranges = df0[['start','width']].values
		colors = df0['colors']#.values
		yrange = (ypos,height)
		bbhc = BrokenBarHCollection(xranges, yrange, facecolors=colors, **kwargs)
		xmin = df0.start.min()
		xmax = df0.end.max()
		return bbhc, xmin, xmax
	
	def draw_one(self, chrom, ax, regions=[], drawname=False, **kwargs):
		h = 0.4
		y = 0.2
		bbhc,xmin,xmax = self.make_collection(self.ideo, chrom, y, h, **kwargs)
		ax.add_collection(bbhc)
		ax.set_xlim(xmin-1e7,xmax+1e7)
		y1 = 0.05
		h1 = y+h+(y-y1)
		for region in regions:
			x = region[0]
			w = max(region[1] - x, xmax*0.02)
			ax.add_patch(MP.Rectangle((x,y1),w,h1, facecolor='red', linewidth=0, alpha=0.8))
		if drawname:
			xmid = (xmin+xmax)/2.
			ypos = y1+h1
			ax.text(xmid,ypos,chrom,ha='center',va='bottom')
		PP.setp(ax, xticks=[], yticks=[], frame_on=False)
		return ax, xmin, xmax, y1+h1
		
class Panel(object):
	"""
	Holds tracks.
	"""
	def __init__(self, tracks=[], figsize=(3,4)):
		self.tracks = tracks
		self.figsize = figsize
		
	def draw(self, ax=None, frameon=False):
		if ax is None:
			fig,ax = PP.subplots(1,1,figsize=self.figsize)
			PP.setp(ax, xticks=[], yticks=[], frame_on=False)
		# make axes for each track and dispatch draw command
		subaxes = self.make_subaxes(ax, frameon)
		for sax, track in zip(subaxes, self.tracks):
			PP.setp(sax, xticks=[], yticks=[])
			track.draw(ax=sax)
		self.ax = ax
		return ax
		
	def make_subaxes(self, ax, frameon=True):
		fig = ax.get_figure()
		tracks = self.tracks
		n = len(tracks)
		# simplest for now: divide equally in vertical direction
		tot_h = float(N.sum([x.h for x in tracks]))
		self.hs = hs = [x.h/tot_h for x in tracks]
		self.ys = ys = list(N.cumsum(hs[::-1])[::-1])[1:]+[0.]
		rects = [[0.,y,1.,0.9*h] for h,y in zip(hs,ys)]
		#h = 1./n
		#rects = [[0.,h*i,1.,0.9*h] for i in range(n)][::-1]
		args = dict(axisbg='w',frameon=frameon, xticks=[], yticks=[])
		subaxes = [PU.add_subplot_axes(ax,r,args,noticks=True) for r in rects]
		return subaxes
	
class Track(object):
	
	def __init__(self, name):
		self.name = name
		self.h = 1.
	
	def draw(self, ax):
		ax.text(0.5,0.5,self.name,ha='center',va='center')
	 
def locus2pos(locus):
	chrom, tmp = locus.split(':')
	st, ed = map(int, tmp.split('-'))
	return chrom,st,ed
def pos2locus(pos):
	return '{0}:{1:d}-{2:d}'.format(*pos)

class Ideogram(Track):
	
	def __init__(self, pos, ideofile, h=0.7, fontsize=FONTSIZE):
		"""
		Args:
			pos: (chrom, st, ed)
			ideofile: path to ideogram file
			fontsize: font size
		"""
		self.h = h
		chrom, st, ed = pos
		self.chrom = chrom
		self.pos = (st,ed)
		self.name = pos2locus(pos)
		self.ideo = Ideograms(fpath=ideofile)
		self.fontsize = fontsize
	
	def draw(self, ax):
		ax,xmin,xmax,ymax = self.ideo.draw_one(self.chrom, ax, [self.pos])
		xmid = (xmin+xmax)/2.
		ax.text(xmid,ymax,self.name,ha='center',va='bottom',fontsize=self.fontsize)
	
class Bed12Gene(Track):
	
	def __init__(self, pos, bedline=None, attr=None, color='k', **kwargs):
		if attr is not None:
			self.attr = attr
		else:
			assert(bedline is not None)
			self.attr = dict(zip(BEDCOLS, bedline.split('\t')))
			for c in ['esizes','estarts']:
				self.attr[c] = map(int, self.attr[c].strip()[:-1].split(','))
			for c in ['st','ed','tst','ted','nexons']:
				self.attr[c] = int(self.attr[c])
		self.name = self.attr['name']
		
		self.chrom, self.st, self.ed = pos
		assert (self.st<self.ed)

		self.color=color
		self.kwargs = kwargs
		
	def draw(self, ax):
		a = self.attr
		st,ed = a['st'],a['ed']
		xmin = self.st #- width*self.margin
		xmax = self.ed #+ width*self.margin
		width = xmax-xmin #ed-st
		#print 'gene,st,ed,width,margin,xmin,xmax',st,ed,width,self.margin,xmin,xmax
		xmid = (xmin+xmax)/2.
		ymax=0.8
		h = 0.4
		ymin = ymax - h
		ymid = ymax-h/2.
		# draw genome line
		ax.plot([xmin,xmax], [ymid,ymid], 'grey') # base
		# draw arrows in introns
		estsi = zip(a['estarts'],a['esizes'])
		ists = [st+est+esi for est,esi in estsi[:-1]]
		ieds = [st+est for est,esi in estsi[1:]]
		if a['strand']=='+':
			dx,shape,hw = 1,'right',0.25
		else:
			dx,shape,hw = -1,'right',0.4
		hl = 0.05*width
		arrowargs = dict(y=ymid, dx=dx, dy=0, shape=shape,
						 fc='grey', linewidth=0,
						 head_width=hw, head_length=hl)
		#print estsi
		#print zip(ists,ieds)
		for ist, ied in zip(ists, ieds):
			iwidth = ied-ist
			imid = (ist+ied)/2.
			if iwidth<0.15*width:
				continue
			#print imid
			ax.arrow(imid-dx*(hl+1)/2.,**arrowargs)
		# put arrow at xmin and xmax
		if st>xmin:
			if dx<0:
				axmin = xmin+hl-dx
				axmax = xmax-dx
			else:
				axmin = xmin-dx
				axmax = xmax-hl-dx
			ax.arrow(axmin, **arrowargs)
			ax.arrow(axmax, **arrowargs)
		# draw exons => BrokenBarHCollection
		# draw st=>TSS & TSE=>ed
		tss = a['tst'] - st # exon coords are st based 
		tse = a['ted'] - st
		if tss==tse:
			estsi1 = estsi3 = []
			estsi2 = estsi
		else:
			estsi1 = [(x,y) for x,y in estsi if x<tss]
			estsi2 = [(x,y) for x,y in estsi if ((x+y)>=tss)&(x<=tse)]
			estsi3 = [(x,y) for x,y in estsi if (x+y)>tse]
			if estsi1:
				x0,y0 = estsi1[-1] # last size needs fixing
				if (x0+y0)>tss:
					estsi1[-1] = (x0, tss-x0)
			if estsi2:
				x0,y0 = estsi2[0]
				if x0<tss: # first start and size need fixing
					estsi2[0] = (tss, y0-(tss-x0))
				x0,y0 = estsi2[-1]
				if (x0+y0)>tse:
					estsi2[-1] = (x0,tse-x0)
			if estsi3:
				x0,y0 = estsi3[0]
				if x0<tse:
					estsi3[0] = (tse,y0-(tse-x0))
		c = self.color
		cargs = dict(facecolor=c, edgecolor=c)#, alpha=0.8)
		cargs.update(self.kwargs)
		#print 'estsi1',estsi1
		#print 'estsi2',estsi2
		#print 'estsi3',estsi3
		# draw UTR
		if (len(estsi1)+len(estsi3))>0:
			yrange = (ymid-h/4., h/2.)
			xranges = [(st+x,y) for x,y in estsi1]+[(st+x,y) for x,y in estsi3]
			bbhc = BrokenBarHCollection(xranges, yrange, **cargs)
			ax.add_collection(bbhc)        
			#print '1,3 xranges', xranges
			#print '1,3 yrange', yrange
		# draw coding
		yrange = (ymid-h/2., h)
		xranges = [(st+x,y) for x,y in estsi2]
		#print '2 xranges', xranges
		#print '2 yrange', yrange
		bbhc = BrokenBarHCollection(xranges, yrange, **cargs)
		ax.add_collection(bbhc)

		# draw gene name
		txt = '%s (%.1fkb)' % (a['name'], width/1000)
		ax.text(xmid, 0, txt, ha='center', va='bottom', fontsize=FONTSIZE)
		PP.setp(ax, xticks=[], yticks=[], frame_on=False, 
				xlim=(xmin,xmax), ylim=(0,1))
		#print "gene", txt, xmin, xmax
		

def compress(wigs, resolution, th=1000):
	if len(wigs)<th:
		w = wigs
	else:
		def _gen():
			st,ed,h = wigs[0]
			v = h*(ed-st)
			for x1,x2,h in wigs:
				v = (v+h*(x2-ed))
				if (x2-st)>resolution: 
					h = v/(x2-st)
					yield (st,x2,h)
					st,ed,v = x2,x2,0
				else:
					ed = x2
		w = [x for x in _gen()]
	return wig2xy(w)

def wig2xy(wigs):
	x1,x2,h = zip(*wigs)
	z = N.zeros(len(h))
	x = list(chain.from_iterable(izip(x1,x1,x2,x2)))
	y = list(chain.from_iterable(izip( z, h, h, z)))
	return x,y

def compress2(wigs, window, minbins):
	st,x2,h = wigs[0]
	x1,ed,h = wigs[-1]
	nbins = (ed-st)/window
	if nbins<minbins:
		window = max(2, int(float(ed-st)/minbins))
	nbins = (ed-st)/window
	if nbins >= 4*len(wigs):
		return wig2xy(wigs)
	y0 = N.zeros(ed-st)
	for x1,x2,h in wigs:
		y0[x1-st:x2-st] = h
	win = N.ones(window)
	y = subsample(y0, window)
	x = N.arange(st, ed+window, window)[:len(y)]
	return x,y

def subsample(arr, n):
    end =  n * int(len(arr)/n)
    return N.mean(arr[:end].reshape(-1, n), 1)

class BigWig(Track):
	
	def __init__(self, fname, pos, name=None, h=0.6,
				 ymax=5,drawymax=False,drawname=True,
				 color='b',fontsize=FONTSIZE,resolution=100,minbins=100):
		self.h = h
		self.fname = fname
		self.pos = pos
		if name is None:
			name = os.path.basename(fname).replace('_star','').replace('.bw','')
		self.name = name
		self.chrom,self.st, self.ed = pos
		assert (self.st<self.ed)
		#self.margin = margin
		width = self.ed - self.st
		self.xmin = self.st# - margin*width
		self.xmax = self.ed# + margin*width
		self.ymax = ymax
		self.color = color
		self.drawymax=drawymax
		self.drawname=drawname
		self.fontsize = fontsize
		self.resolution=resolution
		self.minbins=minbins
		#print 'st,ed,width,margin,xmin,max', self.st,self.ed,width,margin, self.xmin, self.xmax
		#self.bw = BigWigFile(fname)
		
	def draw(self, ax):
		#bw = self.fname
		#wWigIO.open(bw)
		#wigs = wWigIO.getIntervals(bw, self.chrom,self.st,self.ed)
		#wWigIO.close(bw)
		with open(self.fname) as fobj:
			bw = BigWigFile(fobj)
			wigs = bw.get(self.chrom, self.st, self.ed)
		xmin,xmax = self.xmin,self.xmax
		ymax = self.ymax
		fs = self.fontsize
		# draw track name
		if self.drawname:
			ax.text((xmin+xmax)/2., ymax, self.name, ha='center',va='top',fontsize=fs)
		# draw ymax indicator
		if self.drawymax:
			ax.text(xmax,ymax,'%g'% ymax, ha='right',va='top',fontsize=fs)
		if len(wigs)>0:
			#x,y = compress(wigs, self.resolution)
			x,y = compress2(wigs,self.resolution,self.minbins)
			#print "%s #wigs %d #x/4 %d" % (self.name, len(wigs), len(x)/4)
			ax.fill_between(x,0,y,linewidth=0,color=self.color)
		# plot baseline
		ymin = -(ymax/30.)
		ax.plot([xmin,xmax],[ymin,ymin],'grey')
		PP.setp(ax, xlim=(xmin,xmax), ylim=(2*ymin,ymax), frame_on=False)


class Gene(Track):
	"""Track for plotting a gene using :class:SpliceFig class in :module:plotgene
	"""
	def __init__(self, pos, ex, sj, cmap='R', collapsey=False, h=1,
				compress=False, fontsize=FONTSIZE, **kwargs):
		if collapsey:
			self.h = h
		else:
			self.h = h*1.5
		self.sf = PG.SpliceFig(ex,sj,compress=compress,fontsize=fontsize,**kwargs)
		self.compress = compress
		self.chrom,self.st, self.ed = pos
		self.xlim = (self.st,self.ed)
		self.cmap = cmap
		self.collapsey = collapsey
		
	def draw(self, ax):
		sf = self.sf
		if self.collapsey:
			sf.ex['ey'] = 0
		if self.compress:
			sf.draw(ax, xlim=None, cm=self.cmap)
		else:
			sf.draw(ax, xlim=self.xlim, cm=self.cmap)

class Line(Track):
	"""Track for plotting line"""

	def __init__(self, y, ylim=None, h=1.,fontsize=FONTSIZE, **kw):
		self.h = h
		self.y = y
		self.ylim = ylim
		self.kw = kw
		self.fs = fontsize

	def draw(self, ax):
		ax.plot(self.y, **self.kw)
		if self.ylim:
			ax.set_ylim(self.ylim)
		xmin,xmax = ax.get_xlim()
		ymin,ymax = ax.get_ylim()
		ax.text(xmax,ymax,'{0:.1g}'.format(ymax), ha='right',va='top',fontsize=self.fs)


class Image(Track):
	"""Track for plotting image"""

	def __init__(self, z, h=0.2, zlim=None, cm='jet', **kw):
		self.h = h
		self.z = z
		self.zlim = zlim
		self.cm = cm
		self.kw = kw

	def draw(self, ax):
		if self.zlim is not None:
			vmin,vmax=self.zlim
			ax.imshow(self.z, aspect='auto', interpolation='nearest', cmap=self.cm, vmin=0, vmax=vmax)
		else:
			ax.imshow(self.z, aspect='auto', interpolation='nearest', cmap=self.cm)
		
class TrackWrapper(Track):

	def __init__(self, func, **args):
		self.func = func
		self.args = args

	def draw(self, ax):
		args = self.args
		args['ax'] = ax
		self.func(**args)
		


