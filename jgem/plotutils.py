import matplotlib.pyplot as P
import pandas as PD
import numpy as N
from scipy.stats import kde
from matplotlib.colors import LogNorm, ListedColormap

hls10 = [(0.86, 0.37119999999999997, 0.33999999999999997),
         (0.86, 0.68320000000000003, 0.33999999999999997),
         (0.72479999999999989, 0.86, 0.33999999999999997),
         (0.41279999999999994, 0.86, 0.33999999999999997),
         (0.33999999999999997, 0.86, 0.57920000000000016),
         (0.33999999999999997, 0.82879999999999987, 0.86),
         (0.33999999999999997, 0.51679999999999948, 0.86),
         (0.47520000000000029, 0.33999999999999997, 0.86),
         (0.7871999999999999, 0.33999999999999997, 0.86),
         (0.86, 0.33999999999999997, 0.62079999999999991)]

hls8 =  [(0.86, 0.37119999999999997, 0.33999999999999997),
         (0.86, 0.7612000000000001, 0.33999999999999997),
         (0.56880000000000008, 0.86, 0.33999999999999997),
         (0.33999999999999997, 0.86, 0.50120000000000009),
         (0.33999999999999997, 0.82879999999999987, 0.86),
         (0.33999999999999997, 0.43879999999999986, 0.86),
         (0.63119999999999976, 0.33999999999999997, 0.86),
         (0.86, 0.33999999999999997, 0.69879999999999964)]

hls6 =  [(0.86, 0.37119999999999997, 0.33999999999999997), 
        (0.82879999999999987, 0.86, 0.33999999999999997), 
        (0.33999999999999997, 0.86, 0.37119999999999997), 
        (0.33999999999999997, 0.82879999999999987, 0.86), 
        (0.37119999999999997, 0.33999999999999997, 0.86), 
        (0.86, 0.33999999999999997, 0.82879999999999987)]

hls5 =  [(0.86, 0.37119999999999997, 0.33999999999999997),
         (0.72479999999999989, 0.86, 0.33999999999999997),
         (0.33999999999999997, 0.86, 0.57920000000000016),
         (0.33999999999999997, 0.51679999999999948, 0.86),
         (0.7871999999999999, 0.33999999999999997, 0.86)]

hls4 = [(0.86, 0.37119999999999997, 0.33999999999999997),
         (0.56880000000000008, 0.86, 0.33999999999999997),
         (0.33999999999999997, 0.82879999999999987, 0.86),
         (0.63119999999999976, 0.33999999999999997, 0.86)]

set2 = [(0.40000000596046448, 0.7607843279838562, 0.64705884456634521),
         (0.98131487965583808, 0.55538641635109398, 0.38740485135246722),
         (0.55432528607985565, 0.62711267120697922, 0.79595541393055635),
         (0.90311419262605563, 0.54185316071790801, 0.76495195557089413),
         (0.65371782148585622, 0.84708959004458262, 0.32827375098770734),
         (0.9986312957370983, 0.85096502233954041, 0.18488274134841617),
         (0.89573241682613591, 0.76784315109252932, 0.58182240093455595),
         (0.70196080207824707, 0.70196080207824707, 0.70196080207824707)]


hls10cm = ListedColormap(hls10, name='hls10')
hls8cm = ListedColormap(hls8, name='hls8')
hls6cm = ListedColormap(hls6, name='hls6')
hls4cm = ListedColormap(hls4, name='hls4')
set2cm = ListedColormap(set2, name='set2')


####################### util for axes within subplots
# from http://stackoverflow.com/questions/17458580/embedding-small-plots-inside-subplots-in-matplotlib
def add_subplot_axes(ax,rect,axesargs={'axisbg':'w'},noticks=True):
    # rect: [x,y,width,height]
    fig = ax.get_figure()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]
    subax = fig.add_axes([x,y,width,height],**axesargs)
    if not noticks:
        x_labelsize = subax.get_xticklabels()[0].get_size()
        y_labelsize = subax.get_yticklabels()[0].get_size()
        x_labelsize *= rect[2]**0.5
        y_labelsize *= rect[3]**0.5
        subax.xaxis.set_tick_params(labelsize=x_labelsize)
        subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax
#######################


def broken_yaxes(ax, size1=0.8, space=0.03, d=.015):
    # http://matplotlib.org/examples/pylab_examples/broken_axis.html
    h1 = size1 - space/2.
    y2 = h1 + space
    h2 = 1.-y2
    ax1 = add_subplot_axes(ax, [0,0,1,h1],noticks=False)
    ax2 = add_subplot_axes(ax, [0,y2,1,h2],noticks=False)
    ax2.spines['bottom'].set_visible(False)
    ax1.spines['top'].set_visible(False)
    ax2.xaxis.tick_top()
    ax2.tick_params(labeltop='off')
    ax1.xaxis.tick_bottom()
    #d = .015  # how big to make the diagonal lines in axes coordinates
    # arguments to pass plot, just so we don't keep repeating them
    kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
    ax1.plot((-d, +d), (1 - d, 1 + d), **kwargs)  # bottom-left diagonal
    ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)  # bottom-right diagonal    
    kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
    tmp1 = ax1.transAxes.transform(([0,0],[d,d])) # in display coord
    tmp2 = ax2.transAxes.inverted().transform(tmp1) # in ax2 coordinate
    dy = tmp2[1][1]-tmp2[0][1]
    ax2.plot((-d, +d), (-dy, +dy), **kwargs)        # top-left diagonal
    ax2.plot((1 - d, 1 + d), (-dy, +dy), **kwargs)  # top-right diagonal
    P.setp(ax, frame_on=False, xticks=[], yticks=[])
    return ax1, ax2

def plot_stacked(df, ax=None, 
                 plotlegend=True, 
                 cls=hls8, 
                 percent=False, 
                 title='',
                 grid=False,
                 dx=0.5,
                 ylim=None,
                 kind='bar',
                 labelspacing=0.15):
    n0,n1 = df.shape
    if ax is None:
        fig,ax = P.subplots(figsize=(dx*n1,3))
    if percent:
        tmpn = 100*(df.div(df.sum(axis=0),axis=1))
    else:
        tmpn = df
    
    tmpn.T.plot(kind=kind,stacked=True,legend=False,width=1,color=cls,ax=ax,grid=grid)
    if kind=='bar':
        if percent:
            ax.set_ylim([0,100])
        elif ylim is not None:
            if ylim=='max':
                ylim = (0,df.sum().max())
            print('ylim={0}'.format(ylim))
            ax.set_ylim(ylim)
        ax.set_xlim([0,n1-1])
    else:
        if percent:
            ax.set_xlim([0,100])
        ax.set_ylim([-1,n1-1])
    if plotlegend:
        h,l = ax.get_legend_handles_labels()
        ax.legend(h[::-1], l[::-1], 
                  loc='center left', 
                  bbox_to_anchor=(1.0,0.5), 
                  labelspacing=labelspacing)
    ax.set_title(title)
    return ax

def take_top_n(df, n):
    df0 = df.iloc[:n]
    df0.ix['other'] = df.iloc[n:].sum(axis=0)
    print(df0.sum())
    return df0


def plot_panel(panel, usize=1, takelog=True):
    n0,n1,n2=panel.shape # ths, tcodes, scodes
    fig,axr = P.subplots(n2, n1, figsize=(usize*n1,usize*n2),sharex=True,sharey=True)
    P.subplots_adjust(wspace=0.15,hspace=0.15)
    xs = panel.items.values # thresholds
    yls = panel.major_axis.values # tcodes
    zls = panel.minor_axis.values # scodes
    for i,s in enumerate(zls): # n2 (top to down)
        for j,t in enumerate(yls): # n1 (left to right)
            ax = axr[i][j]
            ys = panel[xs,t,s]
            if takelog:
                ys = N.log10(ys+1)
            ax.plot(xs,ys,'.-')
            if i==n2-1:
                ax.set_xlabel(t)
            if j==0:
                ax.set_ylabel(s)
    return fig


def plot_density(xv,yv,ax=None,xlim=None,ylim=None,nbins=50, alpha=0.9,
                 log=False, vmin=None, vmax=None, contour=True, n=6):
    if ax is None:
        fig,ax = P.subplots(figsize=(3,3))
    c = N.corrcoef(xv,yv)
    print('corr.coef.: {0}'.format(c[0][1]))
    k = kde.gaussian_kde([xv,yv])
    if xlim is None:
        xlim = N.min(xv),N.max(xv)
    if ylim is None:
        ylim = N.min(yv),N.max(yv)
    xmi,xma = xlim
    ymi,yma = ylim
    xi, yi = N.mgrid[xmi:xma:nbins*1j, ymi:yma:nbins*1j]
    zi = k(N.vstack([xi.flatten(), yi.flatten()]))
    if log:
        if vmin is None:
            vmin = N.min(zi)
        if vmax is None:
            vmax = N.max(zi)
        print( N.min(zi), N.max(zi))
        zi = zi.reshape(xi.shape)
        ax.pcolor(xi, yi, zi, norm=LogNorm(vmin=vmin,vmax=vmax),alpha=alpha)
        if contour:
            ax.contour(xi,yi,zi,n)
    else:
        zi = zi.reshape(xi.shape)
        ax.pcolormesh(xi, yi, zi, alpha=alpha)
        if contour:
            ax.contour(xi,yi,zi,n)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    