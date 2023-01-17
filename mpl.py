#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------
  Custom Matplotlib Presets
-----------------------------

Created on Fri Feb 14 18:14:18 2020
@author: chris
"""



import matplotlib
#matplotlib.use('Qt5Agg')

rcp = matplotlib.rcParams
keys = list(rcp.keys())
#default = [rcp[key] for key in keys]

def set_fontsize(pt):
    rcp['font.size'] = pt
    rcp['legend.fontsize'] = pt-1
    rcp['xtick.labelsize'] = pt-1
    rcp['ytick.labelsize'] = pt-1

def grid(boolean=True):
    rcp['axes.grid'] = boolean

def simple_legend():
    rcp['legend.fancybox'] = False
    rcp['legend.facecolor'] = 'none'
    rcp['legend.frameon'] = False
    rcp['legend.edgecolor'] = 'none'

#%% ----- Global settings ----- %%#

# Default (non-LaTeX) fonts
rcp['font.sans-serif'][0] = 'Lucida Grande'
rcp['font.sans-serif'][3] = 'DejaVu Sans'
#print(rcp['font.sans-serif'])

# Default font sizes
set_fontsize(11)
rcp['axes.titlesize'] = 'medium'

# Grid settings
grid(True)
rcp['grid.linestyle'] = ':'
rcp['grid.linewidth'] = 0.5

# plot settings
rcp['lines.linewidth'] = 1

# Fig and box sizes for Din A4 paper
sizeSingle = [4.8, 3.4]
sizeDual = [2.8, 2.1]
#boxdefault = {1: (4.2,2.55), 2: (2.1, 1.3)}
boxdefault = {1: (4.2,3.0), 2: (2.38, 1.7)}




#%% ----- Utility functions ----- %%#

def find(s,boolean=False):
    """ find rcParam attribute names which contain a keyword

    lists all instances of matplotlib's rcParams that contain the string <s>.

    returns: (if <boolean>=True) a list of all found instances, else nothing.
    """

    WIDTH1 = 32
    WIDTH2 = 16
    WIDTH3 = 48
    N = 0
    f = []
    
    for name in keys:
        if s in name:
            val = rcp[name]
            valtype = str(type(val))
            f.append(name)
            # format output
            if len(name) < WIDTH1: name = name + ':' + (WIDTH1 - 1 - len(name))*' '
            else: name = name + ':'
            if len(valtype) < WIDTH2: valtype = valtype + (WIDTH2 - len(valtype))*' '
            val = str(val)
            if len(val) > WIDTH3+3: val = val[:WIDTH3] + '...'
            # print output
            print('%2i %s %s %s' % (N, name, valtype, val))
            N += 1
    if boolean: return(f)



# set of 10 good-looking markers styles
m0 = {'marker' : 'o',
      'ms' : 5}
m1 = {'marker' : 's',
      'ms' : 4}
m2 = {'marker' : 'D',
      'ms' : 4}
m3 = {'marker' : '^',
      'ms' : 5}
m4 = {'marker' : 'v',
      'ms' : 5}
m5 = {'marker' : '<',
      'ms' : 5}
m6 = {'marker' : '>',
      'ms' : 5}
m7 = {'marker' : 'x',
      'ms' : 6,
      'mew' : 1.5}
m8 = {'marker' : '+',
      'ms' : 7,
      'mew' : 2}
m9 = {'marker' : '*',
      'ms' : 8,
      'mew' : 0.5}

#import itertools
#marker = itertools.cycle((m0,m1,m2,m3,m4,m5,m6,m7,m8,m9))

def marker(i):
    markers = [m0,m1,m2,m3,m4,m5,m6,m7,m8,m9]
    return(markers[i%10]["marker"])

def pt(i):
    """ return a plot marker style in dict format
    
    by default, matplotlib does not cycle through different marker/point types
    when plotting multiple lines.
    
    returns: a dict with the attributes of marker type <i>%10, as defined above.
    """

    pts = [m0,m1,m2,m3,m4,m5,m6,m7,m8,m9]
    return(pts[i%10])

def color(i):
    """
    same as pt(i) but for colors in the cycle
    """
    colors = rcp['axes.prop_cycle'].by_key()['color']
    return(colors[i%10])

def line(i):
    """
    same as pt(i) but for linestyle in the cycle
    """
    lines = ['solid', 'dashed', 'dashdot', 'dotted', (0, (4, 1, 1, 1, 1, 1)), 
             'solid', 'dashed', 'dashdot', 'dotted', (0, (4, 1, 1, 1, 1, 1))]
    return(lines[i%10])

def cycle():
    """ automatically cycle through plot markers and dash types

    if called, matplotlib behaves similar to gnuplot, where different linestyles 
    are automatically cycled through when plotting multiple graphs. 
    """

    from cycler import cycler
    rcp['axes.prop_cycle'] = cycler(
        color = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                 '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'],
        marker = ['o', 's', 'D', '^', 'v', '<', '>', 'x', '+', '*'],
        markersize = [5, 4, 4, 5, 5, 5, 5, 6, 7, 8],
        markeredgewidth = [1, 1, 1, 1, 1, 1, 1, 1.5, 2, 0.5],
        linestyle = ['solid', 'dashed', 'dashdot', 'dotted', (0, (4, 1, 1, 1, 1, 1)), 
                     'solid', 'dashed', 'dashdot', 'dotted', (0, (4, 1, 1, 1, 1, 1))]
        )

def label(obj, lab, xoffset=0, yoffset=0, **kwargs):
    #if offset == 'auto': ypos = 1.2*obj.subplotpars.top - 0.2
    #elif isinstance(offset,float): ypos = 1-offset
    #else: ypos = 0.995
    if isinstance(obj, matplotlib.figure.Figure):
        ax = obj.get_axes()[0]
        ax.text(0.005+xoffset, 0.995-yoffset, '\\Large{'+lab+'}',
                ha='left', va='top', transform=obj.transFigure, **kwargs)
    elif isinstance(obj, matplotlib.axes.Axes):
        obj.text(0.03+xoffset, 0.96-yoffset, lab,
                ha='left', va='top', transform=obj.transAxes, **kwargs)


def setsci(axis, xy='y', which='major'):
    formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(True) 
    formatter.set_powerlimits((-2,2)) 
    if ('x' in xy) or ('X' in xy):
        if which=='major' or which=='both':
            axis.xaxis.set_major_formatter(formatter) 
        if which=='minor' or which=='both':
            axis.xaxis.set_minor_formatter(formatter) 
    if ('y' in xy) or ('Y' in xy):
        if which=='major' or which=='both':
            axis.yaxis.set_major_formatter(formatter) 
        if which=='minor' or which=='both':
            axis.yaxis.set_minor_formatter(formatter) 

def setnonsci(axis, xy='y', which='major'):
    formatter = matplotlib.ticker.ScalarFormatter(useMathText=True)
    formatter.set_scientific(False) 
    if ('x' in xy) or ('X' in xy):
        if which=='major' or which=='both':
            axis.xaxis.set_major_formatter(formatter) 
        if which=='minor' or which=='both':
            axis.xaxis.set_minor_formatter(formatter) 
    if ('y' in xy) or ('Y' in xy):
        if which=='major' or which=='both':
            axis.yaxis.set_major_formatter(formatter) 
        if which=='minor' or which=='both':
            axis.yaxis.set_minor_formatter(formatter) 


def move_sn_y(axis, offs=0, dig=0, side='left', omit_last=False):
    """Move scientific notation exponent from top to the side.
    
    Additionally, one can set the number of digits after the comma
    for the y-ticks, hence if it should state 1, 1.0, 1.00 and so forth.
    From: https://werthmuller.org/blog/2014/move-scientific-notation/

    Parameters
    ----------
    offs : float, optional; <0>
        Horizontal movement additional to default.
    dig : int, optional; <0>
        Number of decimals after the comma.
    side : string, optional; {<'left'>, 'right'}
        To choose the side of the y-axis notation.
    omit_last : bool, optional; <False>
        If True, the top y-axis-label is omitted.

    Returns
    -------
    locs : list
        List of y-tick locations.

    Note
    ----
    This is kind of a non-satisfying hack, which should be handled more
    properly. But it works. Functions to look at for a better implementation:
    ax.ticklabel_format
    ax.yaxis.major.formatter.set_offset_string
    """

    # Get the ticks
    locs, _ = axis.get_yticks()

    # Put the last entry into a string, ensuring it is in scientific notation
    # E.g: 123456789 => '1.235e+08'
    llocs = '%.3e' % locs[-1]

    # Get the magnitude, i.e. the number after the 'e'
    # E.g: '1.235e+08' => 8
    yoff = int(str(llocs).split('e')[1])

    # If omit_last, remove last entry
    if omit_last:
        slocs = locs[:-1]
    else:
        slocs = locs

    # Set ticks to the requested precision
    form = r'$%.'+str(dig)+'f$'
    axis.set_yticks(locs, list(map(lambda x: form % x, slocs/(10**yoff))))

    # Define offset depending on the side
    if side == 'left':
        offs = -.18 - offs # Default left: -0.18
    elif side == 'right':
        offs = 1 + offs    # Default right: 1.0
        
    # Plot the exponent
    axis.text(offs, .98, r'$\times10^{%i}$' % yoff, transform =
            axis.transAxes, verticalalignment='top')

    # Return the locs
    return locs


def format(fig, boxsize=None, ylabel=1, xlabel=1, yticl=1, xticl=1, title=False,
           left=0, right=0, bottom=0, top=0):
    # convert to inches
    fontheight = 0.0138889*rcp['font.size']
    fontwidth = 0.0138889*rcp['font.size']

    # determine size of spine box
    if boxsize is None: boxsize = 1
    if isinstance(boxsize, int): 
        boxsize = boxdefault[boxsize]
    elif (not isinstance(boxsize,tuple)) and (not isinstance(boxsize,list)):
        raise ValueError('boxsize needs to be of type int or tuple/list.')

    # padding for safety
    if bottom==0: bottom = 0.25*fontheight
    if top==0:    top    = 0.25*fontheight
    if left==0:   left   = 0.25*fontheight
    if right==0:  right  = 0.25*fontheight
    # padding for everything outside of spine box
    left   += 3.0*fontheight*yticl
    bottom += 1.2*fontheight*xticl
    left   += 1.2*fontheight*ylabel
    bottom += 1.2*fontheight*xlabel
    top    += 1.2*fontheight*title
    # nice centering of full-page figure
    if boxsize[0]>3.5: right += 0.65*left

    # transform to figure coordinates
    left   = round(left,2)
    right  = round(right,2)
    top    = round(top,2)
    bottom = round(bottom,2)
    figsize = (boxsize[0]+left+right, boxsize[1]+bottom+top)
    bottom = bottom/figsize[1]
    top    = 1 - top/figsize[1]
    left   = left/figsize[0]
    right  = 1 - right/figsize[0]

    # update figure
    fig.set_size_inches(figsize)
    fig.subplots_adjust(bottom=bottom, top=top, left=left, right=right)
    




#%% ----- Presets ----- %%#

def Surf():
    """ Settings that look good on Microsoft Surface Pro 7 """

    rcp['figure.figsize'] = sizeSingle
    rcp['figure.dpi'] = 240
    rcp['axes.linewidth'] = 0.8
    rcp['lines.linewidth'] = 1.5
    # Space
    rcp['figure.subplot.top'] = 0.92


def Desk():
    """ Settings that look good on a 1440p screen """

    rcp['figure.figsize'] = sizeSingle
    rcp['figure.dpi'] = 180
    # Space
    rcp['figure.subplot.top'] = 0.92


def Work():
    """ Settings that look good on a 1080p screen """

    rcp['figure.figsize'] = sizeSingle
    rcp['figure.dpi'] = 150
    # Space
    rcp['figure.subplot.top'] = 0.92
    # Font
    rcp['text.usetex'] = True
    rcp['font.family'] = 'serif'


def Server():
    """ Settings that look good when plotting on the small cluster """

    rcp['figure.figsize'] = sizeSingle
    rcp['figure.dpi'] = 150
    # Space
    rcp['figure.subplot.top'] = 0.92
    # Font
    rcp['text.usetex'] = False
    rcp['font.family'] = 'serif'


def Publ(num=1,ylabel=False,title=False):
    """ Settings that look good in publications

    <num> is the number of figures that need to fit side by side.
    <num> > 1 might require adjustments depending on axis labels.
    """

    if num == 1: 
        rcp['figure.figsize'] = [5.2, 3.4] # single
        rcp['lines.linewidth'] = 1.5
    elif num==2: 
        rcp['figure.figsize'] = [2.6, 2.2] # double
        rcp['lines.linewidth'] = 1.
        set_fontsize(10)
    rcp['figure.dpi'] = 180
    # Space
    if num==1:
        rcp['figure.subplot.top'] = 0.92
        rcp['figure.subplot.left'] = 0.15
        rcp['figure.subplot.right'] = 0.90
        rcp['figure.subplot.bottom'] = 0.12
    elif num==2:
        rcp['figure.subplot.top'] = 0.98
        if ylabel: rcp['figure.subplot.left'] = 0.21
        else: rcp['figure.subplot.left'] = 0.1
        rcp['figure.subplot.right'] = 0.97
        rcp['figure.subplot.bottom'] = 0.18
    # Font
    rcp['text.usetex'] = True
    rcp['font.family'] = 'serif'
    # Ticks
    rcp['xtick.direction'] = 'in'
    rcp['ytick.direction'] = 'in'


def RevTex(num=1, ylabel=False, simple=True):
    if simple: simple_legend()
    if num==1:
        rcp['figure.figsize'] = [3.2, 2.2] # double
        rcp['lines.linewidth'] = 1.
        rcp['figure.subplot.top'] = 0.97
        rcp['figure.subplot.left'] = 0.17
        rcp['figure.subplot.right'] = 0.97
        rcp['figure.subplot.bottom'] = 0.17
    elif num==2:
        rcp['figure.figsize'] = [1.6, 1.0] # double
        rcp['lines.linewidth'] = 1.
        rcp['figure.subplot.top'] = 0.97
        rcp['figure.subplot.left'] = 0.30
        rcp['figure.subplot.right'] = 0.97
        rcp['figure.subplot.bottom'] = 0.34
    # Font
    rcp['text.usetex'] = True
    rcp['font.family'] = 'serif'
    # Ticks
    rcp['xtick.direction'] = 'in'
    rcp['ytick.direction'] = 'in'
    rcp['xtick.minor.visible'] = True
    rcp['ytick.minor.visible'] = True
    # font sizes
    rcp['font.size'] = 10
    rcp['legend.fontsize'] = 9
    rcp['xtick.labelsize'] = 9
    rcp['ytick.labelsize'] = 9
    # legend
    rcp['legend.labelspacing'] = 0.2
    rcp['legend.framealpha'] = 1.0
    # grid
    rcp['axes.grid'] = False
