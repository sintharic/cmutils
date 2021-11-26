# -*- coding: utf-8 -*-
"""
-----------------------------
  Custom Matplotlib Presets
-----------------------------

Created on Fri Feb 14 18:14:18 2020
@author: chris
"""

import matplotlib

rcp = matplotlib.rcParams
keys = list(rcp.keys())
#default = [rcp[key] for key in keys]





#%% ----- Global settings ----- %%#

# Default (non-LaTeX) fonts
rcp['font.sans-serif'][0] = 'Lucida Grande'
rcp['font.sans-serif'][3] = 'DejaVu Sans'
#print(rcp['font.sans-serif'])

# Default font sizes
rcp['font.size'] = 11
rcp['legend.fontsize'] = 10
rcp['xtick.labelsize'] = 10
rcp['ytick.labelsize'] = 10

# Grid settings
rcp['axes.grid'] = True
rcp['grid.linestyle'] = ':'

# Fig sizes for Din A4 paper
sizeSingle = [4.8, 3.4]
sizeDual = [2.8, 2.1]





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
            if len(name) < WIDTH1: name = name + ":" + (WIDTH1 - 1 - len(name))*" "
            else: name = name + ":"
            if len(valtype) < WIDTH2: valtype = valtype + (WIDTH2 - len(valtype))*" "
            val = str(val)
            if len(val) > WIDTH3+3: val = val[:WIDTH3] + "..."
            # print output
            print("%2i %s %s %s" % (N, name, valtype, val))
            N += 1
    if boolean: return(f)



def grid(boolean=True):
    rcp['axes.grid'] = boolean



# set of 10 good-looking markers styles
m0 = {"marker" : "o",
      "ms" : 5}
m1 = {"marker" : "s",
      "ms" : 4}
m2 = {"marker" : "D",
      "ms" : 4}
m3 = {"marker" : "^",
      "ms" : 5}
m4 = {"marker" : "v",
      "ms" : 5}
m5 = {"marker" : "<",
      "ms" : 5}
m6 = {"marker" : ">",
      "ms" : 5}
m7 = {"marker" : "x",
      "ms" : 6,
      "mew" : 1.5}
m8 = {"marker" : "+",
      "ms" : 7,
      "mew" : 2}
m9 = {"marker" : "*",
      "ms" : 8,
      "mew" : 0.5}

#import itertools
#marker = itertools.cycle((m0,m1,m2,m3,m4,m5,m6,m7,m8,m9))


def pt(i):
    """ return a plot marker style in dict format
    
    by default, matplotlib does not cycle through different marker/point types
    when plotting multiple lines.
    
    returns: a dict with the attributes of marker type <i>%10, as defined above.
    """

    pts = [m0,m1,m2,m3,m4,m5,m6,m7,m8,m9]
    return(pts[i%10])


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
        linestyle = ["solid", "dashed", "dashdot", "dotted", (0, (4, 1, 1, 1, 1, 1)), 
                     "solid", "dashed", "dashdot", "dotted", (0, (4, 1, 1, 1, 1, 1))]
        )





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


def Publ(num=1,ylabel=False):
    """ Settings that look good in publications

    <num> is the number of figures that need to fit side by side.
    <num> > 1 might require adjustments depending on axis labels.
    """

    if num == 1: rcp['figure.figsize'] = [5.2, 3.4] # single
    elif num==2: rcp['figure.figsize'] = [2.6, 2.2] # double
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
        rcp['figure.subplot.bottom'] = 0.16
    # Font
    rcp['text.usetex'] = True
    rcp['font.family'] = 'serif'
    # Ticks
    rcp['xtick.direction'] = 'in'
    rcp['ytick.direction'] = 'in'