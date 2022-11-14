#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
------------------------------------
  Post-Processing of contMech Data
------------------------------------

@author: thescientist
"""



#%% ----- User settings ----- %%#

XLIM = (-0.5,0.5)   # default xrange
YLIM = (-0.5,0.5)   # default yrange
ZLIM = "auto"       # default zrange
SHOW = True         # show plot windows
FMT  = "%.6g"       # number format used to write data to files

WARN = True         # display non-critical warnings and reminders





#%% ----- Plotting ----- %%#

import numpy as np
from numba import njit, double
import matplotlib.pyplot as plt
#from tqdm import tqdm # tqdm 3rdPartyModule


def show(): plt.show()


def plotImg(array, clim=ZLIM, title=False, aspect="true", cmap="turbo", 
            alpha=1.0, axis=None, pad=0.01, cba=0.85, rot=False, **kwargs):
    """ plot array as an image
    
    plots 2D numpy.ndarray <array> as a colored image, similar to gnuplot's 
    "splot ... with pm3d".
    <clim> must be of form (v_min,v_max), the value range to which colors are mapped.

    returns: matplotlib.axes.Axes object containing the plot
    """
    if rot: 
        origin = "lower"
        arrayLoc = array.transpose()
    else: 
        origin = "upper"
        arrayLoc = array

    if SHOW: plt.ion(); plt.show()
    if aspect == "true": aspect = (XLIM[1]-XLIM[0])/(YLIM[1]-YLIM[0])
    if axis==None:
        plt.figure()
        plt.grid(False)
        plt.imshow(arrayLoc,cmap=cmap,aspect=aspect,alpha=alpha,origin=origin,**kwargs)
        if clim!="auto": plt.clim(clim[0],clim[1])
        cbar = plt.colorbar(shrink=cba,aspect=20*cba,pad=0.01)
        cbar.formatter.set_powerlimits((-2, 2))
        cbar.ax.yaxis.set_offset_position('left')
    else:
        im = axis.imshow(arrayLoc,cmap=cmap,aspect=aspect,alpha=alpha,origin=origin,**kwargs)
        if clim!="auto": im.set_clim(clim[0],clim[1])
    
    if title: plt.title(title)
    if SHOW: plt.pause(0.001)
    return(plt.gca())


def plotLines(array, lines, axis=None, dim=0, ls="-", **kwargs):
    """ plot line scan(s)
    
    plots line(s) along dimension <dim> from numpy.ndarray <array>.
    <lines> can be an int or a list/tuple of ints with desired line indices.
    <axis>, the matplotlib.axes.Axes object on which to plot, can be specified.
    as well as typical kwargs for plots.

    returns: matplotlib.axes.Axes object containing the plot
    """

    if dim==0: xlab="Y"
    else: xlab="X"
    if type(lines) == int: lines = [lines]
    
    if SHOW: plt.ion(); plt.show()
    if not axis: plt.figure()
    for line in lines:
        if dim==0:
            z = array[line,:]
            x = np.linspace(YLIM[0],YLIM[1],z.size)
        else:
            z = array[:,line]
            x = np.linspace(XLIM[0],XLIM[1],z.size)
        
        if axis: axis.plot(x,z,ls,label="line "+str(line),lw=0.9,**kwargs)
        else: plt.plot(x,z,ls,label="line "+str(line),lw=0.9,**kwargs)
    if not axis:
        plt.xlabel(xlab)
        plt.ylabel("Height")
    plt.legend()
    if SHOW: plt.pause(0.001)
    return(plt.gca())


def plotSurf(array, kind="color", title=False, clim=ZLIM, axis=None, 
             figsize=None, stride=4, cmap="turbo", lw=0.5, lc="gray", aa=False,
             alpha=1, label=None):
    """ plot 3D view of surface
    
    plots 2D numpy.ndarray <array>, using the data as z values, in a 3D view.
    
    for <kind>="color", the color map <cmap> can be specified, as well as the
    color range <clim> in the form (v_min,v_max).
    <stride> specifies the discretization increment, lower being higher quality.

    for <kind>="wire", the linewidth <lw> and linecolor <lc> can be specified.

    returns: matplotlib.axes.Axes object containing the plot
    """
    
    from mpl_toolkits.mplot3d import Axes3D # noqa: F401 unused import

    if SHOW: plt.ion(); plt.show()
    if axis==None:
        if figsize: fig = plt.figure(figsize=figsize)
        else: fig = plt.figure(figsize=figsize)
        ax = fig.gca(projection='3d')
    else: ax = axis
    X = np.arange(0,array.shape[0],1)
    Y = np.arange(0,array.shape[1],1)
    X, Y = np.meshgrid(X,Y)
    if clim=="auto": clim = (array.min(), array.max())
    if kind=="color": 
        if type(cmap) != str: 
            col = cmap
            cmap = None
        else: 
            col = None
        surf = ax.plot_surface(X,Y, np.transpose(array), cmap=cmap, color=col,
                               rstride=stride, cstride=stride, vmin=clim[0], vmax=clim[1], 
                               linewidth=0, antialiased=aa, alpha=alpha)
    elif kind=="wire": 
        surf = ax.plot_wireframe(X,Y, np.transpose(array), rstride=stride, cstride=stride,
                                 linewidth=lw, color=lc, antialiased=aa, alpha=alpha)
    ax.axis("off")
    #ax.set_box_aspect((np.ptp(X), np.ptp(Y), np.ptp(array)))
    ax.set_zlim(clim[0],clim[1])
    if title: plt.title(title)
    if kind=="color" and axis==None: 
        cb = fig.colorbar(surf, shrink=0.8, aspect=8)
        if label: 
            cb.ax.set_title(label)
            cb.ax.xaxis.set_label_position('top')
    if SHOW: plt.pause(0.001)
    return(ax)
    




#%% ----- Data input / output / manipulation ----- %%#

def readConfig(filepath, usecols=2):
    """ read contMech konfig file
    
    reads column <usecols> from file <filepath>, and converts it to an (nx,ny)
    numpy.ndarray, assuming the first line in the file is of the form "#nx ny"

    returns: numpy.ndarray
    """

    global XLIM, YLIM

    with open(filepath) as file:
        line1 = file.readline()[1:].split()
        lines = [line for line in file.readlines() if line[0].isalnum()]
    
    nx, ny = ( int(line1[0]), int(line1[1]) )
    array = np.array([float(line.split()[usecols]) for line in lines])
    lengthX = float(lines[-1].split()[0])
    lengthY = 2*float(lines[ny//2].split()[1])

    # default: reshape the (nx+1)*(ny+1) flat array to a (nx, ny) ndarray
    if (array.shape[0] == (nx+1)*(ny+1)):
        array = array.reshape((nx+1,ny+1))[:nx,:ny]
    # frames: might have been downsampled to 512x512
    elif ( (array.shape[0] == 513*513) & (array[0] == array[512]) ):
        nx, ny = (512, 512)
        array = array.reshape((nx+1,ny+1))[:nx,:ny]
    # legacy: cmutils used to export nx*ny flat array rather than (nx+1)*(ny+1)
    else: 
        array = array.reshape((nx,ny))
        lengthX = lengthX*nx/(nx-1)

    # update XLIM, YLIM
    XLIM = (-lengthX/2, lengthX/2)
    YLIM = (-lengthY/2, lengthY/2)
    if WARN: 
        print("[Warn:Dimens] Updating XLIM and YLIM to fit imported data:",flush=True)
        print("  lengthX  = ", lengthX, flush=True)
        print("  lengthY  = ", lengthY, flush=True)

    return(array)
    
    

def readConfigOld(filepath, usecols=2):
    """ LEGACY VERSION: read contMech konfig file
    
    reads column <usecols> from file <filepath>, and converts it to an (nx,ny)
    numpy.ndarray, assuming the first line in the file is of the form "#nx ny"

    returns: numpy.ndarray
    """

    global XLIM, YLIM
    array = np.loadtxt(filepath, usecols=usecols)
    
    # Read nx, ny and reshape array accordingly
    fid = open(filepath,"rb")
    line1 = fid.readline().decode('UTF-8')
    line1 = line1[1:].split()
    nx, ny = ( int(line1[0]), int(line1[1]) )
    if ( (array.shape[0] == (nx+1)*(ny+1)) ):# & (array[0] == array[ny]) ): # additional lines
        lXfac = 1
        array = array.reshape((nx+1,ny+1))
        array = array[:nx,:ny]
    elif ( (array.shape[0] == 513*513) & (array[0] == array[512]) ): # frames are downsampled to 512x512
        lXfac = 1
        nx, ny = (512, 512)
        array = array.reshape((nx+1,ny+1))
        array = array[:nx,:ny]
    else:
        lXfac = nx/(nx-1)
        array = array.reshape((nx,ny))
    
    # Update YLIM using first 2 lines
    dummy = np.loadtxt(filepath, usecols=(0,1), max_rows=2)
    dyH = (dummy[1,1] - dummy[0,1])/2
    YLIM = (-dyH*ny, dyH*ny)
    
    # Update XLIM using first and last line
    fid.seek(-2, 2) # 2nd arg = 2 --> relative to EOF
    while fid.read(1) == b"\n": fid.seek(-2, 1) # skip trailing empty lines
    fid.seek(-1, 1) # 2nd arg = 1 --> relative to current position
    while fid.read(1) != b"\n": fid.seek(-2, 1)
    lineN = fid.read().decode('UTF-8')
    idx = lineN.find("\t")
    lengthX = ( float(lineN[:idx]) - dummy[0,0] )*lXfac
    XLIM = (-lengthX/2, lengthX/2)
    fid.close()
    #array = np.rot90(array)

    if WARN: 
        print("[Warn:Dimens] Updating XLIM and YLIM to fit imported data:",flush=True)
        print("  lengthX  = ", lengthX, flush=True)
        print("  lengthY  = ", 2*dyH*ny, flush=True)
    return(array)
    

def readImg(filepath):
    """ read image file (without lateral info)
    
    reads file <filepath> with ny numbers per line and nx lines into a (nx, ny)
    numpy.ndarray.

    returns: numpy.ndarray
    """

    if WARN: print("[Warn:Format] contMech format is upside-down! Use img2config() to compare.",flush=True)
    data = np.loadtxt(filepath)
    data = np.transpose(data)
    return(data)


def img2config(array): return(array.max() - array)

def config2img(array): return(-array)


def resample(array, resol):
    """ resample / change array resolution
    
    resamples the (nx,ny) numpy.ndarray <array> to the new resolution <resol>.
    <resol> can be of form (nx_new, ny_new) or a single float, 
    which is translated to <resol> = (round(<resol>*nx), round(<resol>*ny)).

    returns: resampled numpy.ndarray
    """
    
    from scipy import signal

    if WARN: 
        print("[Warn:Period] Resampling uses Fourier filter, which assumes periodicity.",flush=True)
        print("              For non-periodic arrays, use bilin_resample() instead.",flush=True)
    if type(resol)==int: newShape = tuple([resol*iShape for iShape in array.shape])
    else: newShape = resol
    if len(newShape) != len(array.shape):
        print("[ERROR] 1st and 2nd argument must have same number of dimensions.")
        return(array)
    result = np.copy(array)
    for dim in range(len(newShape)): result = signal.resample(result, newShape[dim], axis=dim)
    return(result)


@njit
def bilin_resample(array, resol):
    """ resample / change array resolution through bilinear interpolation

    resamples the (nx,ny) numpy.ndarray <array> to the new resolution <resol>.
    <resol> can be of form (nx_new, ny_new) or a single int, 
    which is translated to <resol> = (round(nx/<resol>), round(ny/<resol>)).

    returns: resampled numpy.ndarray
    """
    nxNew = resol[0]
    nyNew = resol[1]
    resolX = array.shape[0]/nxNew
    resolY = array.shape[1]/nyNew
  
    result = np.zeros((nxNew,nyNew),dtype=double)
    
    if (array.shape[0]/nxNew == int(array.shape[0]/nxNew)) and (array.shape[1]/nyNew == int(array.shape[1]/nyNew)):
        xIdx = np.arange(nxNew+1)*resolX
        xIdx = xIdx.astype(np.uint16)
        yIdx = np.arange(nyNew+1)*resolY
        yIdx = yIdx.astype(np.uint16)
        
        for ix in range(nxNew):
            for iy in range(nyNew):
                result[ix,iy] = np.median(array[xIdx[ix]:xIdx[ix+1], yIdx[iy]:yIdx[iy+1]])
    
    else:
        xIdx = np.arange(nxNew)*resolX
        yIdx = np.arange(nyNew)*resolY
    
        for ix in range(nxNew):
            for iy in range(nyNew):
                # neighbors
                xL = int(xIdx[ix])
                xR = xL + 1
                yB = int(yIdx[iy])
                yT = yB + 1
        
                # weights: bilinear interpolation
                dx = xIdx[ix] - xL
                dy = yIdx[iy] - yB
                wTR = dx*dy
                wTL = (1-dx)*dy
                wBR = dx*(1-dy)
                wBL = (1-dx)*(1-dy)
                val = wBL*array[xL,yB] + wBR*array[xR,yB] + wTL*array[xL,yT] + wTR*array[xR,yT]
        
                result[ix,iy] = val
  
    return(result)


def reduce(array, resol):
    """ resample / change array resolution

    resamples the (nx,ny) numpy.ndarray <array> to the new resolution <resol>.
    <resol> can be of form (nx_new, ny_new) or a single int, 
    which is translated to <resol> = (round(nx/<resol>), round(ny/<resol>)).

    returns: resampled numpy.ndarray
    """

    if WARN: print("[Warn:Deprec] reduce() is deprecated and will be removed in a future release. Use bilin_resample() instead.",flush=True)
    if type(resol)==int:
        nxNew = round(array.shape[0]/resol)
        nyNew = round(array.shape[1]/resol)
        return resample(array, (nxNew,nyNew))
    else:
        nxNew, nyNew = resol
        return resample(array, (nxNew,nyNew))


def slopeX(array):
    dx = (XLIM[1]-XLIM[0])/array.shape[0]
    return ( np.diff(array, append=array[:1,:], axis=0) + np.diff(array, prepend=array[-1:,:], axis=0) ) / (2*dx)

def slopeY(array):
    dy = (YLIM[1]-YLIM[0])/array.shape[1]
    return ( np.diff(array, append=array[:,:1], axis=1) + np.diff(array, prepend=array[:,-1:], axis=1) ) / (2*dy)


def corners(array, N=256):
    """ meeting point of <array>'s corners

    constructs the 2<N> x 2<N> area around the point where the four corners meet 
    if np.ndarray <array> is periodically repeated in the plane.

    returns: np.ndarray
    """

    result = np.zeros((2*N,2*N))
    result[:N,:N] = array[-N:,-N:]
    result[N:,N:] = array[:N,:N]
    result[N:,:N] = array[:N,-N:]
    result[:N,N:] = array[-N:,:N]
    return(result)


def smoothPBC(array, overlap=None):
    """ make non-periodic image smooth at periodic boundaries

    all points within <overlap> from the edges of np.ndarray <array> are 
    reassigned values that gradually approach the value of the opposite edge.

    returns: np.ndarray
    """

    if overlap is None: overlap = min(array.shape)//10
    result = np.copy(array)

    # array axis 0: use original edges 
    for idx in range(overlap):
        wt_idx = 0.75 - 0.25*np.cos(np.pi*idx/overlap)
        wt_edge = 1 - wt_idx # 0.25 + 0.25*np.cos(np.pi*idx/overlap)
        result[idx,:] = wt_idx*array[idx,:] + wt_edge*array[-1,:]
        result[-(idx+1),:] = wt_idx*array[-(idx+1),:] + wt_edge*array[0,:]

    # array axis 1: use current edges
    edge0 = result[:,0]
    edge1 = result[:,-1]
    for idx in range(overlap):
        wt_idx = 0.75 - 0.25*np.cos(np.pi*idx/overlap)
        wt_edge = 1 - wt_idx # 0.25 + 0.25*np.cos(np.pi*idx/overlap)
        result[:,idx] = wt_idx*result[:,idx] + wt_edge*edge1
        result[:,-(idx+1)] = wt_idx*result[:,-(idx+1)] + wt_edge*edge0

    return(result)


def dumpConfig(arrays, filepath="konfig0py.real", Lx=None, Ly=None):
    """ dump arrays as contMech konfig

    writes list <arrays> of numpy.ndarrays to a file in the typical format 
    that can be plotted in gnuplot or imported into a contMech simulation.
    """
    if type(arrays)==np.ndarray: arrays = [arrays]

    if WARN: print("[Warn:Format] for imported microscope images, use img2config() first.",flush=True)
    (nx,ny)=arrays[0].shape #array.shape
    if not Lx: dx = (XLIM[1]-XLIM[0])/nx
    else: dx = Lx/nx
    if not Ly: dy = (YLIM[1]-YLIM[0])/ny
    else: dy = Ly/ny
    out = open(filepath,"w")
    out.write("#%i\t%i\n\n" % (nx,ny))
    #for ix in tqdm(range(nx)): # tqdm 3rdPartyModule
    for ix in range(nx+1):
        for iy in range(ny+1):
            s = (FMT+"\t"+FMT) % (ix*dx, iy*dy)
            for array in arrays: s += ("\t"+FMT) % array[ix%nx,iy%ny]
            s += "\n"
            out.write(s)
            #out.write( (FMT+"\t"+FMT+"\t"+FMT+"\n") % (ix*dx, iy*dy, array[ix%nx,iy%ny]) )
        out.write("\n") # divider for gnuplot
    out.close()


def dumpImg(array, filename):
    """ dump array as image text file

    writes 2D (nx, ny) numpy.ndarray <array> to a text file with ny values per
    line and nx lines.
    """

    out = open(filename,"w")
    #for ix in tqdm(range(array.shape[0])): # tqdm 3rdPartyModule
    for ix in range(array.shape[0]):
        for iy in range(array.shape[1]-1):
            out.write( (FMT+"\t") % (array[ix,iy]) )
        out.write( (FMT+"\n") % (array[ix,array.shape[1]-1]) )
    out.close()


def dumpLines(array, lines, filename="", dim=0):
    """ Dump line scan(s)
    
    exports individual lines along dimension <dim> from 2D numpy.ndarray <array> 
    into <filename>. 
    <lines> can be an int or a list/tuple of ints with desired line indices to
    be exported.
    """

    if type(lines)==int: lines = [lines]
    if dim==0:
        z = array[lines,:]
        #print(z.shape)
        x = np.linspace(YLIM[0],YLIM[1],z.shape[1])
        label = ["Y","X"]
    else:
        z = array[:,lines]
        z = np.transpose(z)
        #print(z.shape)
        x = np.linspace(XLIM[0],XLIM[1],z.shape[1])
        label = ["X","Y"]
    if filename=="": filename = "lineScans%i.dat"%dim
    out = open(filename,"w")
    out.write("# "+label[0])
    for iLine in lines: out.write("\tZ("+label[1]+str(iLine)+")")
    out.write("\n")
    for i in range(x.size):
        out.write( (FMT) % (x[i]) )
        for zVal in z[:,i]:
            out.write( ("\t"+FMT) % (zVal) )
        out.write("\n")
    out.close()
    

def circularMask(w, h, center=None, rMax=None, rMin=None):
    if center is None: # use the middle of the image
        center = (int(w/2), int(h/2))
    if rMax is None: # use the smallest distance between the center and edges
        rMax = min(center[0], center[1], w-center[0], h-center[1])
    Y, X = np.ogrid[:h, :w]
    dist2 = (X - center[0])**2 + (Y-center[1])**2
    if rMin is None: mask = dist2 <= rMax**2
    else: mask = np.logical_and(dist2 <= rMax**2, dist2 >= rMin**2)
    return mask


def XY(shape, dtype=np.uint16):
    y = np.arange(shape[1],dtype=dtype)
    x = np.arange(shape[0],dtype=dtype)
    return np.tile(x,(shape[1],1)).transpose(), np.tile(y,(shape[0],1))

def center(array):
    x,y = XY(array.shape)
    norm = array.sum()
    xCOM = (x*array).sum()/norm
    yCOM = (y*array).sum()/norm
    return((xCOM,yCOM))

def untilt(array, avg=None):
    """ subtract the tilt from a 2D surface

    subtract the macroscopic tilt from np.ndarray <array>, assuming it to be 
    linear in both directions.

    returns: the untilted numpy.ndarray
    """
    nx,ny = array.shape
    nmin = min(nx,ny)
    if avg==None: avg = nmin//32
    if avg < 1: avg = 1
    X, Y = np.ogrid[:nx, :ny]
    tl = array[:avg,:avg].mean()
    tr = array[:avg,-avg:].mean()
    bl = array[-avg:,:avg].mean()
    br = array[-avg:,-avg:].mean()
    slopeX = (bl-tl + br-tr)/2
    slopeY = (tr-tl + br-bl)/2
    X = X*slopeX/nx
    Y = Y*slopeY/ny
    return(array - X - Y)


def psd(array, output=""):
    """ compute Power Spectral Density (PSD) of an array
    
    calculates the PSD of numpy.ndarray <array> with auto quasi-log q space.
    writes the PSD to file <output> if given.

    returns: numpy.ndarray with 2 columns containing q and C(q) values.
    """

    if WARN: print("[Warn:Dimens] Assuming dx=dy and ny<=nx.",flush=True)
    if WARN: print("[Warn:Period] psd() uses FFT, which assumes periodicity.",flush=True)
    arrayF = np.fft.rfft2(array)
    nxF = arrayF.shape[0]; nyF = arrayF.shape[1]
    idx = np.linspace(0, np.log(nyF-nyF//8), nyF//8)
    idx = np.exp(idx).astype(int) + np.arange(nyF//8)
    result = np.zeros((len(idx)-1,2))
    result[:,0] = np.pi/(YLIM[1]-YLIM[0]) * (idx[:-1] + idx[1:])
    
    abs2F = np.real(arrayF*np.conjugate(arrayF))
    for ir in range(0,len(idx)-1):
        rMax = idx[ir+1]; rMin = idx[ir]
        mask = circularMask(nyF, nxF, (0,0), rMax, rMin)
        mask = np.logical_or(mask, circularMask(nyF, nxF, (0,nxF), rMax, rMin))
        result[ir,1] = np.mean(abs2F[mask])

    # normalization: Tevis D B Jacobs et al 2017 Surf. Topogr.: Metrol. Prop. 5 013001
    area = (XLIM[1]-XLIM[0])*(YLIM[1]-YLIM[0])
    nxny = array.shape[0]*array.shape[1]
    result[:,1] = result[:,1]*area/(nxny**2) 

    if output: np.savetxt(output, result, fmt="%.5e", header="q\tPSD(2D-iso)(q)")
    return(result)


def plotPSD(array, axis=None):
    """ plot PSD of array
    
    plots the PSD of numpy.ndarray <array>.
    <axis>, the matplotlib.axes.Axes object on which to plot, can be specified.

    returns: matplotlib.axes.Axes object containing the plot
    """

    psdData = psd(array)
    if not axis:
        plt.figure()
        plt.loglog()
        plt.xlabel("q")
        plt.ylabel("PSD")
        axis = plt.gca()
    axis.plot(psdData[:,0], psdData[:,1])
    return(axis)


def createBorder(array, nxNew, nyNew):
    """ create border around surface 

    generates a (<nxNew>, <nyNew>) numpy.ndarray that contains the (nx, ny) 
    numpy.ndarray <array>.
    <nxNew> must be larger than nx and <nyNew> larger than ny.

    returns: the new numpy.ndarray
    """

    (nxOld,nyOld) = array.shape
    if WARN: print("[Warn:Border] createBorder() assumes array to be in equilPos format.",flush=True)
    if (nxNew < nxOld or nyNew < nyOld):
        print("[ERROR] surface larger than bordered surface!")
        return(array)
    nxStart = (nxNew-nxOld)/2
    if (nxStart != int(nxStart)) and WARN: print("[Warn:Border] nxNew-nxOld is not a multiple of 2.",flush=True)
    nxStart = int(nxStart)
    nyStart = (nyNew-nyOld)/2
    if (nyStart != int(nyStart)) and WARN: print("[Warn:Border] nyNew-nyOld is not a multiple of 2.",flush=True)
    nyStart = int(nyStart)
    nxEnd = nxStart + nxOld
    nyEnd = nyStart + nyOld
    result = array.max()*np.ones((nxNew,nyNew),dtype=float)
    result[nxStart:nxEnd,nyStart:nyEnd] = array
    return(result)





#%% ----- Other: mostly for legacy contMech versions ----- %%#
    
def logsmooth(infilename, outfilename="", nbin=100, usecols=(0,1)):
    if outfilename == "": outfilename = infilename + "-log"
    psd = np.loadtxt(infilename,usecols=usecols)
    fid = open(outfilename,"w")
    q = psd[:,0]
    # logarithmically spaced bin edges for equally spaced log plot
    if q[1] > 0: step = (1.001*q.max()/q[1])**(1./nbin)
    binedge = 0.999*q[1]*np.power(step, np.arange(nbin+1) )
    # linearly spaced bin edges
    #counts,binedge = np.histogram(q,nbin)
    qval = 0.5*(binedge[0:nbin] + binedge[1:(nbin+1)])
    spec = np.zeros((nbin,len(usecols)-1))
    for ibin in range(1,nbin+1):
        mask = np.logical_and(q>binedge[ibin-1],q<binedge[ibin])
        if mask.sum() > 0: 
            for iCol in range(1,len(usecols)):
                spec[ibin-1,iCol-1] = np.mean(psd[mask,iCol])
        else: spec[ibin-1,1:] = 0
        fid.write("%.5e" % qval[ibin-1])
        for iCol in range(1,len(usecols)): fid.write("\t%.5e" % spec[ibin-1,iCol-1])
        fid.write("\n")
    fid.close()
    #return(spec)


def contStats(filepath="", nbin=100):
    global SHOW
    if filepath != "":
        import os
        os.chdir(filepath)
    dist = readConfig("equilPos0.dat")
    dataElast = readConfig("elSheet1.dat")
    dist = dist-dataElast
    del dataElast
    # Criterion for being in contact
    incontact = dist <= 0
    dist[incontact] = 0
    contArea = incontact.mean(); meanGap = dist.mean()
    del dist
    # Pressure in contact histogram
    pressCont = readConfig("elSheet1.dat",3)
    pressCont = pressCont[incontact]
    pCounts,binedge = np.histogram(pressCont,nbin)
    pVal = 0.5*(binedge[0:nbin] + binedge[1:(nbin+1)])
    pMean = pressCont.mean(); pStd = pressCont.std()
    # Print
    print("pressure in cont.: %.5e +/- %.5e" % (pMean, pStd) )
    print("rel. contact area: %.4f" % contArea)
    print("mean gap:          %.5e" % meanGap)
    # Dump
    fid = open("params.out","a")
    fid.write("\n0\t# PYTHON post-analysis start\n")
    fid.write("%.5e\t\t# mean pressure in contact\n" % pMean)
    fid.write("%.5e\t\t# std of pressure in contact\n" % pStd)
    fid.write("%.4f\t\t# rel. contact area\n" % contArea)
    fid.write("%.5e\t\t# mean gap\n" % meanGap)
    fid.write("0\t# PYTHON post-analysis end\n")
    fid.close()
    # Plot
    show_orig = SHOW
    SHOW = False
    plotImg(incontact)
    plt.figure()
    plt.xlabel("pressure in contact")
    plt.ylabel("probability")
    plt.plot(pVal,pCounts/pressCont.size)
    SHOW = show_orig
    show()


def int1D(filename, xcol=0, ycol=1):
    array = np.loadtxt(filename)
    if array.shape[0] == 2:
        x = array[0,:]
        y = array[1,:]
    else:
        x = array[:,xcol]
        y = array[:,ycol]
    return(np.trapz(y,x))


def flatPunch(nx, ny, r, h):
    piston = circularMask(nx, ny,None,r)
    piston = piston.astype(float)
    piston = -h*piston + h
    if WARN: print("[Info:Format] faltPunch() returns array in equilPos format.",flush=True)
    return(piston)



#%% ----- Interaction parameters for legacy contMech versions ----- %%#

GAMMA = 50e-3   # surface energy in SI base units
E = 2.7e6       # contact modulus for PDMS in SI base units
def curvEla(dy): return(E * np.pi / (2*dy) ) # max elast. curv for 1D surface, qMax = pi/dy
def curvEla2D(dx,dy): return(E * np.pi * np.sqrt((1./dx)**2 + (1./dy)**2) / 2)
def curvMorse1(sigMorse): return(GAMMA/(sigMorse**2))
def curvMorse2(sigMorse): return(GAMMA/(4*sigMorse**2))
def curvCos1(sigCosine): return(GAMMA*(np.pi)**2/(4*sigCosine**2))
def curvCos2(sigCosine): return(GAMMA*(np.pi)**2/(2*sigCosine**2))
def uMorse1(dy,sigMorse): return(np.sqrt(curvMorse1(sigMorse)*dy/E)) # µ_rho for Morse1
def uMorse2(dy,sigMorse): return(np.sqrt(curvMorse2(sigMorse)*dy/E)) # µ_rho for Morse2
def uCos1(dy,sigCosine): return(np.sqrt(curvCos1(sigCosine)*dy/E)) # µ_rho for Cosine1
def uCos2(dy,sigCosine): return(np.sqrt(curvCos2(sigCosine)*dy/E)) # µ_rho for Cosine2
def sigMorse1(dy,uRho): return(np.sqrt(GAMMA*dy/E)/uRho)
def sigMorse2(dy,uRho): return(np.sqrt(GAMMA*dy/E)/(2*uRho))
def sigCos1(dy,uRho): return(np.sqrt(GAMMA*dy/(4*E))*np.pi/uRho)
def sigCos2(dy,uRho): return(np.sqrt(GAMMA*dy/(2*E))*np.pi/uRho)