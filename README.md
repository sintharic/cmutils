cmutils
=======

Utilities for analyzing, plotting, animating and manipulating data from contMech simulations.

--------------------------------------------------------------------------------

Python Usage
------------
These Python modules require a scientific Python3 environment like e.g. Anaconda (https://www.anaconda.com/products/individual-d). 
Assuming the files are in your working directory or Anaconda/Python path, just import them like any other Python module via  
`import cmutils`  
or  
`import cmanim`

cmutils.py
----------
Miscellaneous functions for input, output, manipulation and plotting of config files.
Very powerful in combinations with numpy's and scipy's standard manipulation and statistics tools.

cmanim.py
---------
An easy-to-use library for animating simulation results.
The focus is on fully automatic usage, where the best possible parameters are read from params.out. 

Example to create and save a contact movie:  
`import cmanim`  
`animation = cmanim.runCont2D()`  
`cmanim.save(animation, "movie.gif")`

Specifying custom parameters is possible by calling individual initialization functions one after the other. 
E.g. you can animate the cross-sections along x-direction from two different simulations in one figure, limit it to frames 40 to 100, and zoom in to a certain xrange:  
`cmanim.init("DispX", [path1, path2] )`  
`cmanim.XLIM = [xmin, xmax]`  
`cmanim.START = 40`  
`cmanim.END = 100`  
`cmanim.NFRAMES = 60`  
`animation = cmanim.run()`  

Essentially, the previously mentioned `runCont2D()`is just a wrapper for `cmanim.init("Cont2D",paths="."); cmanim.run()`, where `init(...)` automatically determines plot parameters, which you can (but most of the time don't need to) overwrite.  

cmjobs.py
---------
Scans for currently running contMech processes and estimates time remaining based on nTime and length of moni file.

stlutils.py
-----------
Conversion of Python arrays and contMech configs to .stl files for 3D printing.
Requires the 3rd party module numpy-stl (https://pypi.org/project/numpy-stl/).

mpl.py
------
Contains presets etc. for plotting with matplotlib. It is really useful to ensure consistent visuals in publications and working with matplotlib's rcParams in general, giving access to very advanced functionality.

Recommendations
---------------
To permanently add these modules to your Anaconda path, you have the options described here: https://stackoverflow.com/questions/37006114/anaconda-permanently-include-external-packages-like-in-pythonpath

If you are often using modules in the IPython console, consider importing them by default. 
Any .py file in the folder `~/.ipython/profile_default/startup` will automatically be loaded at IPython startup.

Furthermore, you can define aliases in your .bashrc file, e.g. `alias cmjobs='python3 ~/git/cmutils/cmjobs.py'`.

--------------------------------------------------------------------------------

C++ usage
---------
Most of the contained C++ programs require a working FFTW3 installation. 
A very brief instruction of how to do that is written in the Makefile.
The Makefile also automatically checks many common FFTW installation paths to find it. 
If you encounter problems, specify those installation paths manually.

Apart from GF2D.exe and GF3D.exe, these are stand-alone applications, which do not require you to mess with the source code anymore. 
Instead, they take command line arguments, including a nice `-h` or `--help` option, which should explain everything you need.

GF2D and GF3D, on the other hand, must be edited in the source code to produce different results. 
This approach was chosen because there is a wide range of possible Green's functions and parameters that might be of interest, but usually only in very situational cases. 
It is most helpful to create reference data once, which can then be compared to simulation results.

Green's functions (GF2D.exe and GF3D.exe)
-----------------------------------------
These compute 2D/3D displacement fields from the Fourier coefficients of Green's functions, as given e.g. by Menga, N. IJSOLSTR 164 (2019) 212–220.

convertImg (convert.exe)
------------------------
Converts a 3D image file, as generated by the INM's confocal Mahr microscope, to a surface file, as readable e.g. by the https://contact.engineering website or by `cmutils.py`'s `readImg()`.

extractRange (range.exe)
------------------------
Extracts a given data range from a surface config file, as e.g. created by the contMech simulation or `cmutils.py`'s `dumpConfig()`. The range is given in true x- and y-coordinates.

statsSurf (surf.exe)
--------------------
Reads a surface file, as e.g. created by `convert.exe` or `cmutils.py`'s `dumpImg()`, computes statistics (RMS height/gradient/curvature) in Fourier space and dumps its PSD to a file.

statsParams (params.exe)
------------------------
Outputs the exact RMS height/gradient/curvature stats that would belong to a given set of contMech's "addSelfAffine" parameters: hurst, lambdaR, lambdaS, ...
