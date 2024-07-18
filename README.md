cmutils
=======

Utilities for analyzing, plotting, animating and manipulating data from and for `contMech` simulations.

`contMech` is a continuum mechanical contact simulation program developed in the group of Prof. Dr. Martin Müser at Saarland University, Germany.
As of 2024, the official `contMech` repository is *not* public. 
However, the undocumented source code of a legacy version has been published [here](https://doi.org/10.5281/zenodo.12726265) in the context of a [PRL](https://doi.org/10.1103/PhysRevLett.131.156201).
The `cmutils` repository was made public to grant access to any `contMech` users who do not have a GitHub account.

Documentation can be found at https://sintharic.github.io/cmutils/, which is work in progress.

--------------------------------------------------------------------------------

Python usage and installation
-----------------------------
These Python modules require a scientific Python3 environment like e.g. Anaconda (https://www.anaconda.com/products/individual-d). 
Assuming the files are in your working directory or Anaconda/Python path (see Recommendations below), you can import them like any Python module (e.g. `import cmutils`)  

cmutils.py
----------
Miscellaneous functions for input, output, manipulation and plotting of `contMech` config files.
Very powerful in combinations with numpy's and scipy's standard manipulation and statistics tools.

Example to import, view, rotate and save a `contMech` config file:  
`import numpy as np`  
`import cmutils as cm`  
`dat0 = cm.readConfig("konfig0E.dat")`  
`cm.plotImg(dat0)`  
`cm.dumpConfig(np.rot90(dat0),"konfig0E.real")`  

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

Essentially, the previously mentioned `runCont2D()` is just a wrapper for `cmanim.init("Cont2D",paths="."); cmanim.run()`, where `init(...)` automatically determines plot parameters, which you can (but most of the time don't need to) overwrite.  

cmjobs.py
---------
Scans for currently running `contMech` processes and estimates time remaining based on nTime and line count of the gMoni.dat file.

Usage from terminal:  
`python3 cmjobs.py`  

It can also be used to list or terminate all processes whose working directory matches a certain pattern. For example, `python3 cmjobs.py kill Hertz`. The program will ask for confirmation before doing so.

stlutils.py
-----------
Converts numpy arrays and `contMech` configs to .stl files for CAD software and 3D printers.
Requires the 3rd party module numpy-stl (https://pypi.org/project/numpy-stl/).

Example to convert a blurred random array to stl:  
`import stlutils`  
`import numpy as np`  
`import scipy.ndimage as ndi`  
`data = np.random.default_rng(1234).uniform(0,1,(128,128))`  
`data = ndi.gaussian_filter(data,2)`  
`stlutils.convertArray(data,"random.stl")`  

cmparams.py
-----------
Reads all the parameters used in a `contMech` simulation, as indicated by its params file.

Example to find out whether sheet with ID=1 is elastic:  
`import cmparams as cp`  
`sim = cp.read("params.in")`  
`print(sim.sheet[1].nElast)`  

modparams.py
------------
Modifies `contMech` params files by adding, multiplying or setting values of parameters.

Example to reduce the time step size in all subfolders of the current directory:  
`import modparams as mp`  
`import os`  
`files = [folder+'/params.in' for folder in os.listdir('.') if os.path.isfile(folder+'/params.in')]`  
`for file in files: mp.multiply(['dTime', 'nTime'], [1./2, 2], file)`  

Recommendations
---------------
To permanently add these modules to your Anaconda path, you have the options described here: https://stackoverflow.com/questions/37006114/anaconda-permanently-include-external-packages-like-in-pythonpath.
My recommendation is `conda develop ~/git/cmutils`.

If you are often using modules in the IPython console, also consider importing them by default. 
Any .py file in the folder `~/.ipython/profile_default/startup` will automatically be loaded at IPython startup.

Furthermore, you can define aliases in your .bashrc or .zshrc file, e.g. `alias cmjobs='python3 ~/git/cmutils/cmjobs.py'`.

--------------------------------------------------------------------------------

C++ usage and installation
--------------------------
See list of programs below and compile them with the respective make command.
Some of these require a working FFTW3 installation. 
A very brief instruction on how to install FFTW on Windows is written in the Makefile.
The Makefile automatically checks the most common FFTW installation paths. 
If you encounter problems, specify your individual installation paths manually.

Apart from GF2D.exe and GF3D.exe, these are stand-alone applications, which do not require you to edit the source code. 
Instead, they take command line arguments, including the usual `--help` (or `-h`) option, which should explain everything you need to know.

GF2D and GF3D, on the other hand, do not take arguments and can only be edited in the source code. They are most helpful to create reference data once, which can then be compared to simulation results.

Green's functions (GF2D.exe and GF3D.exe):
------------------------------------------
Compile with `make green`.

Computes 2D/3D displacement fields from the Fourier coefficients of Green's functions, as given e.g. by Menga, N. IJSOLSTR 164 (2019) 212–220.

convertImg (convert.exe):
-------------------------
Compile with `make convert`.

Converts a 3D image file, as generated by the INM's confocal Mahr microscope, to a surface file, as readable e.g. by the https://contact.engineering website or by `cmutils.py`'s `readImg()`.

extractRange (range.exe):
-------------------------
Compile with `make range`.

Extracts a given data range from a surface config file, as e.g. created by the `contMech` simulation or `cmutils.py`'s `dumpConfig()`. The range is given in true x- and y-coordinates.

statsSurf (surf.exe):
---------------------
Compile with `make surf`.

Reads a surface file, as e.g. created by `convert.exe` or `cmutils.py`'s `dumpImg()`, computes statistics (RMS height/gradient/curvature) in Fourier space and dumps its PSD to a file.

statsParams (params.exe):
-------------------------
Compile with `make params`.

Outputs the exact RMS height/gradient/curvature stats that would belong to a given set of `contMech`'s "addSelfAffine" parameters: hurst, lambdaR, lambdaS, ...

gapconvert (gapconvert.exe):
-------------------------
Compile with `make gap`.

Reads a config file containing the local gap between two surfaces, converts the data to the format necessary for Reynolds flow simulations, and saves it to a file with the correct name.

Shell commands:
---------------
Install with `make install`.

Creates shell scripts (batch files) and installs them in a directory (that is likely to be) in the user's PATH. The corresponding commands will have the same names as above, except without the `.exe` and with a leading `cm_`. This way, e.g. `params.exe` can be accessed from anywhere via the command `cm_params`.

Windows users should check their environment variables.