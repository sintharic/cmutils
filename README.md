cmutils
=======

Utilities for analyzing, plotting, animating and manipulating data from contMech simulations.

Usage
-----
These Python modules require a scientific Python3 environment like e.g. Anaconda (https://www.anaconda.com/products/individual-d). 
Assuming the files are in your working directory or Anaconda/Python path, just import them like any other Python module via  
`import cmutils`  
or  
`import cmanim`

cmutils.py
----------
Miscellaneous functions for input, output, manipulation and plotting.
Very powerful in combinations with numpy's and scipy's standard manipulation and statistics tools.

cmanim.py
---------
An easy-to-use library for animating simulation results.
The focus is on fully automatic usage, where the best possible parameters are read from params.out. 

Example:  
`import cmanim`  
`animation1 = cmanim.runDispX()`  
`cmanim.save(animation1, "movie.gif")`  
`animation2 = cmanim.runCont2D()`

Specifying custom parameters and animating more specific things is possible but currently rather inconvenient.

mpl.py
------
contains presetes etc. for plotting with matplotlib. It is really useful to ensure consistent visuals in publications and working with matplotlib's rcParams in general, giving access to very advanced functionality.

Other recommendations
---------------------
To permanently add these modules to your Anaconda path, you have the options described here: https://stackoverflow.com/questions/37006114/anaconda-permanently-include-external-packages-like-in-pythonpath

If you are often using modules in the IPython console, consider importing them by default. 
Any .py file in the folder `~/.ipython/profile_default/startup` will automatically be loaded at IPython startup.
