cmutils
=======

Utilities for analyzing, plotting, animating and manipulating data from contMech simulations.

Usage
-----
These Python modules require a scientific Python3 environment like e.g. Anaconda (https://www.anaconda.com/products/individual-d). 
Assuming the files are in your working directory or Anaconda/Python path, just import them like any other Python module via:
import cmutils
import cmanim

cmutils.py
----------
Miscellaneous functions for input, output, manipulation and plotting.
Very powerful in combinations with numpy's and scipy's standard manipulation and statistics tools.

cmanim.py
---------
An easy-to-use library for animating simulation results.
The focus is on fully automatic usage, where the best possible parameters are read from params.out. 

Example:
import cmanim
animation1 = cmanim.runDispX()
cmanim.save(animation1, "movie.gif")
animation2 = cmanim.runCont2D()

Specifying custom parameters and animating more specific things is possible but currently rather inconvenient.
