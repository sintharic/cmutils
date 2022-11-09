#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------
  Modification of contMech params files
-----------------------------------------

@author: thescientist
"""

SED = 'gsed'

import os, sys
import cmparams as cp

paramNames = ['vzConstCOM', 'dTime']
paramFacs  = [3./4.75, 4.75/3]


def checkParams():
  if len(paramNames) != len(paramFacs): 
    sys.exit('paramNames, paramFacs and paramTypes must have the same length!')


def modify(file):
  # determine file name
  if not file:
    if len(sys.argv) > 1: file = sys.argv[1]
    elif os.path.isfile('params.in'): file = 'params.in'
    else: sys.exit('no params file name was given.')
  elif not os.path.isfile(file): 
    sys.exit('invalid params file name.')

  # read old param values
  paramValue = [0]*len(paramNames)
  with open(file,'r') as fid:
    for line in fid.readlines():
      for i,param in enumerate(paramNames):
        if param[0] == '#': param = param[2:]
        if param[-1] == '#': param = param[:-2]
        try: ptype = cp.paramType[param]
        except: 
          print(f'WARNING: invalid paramType for _{param}_')
          ptype = 'float'
        param = '# ' + param + ' #'
        if param in line:
          #print(line)#DEBUG
          if ptype == 'int': paramValue[i] = int(line.split()[0])
          else: paramValue[i] =  float(line.split()[0])

  # calculate and replace with new param values
  for param,pval,pfac in zip(paramNames,paramValue,paramFacs):
    if param[0] == '#': param = param[2:]
    if param[-1] == '#': param = param[:-2]
    try: ptype = cp.paramType[param]
    except: 
      print(f'WARNING: invalid paramType for _{param}_')
      ptype = 'float'
    param = '# ' + param + ' #'
    if ptype == 'int': 
      newval = int(round(pfac*pval))
      command = f"{SED} -i '/{param}/c\\{newval}\\t\\t{param}' {file}"
    else: 
      newval = pfac*pval
      command = f"{SED} -i '/{param}/c\\{newval:g}\\t\\t{param}' {file}"

    print(command)#DEBUG
    os.system(command)


def run(file=None):
  checkParams()
  modify(file)


if __name__=='__main__': run()