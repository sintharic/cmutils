#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-----------------------------------------
  Modification of contMech params files
-----------------------------------------

@author: thescientist
"""

DEBUG = False
SED = 'sed'

import os, sys
from cmparams import paramType
from subprocess import run as call



# check the installed version of sed/gsed
testfile = "_modparams_test_"
with open(testfile, "w") as fid: fid.write("original")
try: 
  proc = call(["gsed", "-i", r"/original/c\modified", testfile], capture_output=True)
  returnval = proc.returncode
  if returnval: print(proc.stderr)
except: 
  returnval = 1
if returnval:
  try: 
    proc = call(["sed", "-i", r"/original/c\modified", testfile], capture_output=True)
    returnval = proc.returncode
    if returnval: print(proc.stderr)
  except: 
    returnval = 1
  if returnval:
    os.remove(testfile)
    raise OSError("sed/gsed not working as intended.")
  else: 
    with open(testfile, "r") as fid: content = fid.read()
    if content == "modified\n": SED = "sed"
    else: 
      if DEBUG: print(content)
      os.remove(testfile)
      raise OSError("sed/gsed not working as intended.")
else: 
  with open(testfile, "r") as fid: content = fid.read()
  if content == "modified\n": SED = "gsed"
  else: 
    if DEBUG: print(content)
    os.remove(testfile)
    raise OSError("sed/gsed not working as intended.")
os.remove(testfile)


def sed(command, file):
  proc = call([SED, "-i", command, file], capture_output=True)
  if proc.returncode:
    print(proc.stderr)
    raise proc.check_returncode()


def check(params, values=None):
  if isinstance(params, str): params = [params]
  for i in range(len(params)):
    if params[i][:2] == '# ': params[i] = params[i][2:]
    if params[i][-2:] == ' #': params[i] = params[i][:-2]

  if values is not None:
    if isinstance(values, int) or isinstance(values, float): values = [values]
    if len(values)==1: values = [values[0]]*len(params)
    if len(values)!=len(params):
      raise ValueError(f'parameter names and values must have same length, but have {len(params)} and {len(values)}')
    return params, values
  else: 
    return params


def param_type(param):
  try: ptype = paramType[param]
  except: 
    print(f'WARNING: Unknown paramType for _{param}_. Assuming float.')
    ptype = 'float'
  return ptype


def get_values(params, file):
  # sanity checks
  if not os.path.isfile(file): 
    raise ValueError('invalid params file: %s' % file)
  params = check(params)

  values = [0]*len(params)
  with open(file,'r') as fid:
    for line in fid.readlines():
      for i,param in enumerate(params):
        ptype = param_type(param)
        if f'# {param} #' in line:
          #print(line)#DEBUG
          if ptype == 'int': values[i] = int(line.split()[0])
          else: values[i] = float(line.split()[0])

  return values


def set(params, values, file='params.in'):
  # sanity checks
  if not os.path.isfile(file): 
    raise ValueError('invalid params file: %s' % file)
  params, values = check(params, values)

  for param,pval in zip(params,values):
    ptype = param_type(param)
    if ptype == 'int': 
      newval = int(round(pval))
      command = f"/# {param} #/c\\{newval}\\t\\t# {param} #"
    else: 
      newval = pval
      command = f"/# {param} #/c\\{newval:g}\\t\\t# {param} #"

    if DEBUG: print(f"{SED} -i '{command}' {file}")
    
    sed(command, file)


def add(params, values, after, file='params.in'):
  # sanity checks
  if not os.path.isfile(file): 
    raise ValueError('invalid params file: %s' % file)
  params, values = check(params, values)
  if after[:2] == '# ': after = after[2:]
  if after[-2:] == ' #': after = after[:-2]

  for param,pval in zip(params,values):
    ptype = param_type(param)
    if ptype == 'int': 
      newval = int(round(pval))
      command = f"/# {after} #/a\\{newval}\\t\\t# {param} #"
    else: 
      newval = pval
      command = f"/# {after} #/a\\{newval:g}\\t\\t# {param} #"

    if DEBUG: print(f"{SED} -i '{command}' {file}")
    
    sed(command, file)


def multiply(params, factors, file='params.in'):
  # sanity checks
  params, factors = check(params, factors)

  old_values = get_values(params, file)
  new_values = [val*fac for val,fac in zip(old_values, factors)]
  set(params, new_values, file)


# old functions for compatibility
paramNames = ['nTime']
paramFacs  = [1]
def modify(file):
  multiply(file, paramNames, paramFacs)
def run(file=None):
  if len(paramNames) != len(paramFacs): 
    raise ValueError('paramNames, paramFacs and paramTypes must have the same length!')
  modify(file)