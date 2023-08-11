#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------
  Overview of running contMech Jobs
-------------------------------------

@author: thescientist
"""

executable = "contMech"

DEBUG = False
WARN = False

import time, os, sys
import datetime
import psutil as ps
import getpass

user = getpass.getuser()
now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M")
if DEBUG: print("USERNAME:", user)

nFound = 0
wFire = False
t0 = time.time()

# command line options
do_filt = False
function = ""
pattern = None
if len(sys.argv)>1: 
  function = sys.argv[1]
  if function not in ["list", "kill", "finish"]:
    print("ERROR: Unknown option: %s" % function)
    print("Valid options are 'list', 'finish' and 'kill'.")
    sys.exit()
  else:
    do_filt = True
if len(sys.argv)>2: 
  pattern = sys.argv[2]
  if pattern in ["-a", "--all", "*"]: 
    do_filt = False
    pattern = None



def find_processes(exe):
  # find processes with a given name

  result = []
  for proc in ps.process_iter():
    
    # gather basic process info
    pusr = proc.username()
    if pusr != user: continue;
    pname = proc.name()
    try: pcwd = proc.cwd()
    except (ps.AccessDenied, ps.ZombieProcess) as err:
      if DEBUG: print(pname, ": AccessDenied or ZombieProcess.")
      continue;
    pPID = proc.pid
    
    # filter contMech processes
    if exe in pname:
      result.append(proc)

  return result




print("-----------------------------------------------")
print(" ProcID | prog | remaining | working directory ")
print("-----------------------------------------------")

processes = find_processes(executable)
if do_filt and (pattern is not None):
  processes = [proc for proc in processes if pattern in proc.cwd()]

for proc in processes:
  pcwd = proc.cwd()
  pPID = proc.pid

  # gather info from params file
  nTime = 0; nx = 0; ny = 0
  fFire = 0; fLogMeasure = 0; dTime = 0;
  paramsfile = open(pcwd+"/params.out","r")
  for line in paramsfile:
    if "# nxGlobal" in line: nx = int(line.split("#")[0])
    elif "# nyGlobal" in line: ny = int(line.split("#")[0])
    elif ("# nTime" in line): 
      if not (("# nTimeOn" in line) or ("# nTimeOff" in line)): 
        nTime = int(line.split("#")[0])
    elif "# fFire" in line: fFire = int(line.split("#")[0])
    elif "# dTime" in line: dTime = float(line.split("#")[0])
    elif "# fLogMeasure" in line: fLogMeasure = int(line.split("#")[0])
  if ny == 0: ny = nx
  paramsfile.close()

  # count lines of moni file
  iTime = 0
  if not fLogMeasure:
    try: 
      with open(pcwd+"/gMoni.dat","r") as monifile:
        for iTime,_ in enumerate(monifile): pass
        iTime = iTime - 1
        # since iTime starts at 0, the header is automatically neglected, however time step 0 is not
    except FileNotFoundError: iTime = -nTime
    except BaseException as e: raise e
  else:
    if not fFire:
      try:
        with open(pcwd+"/gMoni.dat","rb") as monifile:
          monifile.seek(-2, 2) # 2nd arg = 2 --> relative to EOF
          while monifile.read(1) == b"\n": monifile.seek(-2, 1) # skip trailing empty lines
          monifile.seek(-1, 1) # 2nd arg = 1 --> relative to current position
          while monifile.read(1) != b"\n": monifile.seek(-2, 1)
          lineN = monifile.read().decode('UTF-8')
          idx = lineN.find("\t")
          time = float(lineN[:idx])
          iTime = int(round(time/dTime))
      except FileNotFoundError: iTime = -nTime
      except BaseException as e: raise e
    else: iTime = -nTime

  # determine time started and time remaining
  if iTime <= 3: tRemain = "unknown"
  else:
    tRunning = t0 - os.path.getmtime(pcwd+"/params.out")
    tRemain = int(tRunning*(nTime-iTime)/iTime)
    if tRemain < 0: tRemain = "unknown"
    elif tRemain > 24*60*60: 
      tRemain = str(round(tRemain/(24*60*60),1))+" days"
    else:
      tRemain = str(datetime.timedelta(seconds=tRemain))
  if fFire: 
    tRemain = "*" + tRemain
    wFire = True
  tRemain = tRemain.rjust(9)
    
  nFound += 1
  num = str(pPID).rjust(7)
  #"%7i / %7i"%(iTime,nTime)x

  # shorten directory path
  pcwd = pcwd.replace("/home/"+user, "~")
  pcwd = pcwd.replace("/Users/"+user, "~")
  print(num, "| %.2f"%(iTime/nTime), "|", tRemain, "|", pcwd)



# display total number of processes
if do_filt:
  print(nFound, executable, "process(es) matching pattern '%s'."%pattern)
else: 
  print(f'{nFound} {executable} process(es) running ({now}).')
if wFire: print("* simulations using FIRE may finish early.")



# kill processes matching the pattern
if function=="finish":
  ans = input("Are you sure you want to finish these processes? (y/n) ")
  if ans=="y":
    for proc in processes: 
      with open(f"{proc.cwd()}/finish", "w") as file:
        file.write(f"requested {now}\n")




# kill processes matching the pattern
if function=="kill":
  ans = input("Are you sure you want to kill these processes? (y/n) ")
  if ans=="y":
    for proc in processes: proc.kill()
