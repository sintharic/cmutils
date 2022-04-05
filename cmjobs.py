#user = "christian"
executable = "contMech"

DEBUG = False
WARN = False

import time, os
import datetime
import psutil as ps
import getpass

user = getpass.getuser()
if DEBUG: print("USERNAME:", user)

nFound = 0
t0 = time.time()

print("-----------------------------------------------")
print(" ProcID | prog | remaining | working directory ")
print("-----------------------------------------------")

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
  
  # proceed only for contMech processes
  if executable in pname:
    
    # gather info from params file
    nTime = 0; nx = 0; ny = 0
    paramsfile = open(pcwd+"/params.out","r")
    for line in paramsfile:
      if "# nxGlobal" in line: nx = int(line.split("#")[0])
      elif "# nyGlobal" in line: ny = int(line.split("#")[0])
      elif ("# nTime" in line): 
        if not (("# nTimeOn" in line) or ("# nTimeOff" in line)): 
          nTime = int(line.split("#")[0])
    if ny == 0: ny = nx
    paramsfile.close()

    # count lines of moni file
    readMoni = True
    iTime = 0
    try: monifile = open(pcwd+"/moni1-"+str(ny).zfill(4)+".dat","r")
    except IOError:
      try: monifile = open(pcwd+"/moni0-"+str(ny).zfill(4)+".dat","r")
      except:
        if WARN: print("[Warn:Moni] No moni file found for", pcwd, ny)
        readMoni = False
    if readMoni: 
      for iTime,_ in enumerate(monifile): pass
      # since iTime starts at 0, the header is automatically neglected
      monifile.close()

    # gather time started
    if iTime <= 3: tRemain = " unknown "
    else:
      tRunning = t0 - os.path.getmtime(pcwd+"/params.out")
      tRemain = int(tRunning*(nTime-iTime)/iTime)
      if tRemain < 0: tRemain = " unknown "
      elif tRemain > 24*60*60: 
        tRemain = str(round(tRemain/(24*60*60),1))+" days"
      else:
        tRemain = str(datetime.timedelta(seconds=tRemain))
      tRemain = tRemain.rjust(9)
      
    nFound += 1
    num = str(pPID).rjust(7)
    #"%7i / %7i"%(iTime,nTime)
    print(num, "| %.2f"%(iTime/nTime), "|", tRemain, "|", pcwd)

print(nFound, executable, "process(es) running.")
