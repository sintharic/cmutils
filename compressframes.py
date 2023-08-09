# ----- GLOBAL variables ----- #

# hard-coded
DEBUG = False
FOLDER = "frames"
ARCHIVE = "zip"

# defaults (overwritten bei sys.argv)
PATH = "."
EXTRACT = False



# ----- import modules ----- #

import sys
import os
sys.path.append("/Users/christian/git/fsync")
import shutil
import futil

# parse sys.argv
if len(sys.argv) > 1:
  for arg in sys.argv[1:]:
    if DEBUG: print("sys.arg", arg)
    if arg in ["-x", "-e"]: EXTRACT = True
    elif os.path.isdir(arg): PATH = arg


def main(inpath=PATH, extract=EXTRACT):

  # ----- determine files/folders to operate on ----- #

  paths,files = futil.relDirsFiles(inpath)
  if extract:
    operation = "extracting"
    tail = f"{FOLDER}.{ARCHIVE}"
    paths = [inpath+os.sep+file for file in files if os.path.split(file)[-1]==tail]
  else: 
    operation = "archiving"
    tail = f"{FOLDER}"
    paths = [inpath+os.sep+path for path in paths if os.path.split(path)[-1]==tail]



  # ----- loop through all files/folders ----- #

  nPaths = len(paths)
  if nPaths==0: 
    print("Nothing to do.")
    return
  wd = os.getcwd()
  for i,path in enumerate(paths):
    print(f"{operation} {path} ({i+1}/{nPaths})")
    os.chdir(os.path.split(path)[0])

    # CASE 1: print debug message
    if DEBUG:
      print(os.getcwd())#DEBUG

    # CASE 2: extract archive, remove folder
    elif extract:
      try: shutil.unpack_archive(f"{FOLDER}.{ARCHIVE}", FOLDER)
      except: 
        print("ERROR unpacking '%s' in %s" % (f"{FOLDER}.{ARCHIVE}", os.getcwd()))
        os.chdir(wd)
        continue

      try: os.remove(f"{FOLDER}.{ARCHIVE}")
      except: 
        print("ERROR removing '%s' in %s" % (f"{FOLDER}.{ARCHIVE}", os.getcwd()))
        os.chdir(wd)
        continue

    # CASE 3: create archive, remove folder
    else:
      try: shutil.make_archive(FOLDER, ARCHIVE, FOLDER)
      except: 
        print("ERROR creating '%s' in %s" % (f"{FOLDER}.{ARCHIVE}", os.getcwd()))
        os.chdir(wd)
        continue

      try: shutil.rmtree(FOLDER)
      except: 
        print("ERROR removing '%s' in %s" % (FOLDER, os.getcwd()))
        os.chdir(wd)
        continue

    os.chdir(wd)



if __name__ == '__main__': main()