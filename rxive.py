# ----- GLOBAL variables ----- #

# hard-coded
DEBUG = False
FOLDER = "frames"
ARCHIVE = "zip"

# defaults (overwritten bei sys.argv)
DEFAULT_PATH = "."
DEFAULT_EXTRACT = False



# ----- import modules ----- #

import sys
import os
import shutil

def relDirsFiles(fullpath):
  """ List the whole directory tree and all files contained in a directory

  The file list contains files in the parent and all subdirectories.
  The directorytree can be passed to mkdirtree().

  Parameters
  ----------
  fullpath : str
    path to the folder whose contents are listed

  Returns
  -------
  folderlist : list of str
    the whole directory tree contained in <fullpath>, relative to <fullpath>
  filelist : list of str
    relative paths to all files contained in the directroy tree

  """

  if fullpath[-1] != os.sep: fullpath += os.sep

  # base path
  folderlist = []; filelist = []
  lsdir = sorted(os.listdir(fullpath))
  for obj in lsdir:
    if os.path.isdir(fullpath+obj): folderlist.append(obj)
    else: filelist.append(obj)

  # recursive through sub-directories
  found = True
  newfolders = folderlist.copy()
  while found:
    toadd = []
    for folder in newfolders:
      folder = folder + os.sep
      lsdir = sorted(os.listdir(fullpath+folder))
      for obj in lsdir:
        if os.path.isdir(fullpath+folder+obj): 
          toadd.append(folder+obj)
        else: filelist.append(folder+obj)
    if len(toadd):
      folderlist = folderlist + toadd
      newfolders = toadd.copy()
    else: found = False

  return(folderlist, filelist)


def rxive(inpath=DEFAULT_PATH, extract=DEFAULT_EXTRACT):

  # ----- determine files/folders to operate on ----- #

  paths,files = relDirsFiles(inpath)
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



if __name__ == '__main__': 
  extract = DEFAULT_EXTRACT
  path = DEFAULT_PATH

  # parse sys.argv
  if len(sys.argv) > 1:
    for arg in sys.argv[1:]:
      if DEBUG: print("sys.arg", arg)
      if arg in ["-x", "-e"]: extract = True
      elif os.path.isdir(arg): path = arg
  
  rxive(path, extract)