#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------
  Animation of contMech Simulations
-------------------------------------

An easy-to-use library for animating simulation results.
The focus is on fully automatic usage, where the best possible parameters are read from params.out. 


------------
  Examples
------------

Example to create and save a contact movie:  

.. code-block:: python

    import cmanim

    animation = cmanim.runCont2D()
    cmanim.save(animation, "movie.gif")

Specifying custom parameters is possible by calling individual initialization functions one after the other. 
E.g. you can animate the cross-sections along x-direction from two different simulations in one figure, limit it to frames 40 to 100, and zoom in to a certain xrange:  

.. code-block:: python

    cmanim.init("DispX", [path1, path2] )
    cmanim.XLIM = [xmin, xmax]
    cmanim.START = 40
    cmanim.END = 100
    cmanim.NFRAMES = 60
    animation = cmanim.run()

Essentially, the previously mentioned `runCont2D()` is just a wrapper for `cmanim.init("Cont2D",paths="."); cmanim.run()`, where `init(...)` automatically determines plot parameters, which you can (but most of the time don't need to) overwrite.  


---------------------
  API documentation
---------------------
"""



#%% ----- Global settings ----- %%#

# For testing
testpath = "/home/chris/src8fix/run2"

# Input
SHOW = True
FPS = 5
DPI = 180
LINESTYLE = "-"
MARKERSTYLE = ""
MODE = "DispX"
VERSION = "new" # "new" or "dev"
UNITS = ["arb. u.","arb. u."]

# Data parameters
NX = 0; NY = 0;
XLIM = [0,0]; YLIM = [0,0]
ZLIM = [+1.0e+38,-1.0e+38]
RESMOVIE = [0,0]
XRAMP = None 
YRAMP = None 
INCREMENT = 1
NFRAMES = 1
START = 0
END = -1
SHEET = []
INTER = []
LABELS = []
RAMP = False
GRID = False

IID = 0 # TODO: now always plots inter class with ID 0. Generalize this.

fATTRACT = True # TODO: read this in?

# Add another arbitrary ramp file to be plotted
ADDRAMP = {"plot"   : False,
           "x"      : [],
           "y"      : [],
           "ls"     : "k--",
           "lab"    : [] }




#%% ----- Class and function definitions ----- %%#

import numpy as np
import os
import sys
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import cmutils as cm


class gfmdSheet:
    # Defaults: should match contMech code
    ID = 0
    nx = 0; ny = 0;
    lengthX = 0; lengthY = 0;
    fRough = 0; nElast = 0; fKelvinVoigt = 0
    stiffness0 = 0
    f3dMovie = 0
    # compatibility parameters to join multiple paths
    nTime = 0; nFrames = 0; increment = 1
    # Ramp / vz
    vzConstCOM = 0; dTime = 0
    fSteppedRamp = 0; dzRamp = 0; rampPeriod = 0; nVeloTurnStep = -1
    rampName = ""
    ramp = np.zeros((1,1))

    konfigName = "konfig"
    files = [] # file names
    data = np.zeros((1,1))
    minmaxZ = []
    usecols = [0,1]
    
    def __init__(self, newID): 
        global INCREMENT
        self.ID = newID
        self.increment = INCREMENT
    
    def updateNxNy(self):
        global NX, NY
        if self.nx == 0: self.nx = NX
        if self.ny == 0: self.ny = NY

    def updateFiles(self, path):#, numFrames):
        """ sets up list of frame files belonging to this object """

        global NFRAMES
        
        # konfigName
        if MODE[-2:] == "3D": ending = ".dat"
        else: ending = ".datH"
        if VERSION == "new":
            self.konfigName = "konfig" + str(self.ID)
            if(self.fRough > 0): self.konfigName += "E"
            if(self.nElast > 0): self.konfigName = "frames/" + self.konfigName + "D"
        elif VERSION == "dev":
            if(self.nElast > 0): self.konfigName = "movie/elSheet" + str(self.ID)
            else: self.konfigName = "equilPos" + str(self.ID)
        
        # Individual profile file names
        self.files = []
        if self.nElast == 0:
            file = path + self.konfigName + ending
            if (not os.path.isfile(file)):
                print("[WARNING]", file, "does not exist.",flush=True)
            else: self.files.append(file)
        else:
            iFrame = 0
            while True:
                iFrame += 1
                file = path + self.konfigName + "." + str(iFrame) + ending
                if (not os.path.isfile(file)): break
                self.files.append(file)
            if iFrame < self.nFrames: 
                print("[WARNING]", file, "does not exist.",flush=True)
                self.nTime = int(self.nTime*(iFrame-1)/self.nFrames)
                self.nFrames = iFrame - 1
        
        # Ramp file name
        if self.fSteppedRamp:
            file = path + "ramp" + str(self.ID) + "-" + str(self.ny).zfill(4) + ".dat"
            if os.path.isfile(file):
                self.rampName = file
                #print("ramp period:", self.rampPeriod)
            else: print("[WARNING]", file, "not found.")
        elif self.vzConstCOM != 0:
            file = path + "moni" + str(self.ID) + "-" + str(self.ny).zfill(4) + ".dat"
            if os.path.isfile(file):
                self.rampName = file
                #TODO: check if working
            else: print("[WARNING]", file, "not found.")

        # Sanity check
        if len(self.files) == 0:
            sys.exit("[ERROR] No files found for sheet "+str(self.ID))
    
    def updateCols(self):
        """ determines which columns in the data files must be imported """

        if (MODE == "DispX") or (self.nx==1 and (MODE == "DispY")): 
            self.usecols = [0,1]
        elif MODE == "DispY": self.usecols = [0,2]
        elif MODE == "DispD": 
            if self.nElast: self.usecols = [5,6]
            else: self.usecols = [3,4]
        elif self.nElast>0 and MODE == "PressX": self.usecols = [0,3]
        elif self.nElast>0 and self.nx==1 and (MODE in ["PressX","PressY"]): 
            self.usecols = [0,2]
        elif self.nElast>0 and MODE == "PressY": self.usecols = [0,4]
        elif self.nElast>0 and MODE == "PressD": self.usecols = [5,7]
        elif (MODE=="Dist3D") or (MODE=="Cont3D"): self.usecols = [2]
        elif (MODE=="Disp3D"): self.usecols = [2]
        elif (MODE=="Press3D"):
            if self.nElast>0: self.usecols = [3]
            else: self.usecols = [2]
        elif self.nElast == 0: self.usecols = [0,1]
        else: print("[WARNING] sheet["+str(self.ID)+"] cannot update cols according to MODE.")
        
    def updateData(self,numFrame):
        """ reads and sets up the data file needed to plot the next frame """

        idx = START + self.increment*(numFrame-1-START) # numFrame starts at 1, not 0!
        
        if idx > self.nFrames:
            print("[WARNING] no further frames found for sheet")
            return()

        if MODE[-2:] == "3D":
            #import cmutils as cm
            cm.WARN = False

            if (self.nElast > 0):
                self.data = cm.readConfig(self.files[idx], usecols=self.usecols)
            elif (self.data.size == 1):
                # rigid sheet only updates once
                self.data = cm.readConfig(self.files[0], usecols=self.usecols)
            # TODO: What if elast sheet was downsampled to be smaller than 512x512?
            if (self.nx!=NX) or (self.ny!=NY) or (self.nx>512) or (self.ny>512): 
                self.data = cm.resample(self.data,(NX,NY))
            #WIP
            if (self.data.shape[0]!=len(RESMOVIE[0])):
                self.data = self.data[RESMOVIE[0],:]
                self.data = self.data[:,RESMOVIE[1]]
                #print(numFrame,"elast sheet",self.ID,"after: ", self.data.shape)#DEBUG

        else: 
            if (self.nElast > 0):
                try: self.data = np.loadtxt(self.files[idx], usecols=self.usecols)
                except: print("[WARNING] Konfig file #"+str(numFrame)+" not found")
                if UNITS[0][0]=="µ":
                    self.data = 1e6*self.data
                if UNITS[0]=="mm":
                    self.data = 1e3*self.data
            elif (self.data.size == 1):
                # rigid sheet only updates once
                self.data = np.loadtxt(self.files[0], usecols=self.usecols)
                if UNITS[0][0]=="µ":
                    self.data = 1e6*self.data
                if UNITS[0]=="mm":
                    self.data = 1e3*self.data
                
    
    def updateDims(self,path):
        global NX, NY
        # nx, ny
        if self.nx == 0: self.nx = NX
        if self.ny == 0: self.ny = NY
        # minmaxZ
        if (self.nElast == 0):
            if self.data.size == 1: self.updateData(0)
            self.minmaxZ = [self.data[:,1].min(), self.data[:,1].max()]
        elif (self.nElast > 0):
            monipath = path + "moni"+str(self.ID)+"-"+str(self.ny).zfill(4)+".dat"
            if (not os.path.isfile(monipath)):
                sys.exit("[ERROR] "+monipath+" does not exist.")
            #monifile = open(monipath,"r")
            if MODE[:4] in ["Disp","Cont","Dist"]: monidata = np.loadtxt(monipath, usecols=1)
            elif MODE[:4] == "Pres": 
                if self.fSteppedRamp or self.vzConstCOM != 0:
                    monidata = np.loadtxt(monipath, usecols=6)
                    #print("moni:", monidata.min(), monidata.max()) # DEBUG
                else: 
                    monidata = self.stiffness0 * np.array([-1, 1])
                    print("maxStiff:", monidata.min(), monidata.max()) # DEBUG
            self.minmaxZ = [monidata.min(), monidata.max()]


class interSheet:
    # Defaults: should match contMech code
    ID = 0
    nx = 0; ny = 0 
    resolMovie = 512
    fPotential = 0
    fDumpFrame = 0
    nTime = 0; nFrames = 0; increment = 1

    konfigName = "inter"#TODO: check if necessary
    contFiles = [] # file names
    attrFiles = []
    data = np.zeros((1,1))
    
    def __init__(self, newID): 
        global INCREMENT
        self.ID = newID
        self.increment = INCREMENT

    def updateFiles(self, path):#, numFrames):
        """ sets up list of frame files belonging to this object """

        global NFRAMES
        
        self.nFrames = NFRAMES

        if not self.fDumpFrame:
            print("Nothing to do in inter"+str(self.ID))
            return
        
        # konfigName
        ending = ".dat"
        if VERSION == "dev": sys.exit("interaction class should not happen here.")#DEBUG
        self.konfigName = "frames/inter" + str(self.ID)
        
        # Individual profile file names
        self.contFiles = []
        iFrame = 0
        for iFrame in range(1,self.nFrames+1):
            #iFrame += 1
            file = path + self.konfigName + "contact." + str(iFrame) + ending
            #print("-> ", file) #DEBUG
            if (not os.path.isfile(file)): break
            self.contFiles.append(file)
        if iFrame < self.nFrames: 
            print("[WARNING]", file, "does not exist.",flush=True)
            self.nTime = int(self.nTime*(iFrame-1)/self.nFrames)
            self.nFrames = iFrame - 1
        
        self.attrFiles = []
        if self.fPotential > 1:
            for iFrame in range(1,self.nFrames+1):
                file = path + self.konfigName + "attract." + str(iFrame) + ending
                #print("-> ", file) #DEBUG
                if (not os.path.isfile(file)): 
                    print("[WARNING]", file, "does not exist.",flush=True)
                else: self.attrFiles.append(file)

        # Sanity check
        if len(self.contFiles) == 0:
            sys.exit("[ERROR] No files found for inter "+str(self.ID))
        if self.fPotential > 1:
            if len(self.contFiles) != len(self.attrFiles): 
                sys.exit("[ERROR] "+str(self.nFrames)+" frames, "+str(len(self.contFiles))+" cont and "+str(len(self.attrFiles))+" attr files.")

    def updateData(self,numFrame):
        """ reads and sets up the data file needed to plot the next frame """
        
        idx = START + self.increment*(numFrame-1-START) # numFrame starts at 1, not 0!

        if not self.fDumpFrame:
            print("Nothing to do in inter"+str(self.ID))
            return

        self.data = np.zeros((self.nx,self.ny),dtype=np.int16)
        try: 
            #Adjust: CONVERSION!!!!
            #f = np.loadtxt(self.contFiles[idx]) 
            #LX = XLIM[1] - XLIM[0]; LY = YLIM[1] - YLIM[0]
            #f[:,0] = f[:,0]*NX/LX; f[:,1] = f[:,1]*NY/LY;
            #f = np.ceil(f[1:,:]).astype(np.uint32)
            f = np.loadtxt(self.contFiles[idx], dtype=np.uint32)
            f = f[1:,:]
            self.data[f[:,0],f[:,1]] = 1 # 1 is repulsive
        except: print("[WARNING] Cont file #"+str(numFrame)+" does not contain valid indices.")
        
        if (self.fPotential > 1) and fATTRACT:
            try: 
                #Adjust: CONVERSION!!!!
                #f = np.loadtxt(self.contFiles[idx]) 
                #LX = XLIM[1] - XLIM[0]; LY = YLIM[1] - YLIM[0]
                #f[:,0] = f[:,0]*NX/LX; f[:,1] = f[:,1]*NY/LY;
                #f = np.ceil(f[1:,:]).astype(np.uint32)
                f = np.loadtxt(self.attrFiles[idx], dtype=np.uint32)
                f = f[1:,:]
                self.data[f[:,0],f[:,1]] = -1 # -1 is attractive
            except: print("[WARNING] Attr file #"+str(numFrame)+" does not contain valid indices.")

        #WIP
        self.data = self.data[RESMOVIE[0],:]
        self.data = self.data[:,RESMOVIE[1]]
        #DEBUG
        #print("contact:", (self.data==1).sum()/(NX*NY))
        #print("attract:", (self.data==-1).sum()/(NX*NY))





#%% ----- Initialization ----- %%#

def checkDir(path):
    """ Checks <path> directory for "frames" subfolder and "params.out". """

    if path[-1] != "/": path += "/"
    if (not os.path.isdir(path+"frames")):
        sys.exit("[ERROR] Could not find the subdirectory 'frames'.")
    
    # Check for params.out file
    if (not os.path.isfile(path + "params.out")):
        sys.exit("[ERROR] Could not find 'params.out'.")


def initParams(path):
    """ initialize params.out
        
        reads <path>/params.out to recreate sheet and interaction setup,
        initializes them and checks corresponding files and data to determine
        optimal plot parameters.
    """

    global XLIM, YLIM, ZLIM, NX, NY, NFRAMES, SHEET, RAMP, LABELS, RESMOVIE
    if path[-1] != "/": path += "/"
    paramsfile = open(path + "params.out","r")
    freqFrame = 1; nTime = 2; dTime = 0; nSheet = 0;
    lengthX = 0; lengthY = 0
    reading = "global"
    for line in paramsfile:
        val = line.split("\t")[0]
        # Read global parameters
        if ("# lengthX" in line): 
            lengthX = float(val)
            XLIM[1] = max(XLIM[1],float(val)/2)
        elif ("# lengthY" in line): 
            lengthY = float(val)
            YLIM[1] = max(YLIM[1],float(val)/2)
        elif ("# nxGlobal" in line): NX = max(NX,int(val))
        elif ("# nyGlobal" in line): NY = max(NY,int(val))
        elif ("# nTime" in line): nTime = int(val)
        elif ("# dTime" in line): dTime = float(val)
        elif ("# freqFrame" in line): freqFrame = int(val)
        elif ("# frameInterval" in line): freqFrame = int(val)
        elif ("# nSheet" in line): nSheet = int(val)
        # Read sheet parameters
        elif ("# sheet start" in line): 
            SHEET.append(gfmdSheet(int(val)))
            reading = "sheet"
        elif ("# nx #" in line): SHEET[-1].nx = int(val)
        elif ("# ny #" in line): SHEET[-1].ny = int(val)
        elif ("# nElast" in line): SHEET[-1].nElast = int(val)
        elif ("# stiffness0" in line): SHEET[-1].stiffness0 = float(val)
        elif ("# fRough" in line): SHEET[-1].fRough = int(val) # works for fRoughAdd and fRoughRead
        elif ("# fTopo" in line): SHEET[-1].fRough = int(val) # works for fRoughAdd and fRoughRead
        elif ("# f3dMovie" in line): SHEET[-1].f3dMovie = int(val)
        elif ("# nVeloTurnStep" in line): SHEET[-1].nVeloTurnStep = int(val)
        elif ("# fSteppedRamp" in line): SHEET[-1].fSteppedRamp = int(val)
        elif ("# vzConstCOM" in line): SHEET[-1].vzConstCOM = float(val)
        elif ("# dzRamp" in line): SHEET[-1].dzRamp = float(val)
        elif ("# rampSteps" in line): SHEET[-1].rampPeriod += int(val)
        elif ("# rampRelax" in line): SHEET[-1].rampPeriod += int(val)
        elif ("# fKelvinVoigt" in line): SHEET[-1].fKelvinVoigt = int(val)
        # Read inter parameters
        elif ("# inter start" in line): 
            INTER.append(interSheet(int(val)))
            reading = "inter"
        elif ("# fDumpFrame" in line): INTER[-1].fDumpFrame += int(val)
        elif ("# fPotential" in line): INTER[-1].fPotential += int(val)
        elif ("# resolMovie" in line) and (reading=="inter"): INTER[-1].resolMovie += int(val)

    
    # Update global params
    if NY == 0: NY = NX
    #WIP
    if RESMOVIE[0] == 0: RESMOVIE[0] = range(NX)
    if RESMOVIE[1] == 0: RESMOVIE[1] = range(NY)
    if lengthY == 0: lengthY = lengthX
    nFrames = int(nTime/freqFrame)
    NFRAMES = nFrames # TODO: this needs to be adjusted if joining multiple paths
    if YLIM[1] == 0: YLIM[1] = XLIM[1]
    XLIM[0] = - XLIM[1]; YLIM[0] = - YLIM[1]    
    paramsfile.close()
    
    # Update SHEET params
    minima = []
    maxima = []
    for iSheet in range(-nSheet,0):
        SHEET[iSheet].dTime = dTime
        SHEET[iSheet].nTime = nTime
        SHEET[iSheet].lengthX = lengthX
        SHEET[iSheet].lengthY = lengthY
        SHEET[iSheet].nFrames = nFrames
        SHEET[iSheet].updateNxNy()
        SHEET[iSheet].updateFiles(path)
        SHEET[iSheet].updateCols()
        SHEET[iSheet].updateDims(path)
        if SHEET[iSheet].rampName != "": RAMP = True
        if (MODE[-2:]=="3D") or SHEET[iSheet].nElast > 0:
            minima.append(SHEET[iSheet].minmaxZ[0])
            maxima.append(SHEET[iSheet].minmaxZ[1])

    # Update inter params
    if MODE=="Cont2D":
        for iInter in range(-len(INTER),0):
            INTER[iInter].updateFiles(path)
            print("resolMovie: ", INTER[iInter].resolMovie)#DEBUG
            if INTER[iInter].resolMovie < NX:
                INTER[iInter].nx = NX // ( (NX-1)//INTER[iInter].resolMovie + 1 )
            else: INTER[iInter].nx = NX
            if INTER[iInter].resolMovie < NY:
                INTER[iInter].ny = NY // ( (NY-1)//INTER[iInter].resolMovie + 1 )
            else: INTER[iInter].ny = NY
        RESMOVIE[0] = range(INTER[IID].nx)
        RESMOVIE[1] = range(INTER[IID].ny)
    # Update ZLIM
    if MODE[-2:]=="3D":
        ZLIM[0] = max(maxima) - min(minima) # TODO: what about joining paths in 3D?
        ZLIM[1] = -ZLIM[0]
    else:
        dummy = ZLIM[0]
        ZLIM[0] = 2.0*min(minima+[dummy]) - 1.0*max(maxima+[ZLIM[1]])
        ZLIM[1] = 2.0*max(maxima+[ZLIM[1]]) - 1.0*min(minima+[dummy])
    #TODO What about 3D with more than 2 sheets?
    if (MODE[-2:]=="3D") and len(SHEET)>2:
        sys.exit("[ERROR] 3D currently only supported for 2 sheets.")

    # Update labels
    for iSheet in range(nSheet):
        if SHEET[iSheet].nElast: LABELS.append("Elastic surface "+str(iSheet))
        else: LABELS.append("Rigid surface "+str(iSheet))



def animate2D():
    """ run animation from 2D files """

    global MODE, XLIM, YLIM, ZLIM, RAMP

    # Setup figure and axes
    if RAMP: 
        fig, axes = plt.subplots(1,2,figsize=[10, 4.2],dpi=DPI)
        ax1 = axes[0]; ax2 = axes[1]
        ax2.set_xlabel("normal displacement ("+UNITS[0]+")")
        if UNITS[1][-1] == "N": ax2.set_ylabel("force ("+UNITS[1]+")")
        else: ax2.set_ylabel("pressure ("+UNITS[1]+")")
        ax2.grid(GRID)
        if XRAMP != None: ax2.set_xlim(XRAMP[0], XRAMP[1])
        if YRAMP != None: ax2.set_ylim(YRAMP[0], YRAMP[1])
    else: fig, ax1 = plt.subplots(figsize=[5.6, 4.2],dpi=DPI)

    # Setup coordinate labels and ranges
    if MODE[:4] == "Cont":
        ax1.set_xlabel("X")
        ax1.set_ylabel("Y")
        ax1.grid(False)
    else:
        if   MODE[-1] == "X": axlim = XLIM
        elif MODE[-1] == "Y": axlim = YLIM
        elif MODE[-1] == "D":
            axlim = [0, np.sqrt(XLIM[1]**2 + YLIM[1]**2)]
            axlim[0] = -axlim[1]
        ax1.set_xlim(1.1*axlim[0], 1.1*axlim[1])
        ax1.set_xlabel("lateral coordinate " + MODE[-1] + " ("+UNITS[0]+")")
        if MODE[:4] == "Disp": 
            ax1.set_ylabel("normal displacement ("+UNITS[0]+")")
            ax1.set_ylim(ZLIM[0], ZLIM[1])
        elif MODE[:4] == "Pres": ax1.set_ylabel("Pressure Z ("+UNITS[1]+")")

    # Lists of plot objects
    lines1 = [0]*len(SHEET)
    lines2 = [0]*len(SHEET)
    lineColor = [0]*len(SHEET);
    
    # ax1: Contact movie
    if MODE[:4] == "Cont":
        aspect = SHEET[1].lengthX / SHEET[1].lengthY# * len(RESMOVIE[1]) / len(RESMOVIE[0])
        im = ax1.imshow(np.zeros((INTER[IID].nx,INTER[IID].ny)), cmap="jet",vmin=-1,vmax=1, aspect=aspect)
        ax1.get_images()[0].set_clim(-1, 1)
        lines1 = im,
        for iSheet in range(len(SHEET)):
            lineColor[iSheet] = plt.rcParams['axes.prop_cycle'].by_key()['color'][iSheet]

    # ax1: Profile animation
    else:
        for iSheet in range(len(SHEET)):
            if (SHEET[iSheet].nElast == 0):
                #pobj, = ax1.plot(SHEET[iSheet].data[:,0], SHEET[iSheet].data[:,1], label=LABELS[iSheet])
                #lineColor[iSheet] = pobj.get_color()
                pobj = ax1.fill_between(SHEET[iSheet].data[:,0], SHEET[iSheet].data[:,1], ZLIM[1], facecolor=(0.8,0.8,0.8), linewidth=0)#FILL
                lineColor[iSheet] = (0.8, 0.8, 0.8)
            else:
                pobj, = ax1.plot([],[], ls=LINESTYLE, marker=MARKERSTYLE, label=LABELS[iSheet])
                lineColor[iSheet] = pobj.get_color()
            lines1[iSheet] = pobj
    
    # ax2: Ramp animation
    if RAMP:
        for iSheet in range(len(SHEET)):
            if (SHEET[iSheet].fSteppedRamp != 0):
                SHEET[iSheet].ramp = np.loadtxt(SHEET[iSheet].rampName, usecols=(0,1))
                if UNITS[0][0]=="µ":
                    SHEET[iSheet].ramp[:,0] = SHEET[iSheet].ramp[:,0]*1e6
                if UNITS[0]=="mm":
                    SHEET[iSheet].ramp[:,0] = SHEET[iSheet].ramp[:,0]*1e3
                if UNITS[1][-1]=="N":
                    SHEET[iSheet].ramp[:,1] = SHEET[iSheet].ramp[:,1]*SHEET[iSheet].lengthX*SHEET[iSheet].lengthY
                if UNITS[1][0]=="m":
                    SHEET[iSheet].ramp[:,1] = SHEET[iSheet].ramp[:,1]*1000
                if UNITS[1][0]=="k":
                    SHEET[iSheet].ramp[:,1] = SHEET[iSheet].ramp[:,1]/1000
                #print("-> loaded ramp") #DEBUG
            elif (SHEET[iSheet].vzConstCOM != 0):
                SHEET[iSheet].ramp = np.loadtxt(SHEET[iSheet].rampName, usecols=(1,6))
                #print("-> loaded moni") #DEBUG
                #if SHEET[iSheet].fKelvinVoigt: SHEET[iSheet].ramp[:,0] = SHEET[iSheet].ramp[:,0] - 0.3e-6#QUICKFIX_KV
                if UNITS[1][-1]=="N":
                    area = SHEET[iSheet].lengthX*SHEET[iSheet].lengthY
                    SHEET[iSheet].ramp[:,1] = SHEET[iSheet].ramp[:,1]*area
                #else: "Pa"
                if UNITS[1][0]=="m":
                    SHEET[iSheet].ramp[:,1] = SHEET[iSheet].ramp[:,1]*1000
                if UNITS[1][0]=="k":
                    SHEET[iSheet].ramp[:,1] = SHEET[iSheet].ramp[:,1]/1000
                #SHEET[iSheet].ramp[:,1] = 0.853*SHEET[iSheet].ramp[:,1]#QUICKFIX_KV_DIRTY
            else:
                lines2[iSheet] = 0 
                #print("-> loaded nothing") #DEBUG
                continue
            ax2.plot(SHEET[iSheet].ramp[:,0], SHEET[iSheet].ramp[:,1], 
                     label=LABELS[iSheet], color=lineColor[iSheet])
            pobj, = ax2.plot(SHEET[iSheet].ramp[START,0], SHEET[iSheet].ramp[START,1], 
                            marker="o", color=lineColor[iSheet])
            lines2[iSheet] = pobj
    if ADDRAMP["plot"]:
        ax2.plot(ADDRAMP["x"],ADDRAMP["y"],ADDRAMP["ls"],label=ADDRAMP["lab"])


    if RAMP:
        ax2.legend()
        plt.tight_layout() # looks better with wide plot

    

    def updateFrame(iFrame):
        """ actual animator function

        goes through all sheet and inter objects and reads in an updates 
        corresponding data for the next frame.
        """

        # Update objects
        if MODE[:4]=="Cont": INTER[IID].updateData(iFrame)
        for iSheet in range(len(SHEET)):
            SHEET[iSheet].updateData(iFrame)

        # Contact movie
        if MODE[:4] == "Cont":
            im = lines1[0]
            im.set_array(INTER[IID].data)

        # Profile animation
        else:
            #Adjust: offset so that rigid body is moving, elastic one is still
            """
            offset = 0; numElast = 0
            for iSheet in range(len(SHEET)):
                if SHEET[iSheet].nElast: 
                    numElast += 1
                    offset += SHEET[iSheet].data[0,1]
            offset /= numElast
            #Adjust end """
            for iSheet,line in enumerate(lines1):
                if not SHEET[iSheet].nElast: continue #Adjust
                line.set_xdata(SHEET[iSheet].data[:,0])
                line.set_ydata(SHEET[iSheet].data[:,1])# - offset) #Adjust
                #if not SHEET[iSheet].nElast: continue #Adjust
                line.set_label("Elastic surface "+str(iSheet)+" Frame "+str(iFrame))

        # Ramp animation
        if RAMP:
            for iSheet,line in enumerate(lines2):
                if (SHEET[iSheet].rampName == ""): continue
                # TODO interpolate in case frames and ramp are not in sync:
                idx = START + SHEET[iSheet].increment*int((iFrame-1-START))
                idx *= len(SHEET[iSheet].ramp)/SHEET[iSheet].nFrames
                idx = int(idx)
                #print("->",str(idx))#DEBUG
                if idx >= len(SHEET[iSheet].ramp): continue
                line.set_xdata(SHEET[iSheet].ramp[idx,0])
                line.set_ydata(SHEET[iSheet].ramp[idx,1])

        # Return plot objects
        if RAMP: return lines1, lines2
        else: # Legend updates each time. if RAMP, no second legend needed
            if MODE[:4] != "Cont": ax1.legend(loc="lower left")
            return lines1
    
    
    def updateTime():
        """ counts up and loops frame index """

        t = 0
        while t<NFRAMES:
            t += animObj.direction
            yield (t%NFRAMES)+1+START
    
    
    def keyPress(event):
        """ adds pause/play and reverse/forward during animation playback """

        if event.key.isspace():
            if animObj.running: animObj.event_source.stop()
            else: animObj.event_source.start()
            animObj.running ^= True
        elif event.key == 'left': animObj.direction = -1
        elif event.key == 'right': animObj.direction = +1
        # Manually update the plot
        if event.key in ['left','right']:
            t = animObj.frame_seq.send(1)
            updateFrame(t)
            plt.draw()
            
    fig.canvas.mpl_connect('key_press_event', keyPress)
    
    # Run the animation
    animObj = FuncAnimation(fig, func=updateFrame, frames=updateTime, interval=1000/FPS, save_count=NFRAMES)
    animObj.running = True
    animObj.direction = +1
    return(animObj)



def animate3D():
    """ run animation from 3D files """

    global MODE

    fig = plt.figure(dpi=DPI)
    plt.xlabel("X")
    plt.ylabel("Y")
    aspect = SHEET[1].lengthX / SHEET[1].lengthY# * len(RESMOVIE[1]) / len(RESMOVIE[0])
    im = plt.imshow(np.zeros((len(RESMOVIE[0]),len(RESMOVIE[1]))), cmap="jet", aspect=aspect)
    if MODE=="Cont3D": plt.clim(0,1)
    elif MODE=="Press3D": plt.clim(-SHEET[1].stiffness0/4, SHEET[1].stiffness0/40)
    else: plt.clim(ZLIM[0], ZLIM[1])
    #ax.set_xticks(np.arange(5)*NX/4)
    #ax.set_xticklabels([str((XLIM[1]-XLIM[0])*val/4) for val in range(5)])
    #ax.set_yticks(np.arange(5)*NY/4)
    #ax.set_yticklabels([str((YLIM[1]-YLIM[0])*val/4) for val in range(5)])
    if MODE!="Cont3D": plt.colorbar()
    
    def updateFrame(iFrame):
        """ actual animator function

        goes through all sheet objects and reads in an updates 
        corresponding data for the next frame.
        """

        for iSheet in range(len(SHEET)): SHEET[iSheet].updateData(iFrame)
        #TODO might want to be transposed or rot90 is correct
        if MODE=="Cont3D": array = SHEET[0].data <= SHEET[1].data
        elif MODE=="Press3D": array = SHEET[1].data # != 0 # < -SHEET[1].stiffness0/100
        elif MODE=="Dist3D": array = SHEET[0].data - SHEET[1].data
        elif MODE=="Disp3D": array = SHEET[1].data
        im.set_array(array)
        return im,
    
    
    def updateTime():
        """ counts up and loops frame index """

        t = START
        while t<NFRAMES:
            t += animObj.direction
            yield t%NFRAMES+1+START
    

    def keyPress(event):
        """ adds pause/play and reverse/forward during animation playback """

        if event.key.isspace():
            if animObj.running:
                animObj.event_source.stop()
            else:
                animObj.event_source.start()
            animObj.running ^= True
        elif event.key == 'left':
            animObj.direction = -1
        elif event.key == 'right':
            animObj.direction = +1
        # Manually update the plot
        if event.key in ['left','right']:
            t = animObj.frame_seq.send(1)
            updateFrame(t)
            plt.draw()
            
    fig.canvas.mpl_connect('key_press_event', keyPress)
    
    animObj = FuncAnimation(fig, func=updateFrame, frames=updateTime, interval=1000/FPS, save_count=NFRAMES)
    animObj.running = True
    animObj.direction = +1
    return(animObj)
    

def dataReset(): 
    """ hard reset """

    global SHEET, NX, NY, XLIM,YLIM, ZLIM, NFRAMES
    del SHEET[:], INTER[:]
    NX = 0; NY = 0;
    XLIM = [0,0]; YLIM = [0,0]
    ZLIM = [+1.0e+38,-1.0e+38]
    NFRAMES = 1



def save(animation, filename):
    """ save animation to file 

    For mp4 files, make sure ffmpg and h264 are installed correctly:
    conda install -c conda-forge ffmpeg
    conda install -c conda-forge x264=20131218

    (see: https://stackoverflow.com/questions/13316397/matplotlib-animation-no-moviewriters-available)
    (see: https://github.com/sigsep/sigsep-mus-db/issues/14)
    """

    if filename[-3:]=="mp4": 
        from matplotlib.animation import FFMpegWriter
        ffmpeg = FFMpegWriter(fps=FPS)
        animation.save(filename, writer=ffmpeg)
    else: animation.save(filename, writer='imagemagick', fps=FPS)



def initGlobals(paths):
    """ initializes global parameters
    
    goes through all simulation directories in <paths> and tries to find optimal 
    parameters to plot all simulations into one animation
    """

    global SHOW, FPS, LINESTYLE, MARKERSTYLE, MODE, VERSION, UNITS, NX, NY, XLIM, YLIM, ZLIM, XRAMP, YRAMP, INCREMENT, START, END, NFRAMES, RAMP, IID

    # Init prams file: TODO: before or after the rest of this function?
    if (os.path.isfile("input.py")):
        import input
        print("Using input.py")
        if input.SHOW:      SHOW = input.SHOW
        if input.FPS:       FPS = input.FPS
        if input.LINESTYLE: LINESTYLE = input.LINESTYLE
        if input.MARKERSTYLE: MARKERSTYLE = input.MARKERSTYLE
        if input.MODE:      MODE = input.MODE
        if input.VERSION:   VERSION = input.VERSION
        if input.NX:        NX = input.NX
        if input.NY:        NY = input.NY
        if input.XLIM:      XLIM = input.XLIM
        if input.YLIM:      YLIM = input.YLIM
        if input.XRAMP:     XRAMP = input.XRAMP
        if input.YRAMP:     YRAMP = input.YRAMP
        if input.INCREMENT: INCREMENT = input.INCREMENT
        if input.ZLIM:      ZLIM = input.ZLIM
        if input.START:     START = input.START
        if input.END:       END = input.END
        if input.NFRAMES:   NFRAMES = input.NFRAMES
        if input.RAMP:      RAMP = input.RAMP
        if input.IID:       IID = input.IID
        if input.UNITS:     UNITS = input.UNITS
        return

    # Only relevant for joining multiple paths:
    if len(paths) > 1:
        travelPerFrame = []
        for iSheet in range(len(SHEET)):
            if SHEET[iSheet].nElast == 0: travelPerFrame.append(0)
            elif SHEET[iSheet].vzConstCOM != 0:
                travelPerFrame.append( abs(SHEET[iSheet].vzConstCOM*SHEET[iSheet].dTime*
                                           SHEET[iSheet].nTime/SHEET[iSheet].nFrames) )
            else: travelPerFrame.append( abs(1.*SHEET[iSheet].nTime * SHEET[iSheet].dzRamp / 
                                            (SHEET[iSheet].rampPeriod * SHEET[iSheet].nFrames)) )
        travelRef = max(travelPerFrame)
        nFrames = []
        for iSheet in range(len(SHEET)):
            if travelPerFrame[iSheet] > 0:
                SHEET[iSheet].increment = INCREMENT*int(travelRef/travelPerFrame[iSheet])
                nFrames.append( SHEET[iSheet].nFrames/SHEET[iSheet].increment )
        
        # Update NFRAMES
        if END > 0: NFRAMES = int(min(nFrames+[END]) - START)
        else: NFRAMES = int(min(nFrames) - START)




def dumpGlobals(path):
    """ shows automatically adjusted parameters for user to check """

    print("\nPath:", str(path), "\n")

    print("SHOW =", str(SHOW))
    print("FPS =", FPS)
    print("DPI =", DPI)
    print("GRID = %s\n" % str(GRID))
    print("LINESTYLE = \"%s\"" % LINESTYLE)
    print("MARKERSTYLE = \"%s\"" % MARKERSTYLE)
    print("LABELS = \"%s\"" % LABELS)
    print("MODE = \"%s\"" % MODE)
    print("VERSION = \"%s\"" % VERSION)
    print("UNITS = %s\n" % str(UNITS))

    print("NX = %i" % NX)
    print("NY = %i" % NY)
    print("XLIM = %s" % str(XLIM))
    print("YLIM = %s" % str(YLIM))
    print("ZLIM = %s" % str(ZLIM))
    print("XRAMP = %s" % str(XRAMP))
    print("YRAMP = %s" % str(YRAMP))
    print("INCREMENT = %i" % INCREMENT)
    print("START = %i" % START)
    print("END = %i" % END)
    print("NFRAMES = %i" % NFRAMES)
    print("RAMP = %s" % str(RAMP))
    print("IID = %i\n" % IID)
    print("RESMOVIE = %s\n" % str(RESMOVIE))
    for iSheet in range(len(SHEET)):
        SHEET[iSheet].increment = INCREMENT*SHEET[iSheet].increment
        #if SHEET[iSheet].nElast:#DEBUG
        #    print("SHEET["+str(iSheet)+"] increment: ", SHEET[iSheet].increment)


def init(mode,paths):
    global MODE; MODE = mode
    dataReset()
    if type(paths)==str: paths = [paths]
    for path in paths:
        checkDir(path)
        initParams(path)
    initGlobals(paths)
    dumpGlobals(paths)

def run():
    global MODE
    if MODE[-2:] == "3D": animation = animate3D()
    else: animation = animate2D()
    if SHOW: plt.show()
    return(animation)





#%% ----- Functions for use in iPython console ----- %%#


def runDispX(path=".",*morepaths):
    """ animates displacement along X direction

    with stress-displacement curve if applicable.
    """

    init("DispX",[path]+list(morepaths))
    return( run() )

def runDispY(path=".",*morepaths):
    """ animates displacement along Y direction

    with stress-displacement curve if applicable.
    """
    
    init("DispY",[path]+list(morepaths))
    return( run() )

def runDispD(path=".",*morepaths):
    """ animates displacement along diagonal direction

    with stress-displacement curve if applicable.
    """
    
    init("DispD",[path]+list(morepaths))
    return( run() )

def runCont2D(path=".",*morepaths):
    """ animates contact movie from interSheet frames

    with stress-displacement curve if applicable.
    """
    
    init("Cont2D",[path]+list(morepaths))
    return( run() )

def runCont3D(path=".",*morepaths):
    """ animates contact movie from sheet configs """

    init("Cont3D",[path]+list(morepaths))
    return( run() )

def runDisp3D(path=".",*morepaths):
    """ animates displacement movie from elastic sheet config """

    init("Disp3D",[path]+list(morepaths))
    return( run() )

def runDist3D(path=".",*morepaths):
    """ animates surface distance movie from sheet configs """

    init("Dist3D",[path]+list(morepaths))
    return( run() )



def runPressX(path=".",*morepaths):
    """ animates pressure along X direction

    with stress-displacement curve if applicable.
    """

    print("[WARNING] Pressure animation supported as of contMech version March 2021.")
    init("PressX",[path]+list(morepaths))
    return( run() )
    

def runPressY(path=".",*morepaths):
    """ animates displacement along Y direction

    with stress-displacement curve if applicable.
    """

    print("[WARNING] Pressure animation supported as of contMech version March 2021.")
    init("PressY",[path]+list(morepaths))
    return( run() )
    

def runPressD(path=".",*morepaths):
    """ animates displacement along diagonal direction

    with stress-displacement curve if applicable.
    """

    print("[WARNING] Pressure animation supported as of contMech version March 2021.")
    init("PressD",[path]+list(morepaths))
    return( run() )


def reanimate():
    """ runs the previous animation again.

    helpful after manually changing parameters
    """

    initGlobals(["."])
    animation = animate2D()
    if SHOW: plt.show()
    return(animation)

