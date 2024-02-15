#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
-------------------------------------
  Interface for contMech Parameters
-------------------------------------

@author: thescientist
"""

import os, sys

paramType = {
  # --- SHEET : GENERAL --- #
  'ID' : 'int',
  'nElast' : 'int',
  'stiffness0' : 'float',
  'contactMod' : 'float',#OLD
  'fThickness0' : 'int',
  'poisson0' : 'float',
  'thickness0' : 'float',
  'elastExpnt0' : 'float',
  'stiffness1' : 'float',
  'fThickness1' : 'int',
  'poisson1' : 'float',
  'thickness1' : 'float',
  'elastExpnt1' : 'float',
  'fDispX' : 'int',
  'fDispY' : 'int',
  # mechanical: viscoelasticity
  'fKelvinVoigt' : 'int',
  'tauKV' : 'float',
  'scalKV' : 'float',
  'rKV_LPF' : 'int',
  'fMaxwell' : 'int',#NEW
  'nMaxwell' : 'int',#NEW
  'massGFMD' : 'float',#NEW
  'fDumpMaxwell' : 'int',#NEW
  # mechanical: GFMD dynamics
  'fMassWeightg' : 'int',
  'zeroModeMass' : 'float',

  # --- SHEET : TOPOGRAPHY --- #
  'fRoughRead' : 'int',
  'fRoughAdd' : 'int',
  'fTopoRead' : 'int',#NEW
  'fTopoAdd' : 'int',#NEW
  'dzStep' : 'float',
  # topography: Hertz
  'rXhertz' : 'float',
  'rYhertz' : 'float',
  'hertzExpnt' : 'float',
  # topography: self-affine fractal
  'hurst' : 'float',
  'fRollOff' : 'int',
  'lambdaR' : 'float',
  'lambdaS' : 'float',
  'peklenik' : 'float',
  'fRoughNorm' : 'int',
  'rRoughNorm' : 'float',
  'fBoxMuller' : 'int',
  # topography: single wave roughness
  'fAddSWR' : 'int',
  'nqxAddSWR' : 'int',
  'nqyAddSWR' : 'int',
  'heightSWR' : 'float',
  # topography: flat punch
  'radiusFlatPunch' : 'float',
  'heightFlatPunch' : 'float',
  # topography: sphere
  'rSphere' : 'float',

  # --- SHEET : PRESSURE/DISPLACEMENT CONTROL --- #
  'pressInit' : 'float',
  'pressFinal' : 'float',
  'forceInit' : 'float',#NEW
  'forceFinal' : 'float',#NEW
  'pressPhi' : 'float',#NEW
  'pressTheta' : 'float',#NEW
  'fConstCOM' : 'int',
  'fzOpposite' : 'int',
  'fExtrapolateInf' : 'int',#NEW
  'zConstCOM' : 'float',
  'zOpposite' : 'float',
  'nVeloTurnStep' : 'int',
  'rVeloTurnStep' : 'float',
  'veloTurnPress' : 'float',
  'nVeloTransition' : 'int',#NEW
  'pressTurnaround' : 'float',#NEW
  'forceTurnaround' : 'float',#NEW
  'vzConstCOM' : 'float',
  # control: ramp
  'fSteppedRamp' : 'int',
  'rampSteps' : 'int',
  'rampRelax' : 'int',
  'dzRamp' : 'float',
  'dpRamp' : 'float',
  # control: lateral
  'fLateral' : 'int',
  'vX' : 'float',
  'vY' : 'float',
  'frictRelax' : 'int',
  'fConstCOMx' : 'int',#NEW
  'xConstCOM' : 'float',#NEW
  'vxConstCOM' : 'float',#NEW
  'fConstCOMy' : 'int',#NEW
  'yConstCOM' : 'float',#NEW
  'vyConstCOM' : 'float',#NEW
  # control: steady state sliding
  'fSteadySlide' : 'int',

  # --- SHEET : ON-SITE --- #
  'fOnSitePotential' : 'int',
  'fOnSitePeriod' : 'float',
  'gammaOnSite' : 'float',
  'frictionCoeffOS' : 'float',
  'vXOnSite' : 'float',
  'vYOnSite' : 'float',
  'xOnSite' : 'float',#NEW
  'yOnSite' : 'float',#NEW
  'surfEnergOS' : 'float',
  'frictionTauOS' : 'float',#NEW
  'potCurveRelOS' : 'float',

  # --- SHEET : OBSERVABLES --- #
  'resolMovie' : 'int',
  'f3dMovie' : 'int',
  
  # --- INTERACTION : GENERAL --- #
  'sheetID0' : 'int',
  'sheetID1' : 'int',
  'fConstraint' : 'int',
  'fPotential' : 'int',
  'surfEnerg' : 'float',
  'potRange' : 'float',
  'potCurve' : 'float',
  'potCurveRel' : 'float',
  'relRepulsion' : 'float',
  'fPotentialTest' : 'int',
  'nTimeOn' : 'int', #NEW
  'nTimeOff' : 'int',
  'frictionCoeff' : 'float',
  'frictionTau' : 'float',#NEW
  'relFrictPeak' : 'float',
  'velFrictThresh' : 'float',

  # --- INTERACTION : OBSERVABLES --- #
  'fDumpGap' : 'int',
  'fDumpFrame' : 'int',
  'fDumpLateral' : 'int',
  'resolMovie' : 'int',

  # --- GLOBAL : GENERAL --- #
  'lengthX' : 'float',
  'lengthY' : 'float',
  'nxGlobal' : 'int',
  'nyGlobal' : 'int',
  'nSheet' : 'int',
  'nInter' : 'int',

  'nTime' : 'int',
  'nRelax' : 'int',
  'dTime' : 'float',

  'randSeed' : 'int',
  'dampGlobal' : 'float',

  # --- GLOBAL : FIRE --- #
  'fFire' : 'int',
  'fireDecrmt' : 'float',
  'fireIncrmt' : 'float',
  'fireRedrct' : 'float',

  # --- GLOBAL : TEMPERATURE --- #
  'fLangevin' : 'int',
  'tempInit' : 'float',
  'tempFinal' : 'float',

  # --- GLOBAL : OBSERVABLES --- #
  'freqFrame' : 'int',
  'frameInterval' : 'int',#NEW
  'fLogMeasure' : 'int'#NEW
}

class gfmdSheet:
  konfigName = "konfig"
  ID = 0
  def __init__(self,newID):
    self.ID = newID
    self.konfigName += str(newID)

  # --- MECHANICAL PROPERTIES --- #
  nElast = 0
  contactMod = 1 #OLD
  stiffness0 = 0.5
  fThickness0 = 0
  poisson0 = 0.25
  thickness0 = 1
  elastExpnt0 = 1
  #stiffness1 = 0.
  #fThickness1 = 0
  #poisson1 = 0.25
  #thickness1 = 1
  #elastExpnt1 = 1
  fDispX = fDispY = 0
  # mechanical: viscoelasticity
  fKelvinVoigt = 0
  tauKV = 1
  scalKV = 1000
  rKV_LPF = 4
  fMaxwell = 0
  nMaxwell = 1
  massGFMD = 1.
  fDumpMaxwell = 0
  # mechanical: GFMD dynamics
  fMassWeightg = 1
  zeroModeMass = 1


  # --- TOPOGRAPHY --- #
  fRoughRead = 0; fRoughAdd = 0
  fTopoRead = 0; fTopoAdd = 0 #NEW
  dzStep = 0
  # topography: Hertz
  rXhertz = 0
  rYhertz = -1
  hertzExpnt = 2
  # topography: self-affine fractal
  hurst = 0.8
  fRollOff = 1
  lambdaR = 0.5
  lambdaS = 0.05
  peklenik = 1
  fRoughNorm = 1
  rRoughNorm = 1
  fBoxMuller = 0
  # topography: single wave roughness
  fAddSWR = 0
  nqxAddSWR = nqyAddSWR = 2
  heightSWR = 0.1
  # topography: flat punch
  radiusFlatPunch = 0.25
  heightFlatPunch = 1
  # topography: sphere
  rSphere = 0.1


  # --- PRESSURE/DISPLACEMENT CONTROL --- #
  pressInit = pressFinal = 0
  fConstCOM = fzOpposite = fExtrapolateInf = 0
  zConstCOM = zOpposite = 0
  nVeloTurnStep = 0; rVeloTurnStep = 0
  veloTurnPress = 0
  vzConstCOM = 0
  # control: ramp
  fSteppedRamp = 0
  rampSteps = 150; rampRelax = 350
  dzRamp = dpRamp = 0
  # control: lateral
  fLateral = 0
  vX = vY = 0
  frictRelax = 500 #NEW


  # --- ON-SITE --- #
  fOnSitePotential = 0
  fOnSitePeriod = 0
  gammaOnSite = 0
  frictionCoeffOS = 0
  vXOnSite = vYOnSite = 0


  # --- OBSERVABLES --- #
  resolMovie = 512
  f3dMovie = 0


class interSheet:
  ID = 0
  sheetID0 = 0; sheetID1 = 1
  def __init__(self,newID):
    self.ID = newID

  # --- INTERACTION --- #
  fConstraint = 0
  fPotential = 0
  surfEnerg = 0
  potRange = 0 
  potCurve = 0
  potCurveRel = 0
  fPotentialTest = 0
  nTimeOn = 0 #NEW
  nTimeOff = 0
  frictionCoeff = 0
  frictRelax = 1 #NEW

  # --- OBSERVABLES --- #
  fDumpGap = 0
  fDumpFrame = 0
  fDumpLateral = 0
  resolMovie = 512


class simulation:
  path = ''
  sheet = []
  inter = []
  SHEET = []
  INTER = []
  def __init__(self, path):
    self.path = os.path.split(path)[0]
    self.sheet.clear()
    self.inter.clear()
    self.SHEET = self.sheet # old name for compatibility
    self.INTER = self.inter # old name for compatibility

  lengthX = lengthY = 0
  nxGlobal = nyGlobal = 0
  nSheet = 0
  nInter = 0

  nTime = 0
  dTime = 0
  nRelax = 0

  randSeed = 4712
  dampGlobal = 1.5

  # --- FIRE --- #
  fFire = 0
  fireDecrmt = 1
  fireIncrmt = 1
  fireRedrct = 0

  # --- TEMPERATURE --- #
  fLangevin = 0
  tempInit = tempFinal = 0

  # --- OBSERVABLES --- #
  freqFrame = 0
  frameInterval = 0
  fLogMeasure = 0


def read(file):
  if not os.path.isfile(file): sys.exit('not a valid params file: '+file)

  sim = simulation(file)
  sheet = gfmdSheet(-1)
  inter = interSheet(-1)
  sim_attributes = sim.__dir__()
  sim_attributes = [a for a in sim_attributes if a[:2]!='__']
  sheet_attributes = sheet.__dir__()
  sheet_attributes = [a for a in sheet_attributes if a[:2]!='__']
  inter_attributes = inter.__dir__()
  inter_attributes = [a for a in inter_attributes if a[:2]!='__']

  with open(file, 'r') as fid:
    sheetID = -1
    interID = -1
    
    for iLine,line in enumerate(fid.readlines()):
      # special lines
      if '# nSheet #' in line: 
        sim.nSheet = int(line.split()[0])
        for i in range(sim.nSheet): sim.sheet.append(gfmdSheet(i))
        continue
      if '# nInter #' in line: 
        sim.nInter = int(line.split()[0])
        for i in range(sim.nInter): sim.inter.append(interSheet(i))
        continue
      if '# sheet start' in line: 
        sheetID = int(line.split()[0])
        continue
      if '# sheet end' in line: 
        sheetID = -1
        continue
      if '# inter start' in line: 
        interID = int(line.split()[0])
        continue
      if '# inter end' in line: 
        interID = -1
        continue

      for attribute in sim_attributes:
        if f'# {attribute} #' in line:
          try:
            if paramType[attribute] == 'int': sim.__setattr__(attribute, int(line.split()[0]))
            else: sim.__setattr__(attribute, float(line.split()[0]))
          except:
            if (interID >= 0) or (sheetID >= 0):
              msg = '[WARNING] cannot process line %4i (' % iLine
              msg += ' sheetID='+str(sheetID)+'/'+str(len(sim.sheet))
              msg += ' interID='+str(interID)+'/'+str(len(sim.inter))
              print(msg+' ):')
              print(' '+line[:-1])

      if sheetID >=0:
        for attribute in sheet_attributes:
          if f'# {attribute} #' in line:
            if paramType[attribute] == 'int': sim.sheet[sheetID].__setattr__(attribute, int(line.split()[0]))
            else: sim.sheet[sheetID].__setattr__(attribute, float(line.split()[0]))

      if interID >=0:
        for attribute in inter_attributes:
          if f'# {attribute} #' in line:
            if paramType[attribute] == 'int': sim.inter[interID].__setattr__(attribute, int(line.split()[0]))
            else: sim.inter[interID].__setattr__(attribute, float(line.split()[0]))

  # update default or derived parameters
  if (sim.freqFrame > 0): sim.frameInterval = sim.freqFrame
  if (sim.frameInterval > 0): sim.freqFrame = sim.frameInterval
  for iSheet in range(sim.nSheet):
    if (sim.sheet[iSheet].fRoughAdd > 0): sim.sheet[iSheet].fTopoAdd = sim.sheet[iSheet].fRoughAdd
    if (sim.sheet[iSheet].fRoughRead > 0): sim.sheet[iSheet].fTopoRead = sim.sheet[iSheet].fRoughRead
    if (sim.sheet[iSheet].fTopoAdd > 0): sim.sheet[iSheet].fRoughAdd = sim.sheet[iSheet].fTopoAdd
    if (sim.sheet[iSheet].fTopoRead > 0): sim.sheet[iSheet].fRoughRead = sim.sheet[iSheet].fTopoRead
    if sim.sheet[iSheet].rYhertz == -1: 
      sim.sheet[iSheet].rYhertz = sim.sheet[iSheet].rXhertz
    if sim.sheet[iSheet].fDispY: sim.sheet[iSheet].nElast = 2
    if sim.sheet[iSheet].fDispX: sim.sheet[iSheet].nElast = 3
    if sim.sheet[iSheet].fRoughAdd or sim.sheet[iSheet].fRoughRead:
      sim.sheet[iSheet].konfigName += "E"
    disps = {0 : "", 1 : "D", 2 : "Dzy", 3 : "Dzyx"}
    sim.sheet[iSheet].konfigName += disps[sim.sheet[iSheet].nElast] + ".dat"
  if sim.lengthY == 0: sim.lengthY = sim.lengthX
  if sim.nyGlobal == 0: sim.nyGlobal = sim.nxGlobal

  # old name for compatibility
  sim.SHEET = sim.sheet
  sim.INTER = sim.inter
  
  return sim