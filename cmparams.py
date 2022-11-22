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
  'fDispX' : 'int',
  'fDispY' : 'int',
  # mechanical: viscoelasticity
  'fKelvinVoigt' : 'int',
  'tauKV' : 'float',
  'scalKV' : 'float',
  'rKV_LPF' : 'int',
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
  'fConstCOM' : 'int',
  'fzOpposite' : 'int',
  'fExtrapolateInf' : 'int',#NEW
  'zConstCOM' : 'float',
  'zOpposite' : 'float',
  'nVeloTurnStep' : 'int',
  'rVeloTurnStep' : 'float',
  'veloTurnPress' : 'float',
  'pressTurnaround' : 'float',#NEW
  'vzConstCOM' : 'float',
  # control: ramp
  'fSteppedRamp' : 'int',
  'rampSteps' : 'int',
  'rampRelax' : 'int',
  'dzRamp' : 'float',
  # control: lateral
  'fLateral' : 'int',
  'vX' : 'float',
  'vY' : 'float',
  'frictRelax' : 'int',

  # --- SHEET : ON-SITE --- #
  'fOnSitePotential' : 'int',
  'fOnSitePeriod' : 'float',
  'gammaOnSite' : 'float',
  'frictionCoeffOS' : 'float',
  'vXOnSite' : 'float',
  'vYOnSite' : 'float',

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
  'potCurveRel' : 'float',
  'fPotentialTest' : 'int',
  'nTimeOff' : 'int',
  'frictionCoeff' : 'float',

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
  ID = 0
  def __init__(self,newID):
    self.ID = newID

  # --- MECHANICAL PROPERTIES --- #
  nElast = 0
  stiffness0 = 0.5
  contactMod = 1 #OLD
  fThickness0 = 0
  poisson0 = 0.25
  thickness0 = 1
  elastExpnt0 = 1
  fDispX = fDispY = 0
  # mechanical: viscoelasticity
  fKelvinVoigt = 0
  tauKV = 1
  scalKV = 1000
  rKV_LPF = 4
  # mechanical: GFMD dynamics
  fMassWeightg = 1
  zeroModeMass = 1


  # --- TOPOGRAPHY --- #
  fRoughRead = 0; fRoughAdd = 0
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
  potCurveRel = 0
  fPotentialTest = 0
  nTimeOff = 0
  frictionCoeff = 0

  # --- OBSERVABLES --- #
  fDumpGap = 0
  fDumpFrame = 0
  fDumpLateral = 0
  resolMovie = 512


class simulation:
  path = ''
  SHEET = []
  INTER = []
  def __init__(self, file):
    self.path = file
    self.SHEET.clear()
    self.INTER.clear()

  lengthX = lengthY = 0
  nxGlobal = nyGlobal = 0
  nSheet = 0
  nInter = 0

  nTime = 0
  dTime = 0

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
  fLogMeasure = 0


def read(file):
  if not os.path.isfile(file): sys.exit('not a valid params file: '+file)

  sim = simulation(file)

  with open(file, 'r') as fid:
    sheetID = 666
    interID = 666
    
    for iLine,line in enumerate(fid.readlines()):

      try:
        # --- GLOBAL --- #
        
        # basic
        if '# nxGlobal #' in line: sim.nxGlobal = int(line.split()[0])
        elif '# nyGlobal #' in line: sim.nyGlobal = int(line.split()[0])
        elif '# lengthX #' in line: sim.lengthX = float(line.split()[0])
        elif '# lengthY #' in line: sim.lengthY = float(line.split()[0])
        elif '# nSheet #' in line: 
          sim.nSheet = int(line.split()[0])
          for i in range(sim.nSheet): sim.SHEET.append(gfmdSheet(i))
        elif '# nInter #' in line: 
          sim.nInter = int(line.split()[0])
          for i in range(sim.nInter): sim.INTER.append(interSheet(i))

        elif '# nTime #' in line: sim.nTime = int(line.split()[0])
        elif '# dTime #' in line: sim.dTime = float(line.split()[0])

        elif '# randSeed #' in line: sim.randSeed = int(line.split()[0])
        elif '# dampGlobal #' in line: sim.dampGlobal = float(line.split()[0])
        # FIRE
        elif '# fFire #' in line: sim.fFire = int(line.split()[0])
        elif '# fireDecrmt #' in line: sim.fireDecrmt = float(line.split()[0])
        elif '# fireIncrmt #' in line: sim.fireIncrmt = float(line.split()[0])
        elif '# fireRedrct #' in line: sim.fireRedrct = float(line.split()[0])
        # temperature
        elif '# fLangevin #' in line: sim.fLangevin = int(line.split()[0])
        elif '# tempInit #' in line: sim.tempInit = float(line.split()[0])
        elif '# tempFinal #' in line: sim.tempFinal = float(line.split()[0])
        # observables
        elif '# freqFrame #' in line: sim.freqFrame = int(line.split()[0])
        elif '# frameInterval #' in line: sim.freqFrame = int(line.split()[0])
        elif '# fLogMeasure #' in line: sim.fLogMeasure = int(line.split()[0])

        
        
        # --- gfmdSheet --- #
        
        elif '# sheet start' in line: sheetID = int(line.split()[0])
        elif '# sheet end' in line: sheetID = 666
        
        # mechanical: basic
        elif '# nElast #' in line: sim.SHEET[sheetID].nElast = int(line.split()[0])
        elif '# stiffness0 #' in line: sim.SHEET[sheetID].stiffness0 = float(line.split()[0])
        elif '# contactMod #' in line: 
          sim.SHEET[sheetID].contactMod = float(line.split()[0])
          print('[WARNING] contactMod is deprecated.')
        elif '# fThickness0 #' in line: sim.SHEET[sheetID].fThickness0 = int(line.split()[0])
        elif '# poisson0 #' in line: sim.SHEET[sheetID].poisson0 = float(line.split()[0])
        elif '# thickness0 #' in line: sim.SHEET[sheetID].thickness0 = float(line.split()[0])
        elif '# elastExpnt0 #' in line: sim.SHEET[sheetID].elastExpnt0 = float(line.split()[0])
        elif '# fDispX #' in line: sim.SHEET[sheetID].fDispX = int(line.split()[0])
        elif '# fDispY #' in line: sim.SHEET[sheetID].fDispY = int(line.split()[0])
        # mechanical: viscoelastic
        elif '# fKelvinVoigt #' in line: sim.SHEET[sheetID].fKelvinVoigt = int(line.split()[0])
        elif '# tauKV #' in line: sim.SHEET[sheetID].tauKV = float(line.split()[0])
        elif '# scalKV #' in line: sim.SHEET[sheetID].scalKV = float(line.split()[0])
        elif '# rKV_LPF #' in line: sim.SHEET[sheetID].rKV_LPF = int(line.split()[0])
        # mechanical: GFMD dynamics
        elif '# fMassWeightg #' in line: sim.SHEET[sheetID].fMassWeightg = int(line.split()[0])
        elif '# zeroModeMass #' in line: sim.SHEET[sheetID].zeroModeMass = float(line.split()[0])
        
        # topography: basic
        elif '# fRoughRead #' in line: sim.SHEET[sheetID].fRoughRead = int(line.split()[0])
        elif '# fTopoRead #' in line: sim.SHEET[sheetID].fRoughRead = int(line.split()[0])
        elif '# fRoughAdd #' in line: sim.SHEET[sheetID].fRoughAdd = int(line.split()[0])
        elif '# fTopoAdd #' in line: sim.SHEET[sheetID].fRoughAdd = int(line.split()[0])
        elif '# dzStep #' in line: sim.SHEET[sheetID].dzStep = float(line.split()[0])      
        # topography: Hertz
        elif '# rXhertz #' in line: sim.SHEET[sheetID].rXhertz = float(line.split()[0])
        elif '# rYhertz #' in line: sim.SHEET[sheetID].rYhertz = float(line.split()[0])
        elif '# hertzExpnt #' in line: sim.SHEET[sheetID].hertzExpnt = float(line.split()[0])
        # topography: self-affine fractal
        elif '# hurst #' in line: sim.SHEET[sheetID].hurst = float(line.split()[0])
        elif '# fRollOff #' in line: sim.SHEET[sheetID].fRollOff = int(line.split()[0])
        elif '# lambdaR #' in line: sim.SHEET[sheetID].lambdaR = float(line.split()[0])
        elif '# lambdaS #' in line: sim.SHEET[sheetID].lambdaS = float(line.split()[0])
        elif '# peklenik #' in line: sim.SHEET[sheetID].peklenik = float(line.split()[0])
        elif '# fRoughNorm #' in line: sim.SHEET[sheetID].fRoughNorm = int(line.split()[0])
        elif '# rRoughNorm #' in line: sim.SHEET[sheetID].rRoughNorm = float(line.split()[0])
        elif '# fBoxMuller #' in line: sim.SHEET[sheetID].fBoxMuller = int(line.split()[0])
        # topography: single wave roughness
        elif '# fAddSWR #' in line: sim.SHEET[sheetID].fAddSWR = int(line.split()[0])
        elif '# nqxAddSWR #' in line: sim.SHEET[sheetID].nqxAddSWR = int(line.split()[0])
        elif '# nqyAddSWR #' in line: sim.SHEET[sheetID].nqyAddSWR = int(line.split()[0])
        elif '# heightSWR #' in line: sim.SHEET[sheetID].heightSWR = float(line.split()[0])
        # topography: flat punch
        elif '# radiusFlatPunch #' in line: sim.SHEET[sheetID].radiusFlatPunch = float(line.split()[0])
        elif '# heightFlatPunch #' in line: sim.SHEET[sheetID].heightFlatPunch = float(line.split()[0])
        # topography: sphere
        elif '# rSphere #' in line: sim.SHEET[sheetID].rSphere = float(line.split()[0])

        # control: basic
        elif '# pressInit #' in line: sim.SHEET[sheetID].pressInit = float(line.split()[0])
        elif '# pressFinal #' in line: sim.SHEET[sheetID].pressFinal = float(line.split()[0])
        elif '# fConstCOM #' in line: sim.SHEET[sheetID].fConstCOM = int(line.split()[0])
        elif '# fzOpposite #' in line: sim.SHEET[sheetID].fzOpposite = int(line.split()[0])
        elif '# fExtrapolateInf #' in line: sim.SHEET[sheetID].fExtrapolateInf = int(line.split()[0])#NEW
        elif '# zConstCOM #' in line: sim.SHEET[sheetID].zConstCOM = float(line.split()[0])
        elif '# zOpposite #' in line: sim.SHEET[sheetID].zOpposite = float(line.split()[0])
        elif '# nVeloTurnStep #' in line: sim.SHEET[sheetID].nVeloTurnStep = int(line.split()[0])
        elif '# rVeloTurnStep #' in line: sim.SHEET[sheetID].rVeloTurnStep = float(line.split()[0])
        elif '# veloTurnPress #' in line: sim.SHEET[sheetID].veloTurnPress = float(line.split()[0])#NEW
        elif '# pressTurnaround #' in line: sim.SHEET[sheetID].veloTurnPress = float(line.split()[0])#NEW
        elif '# vzConstCOM #' in line: sim.SHEET[sheetID].vzConstCOM = float(line.split()[0])
        # control: ramp
        elif '# fSteppedRamp #' in line: sim.SHEET[sheetID].fSteppedRamp = int(line.split()[0])
        elif '# rampSteps #' in line: sim.SHEET[sheetID].rampSteps = int(line.split()[0])
        elif '# rampRelax #' in line: sim.SHEET[sheetID].rampRelax = int(line.split()[0])
        elif '# dzRamp #' in line: sim.SHEET[sheetID].dzRamp = float(line.split()[0])
        elif '# dpRamp #' in line: sim.SHEET[sheetID].dpRamp = float(line.split()[0])      
        # control: lateral
        elif '# fLateral #' in line: sim.SHEET[sheetID].fLateral = int(line.split()[0])
        elif '# vX #' in line: sim.SHEET[sheetID].vX = float(line.split()[0])
        elif '# vY #' in line: sim.SHEET[sheetID].vY = float(line.split()[0])
        elif '# frictRelax #' in line: sim.SHEET[sheetID].frictRelax = int(line.split()[0])      

        # OnSite
        elif '# fOnSitePotential #' in line: sim.SHEET[sheetID].fOnSitePotential = int(line.split()[0])
        elif '# fOnSitePeriod #' in line: sim.SHEET[sheetID].fOnSitePeriod = float(line.split()[0])
        elif '# gammaOnSite #' in line: sim.SHEET[sheetID].gammaOnSite = float(line.split()[0])
        elif '# frictionCoeffOS #' in line: sim.SHEET[sheetID].frictionCoeffOS = float(line.split()[0])
        elif '# vXOnSite #' in line: sim.SHEET[sheetID].vXOnSite = float(line.split()[0])
        elif '# vYOnSite #' in line: sim.SHEET[sheetID].vYOnSite = float(line.split()[0])

        # observables
        elif ('# resolMovie #' in line) and sheetID != 666: sim.SHEET[sheetID].resolMovie = int(line.split()[0])
        elif '# f3dMovie #' in line: sim.SHEET[sheetID].f3dMovie = int(line.split()[0])      



        # --- interSheet --- #
        
        elif '# inter start' in line: interID = int(line.split()[0])
        elif '# inter end' in line: interID = 666

        # basic
        elif '# sheetID0 #' in line: sim.INTER[interID].sheetID0 = int(line.split()[0])
        elif '# sheetID1 #' in line: sim.INTER[interID].sheetID1 = int(line.split()[0])
        # interactions
        elif '# fConstraint #' in line: sim.INTER[interID].fConstraint = int(line.split()[0])
        elif '# fPotential #' in line: sim.INTER[interID].fPotential = int(line.split()[0])
        elif '# surfEnerg #' in line: sim.INTER[interID].surfEnerg = float(line.split()[0])
        elif '# potRange #' in line: sim.INTER[interID].potRange = float(line.split()[0])
        elif '# potCurveRel #' in line: sim.INTER[interID].potCurveRel = float(line.split()[0])
        elif '# fPotentialTest #' in line: sim.INTER[interID].fPotentialTest = int(line.split()[0])
        elif '# nTimeOff #' in line: sim.INTER[interID].nTimeOff = int(line.split()[0])
        elif '# frictionCoeff #' in line: sim.INTER[interID].frictionCoeff = float(line.split()[0])
        # observables
        elif '# fDumpGap #' in line: sim.INTER[interID].fDumpGap = int(line.split()[0])
        elif '# fDumpFrame #' in line: sim.INTER[interID].fDumpFrame = int(line.split()[0])
        elif '# fDumpLateral #' in line: sim.INTER[interID].fDumpLateral = int(line.split()[0])
        elif ('# resolMovie #' in line) and interID != 666: sim.INTER[interID].resolMovie = int(line.split()[0])

      except:
        if (interID != 666) or (sheetID != 666):
          msg = '[WARNING] cannot process line %4i (' % iLine
          msg += ' sheetID='+str(sheetID)+'/'+str(len(sim.SHEET))
          msg += ' interID='+str(interID)+'/'+str(len(sim.INTER))
          print(msg+' ):')
          print(' '+line[:-1])

  # update defaults
  for iSheet in range(sim.nSheet):
    if sim.SHEET[iSheet].rYhertz == -1: sim.SHEET[iSheet].rYhertz = sim.SHEET[iSheet].rXhertz
  if sim.lengthY == 0: sim.lengthY = sim.lengthX
  if sim.nyGlobal == 0: sim.nyGlobal = sim.nxGlobal
  
  return sim