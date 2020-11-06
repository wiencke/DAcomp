#!/usr/bin/env python

from fipy import *
import numpy as np

cfgfile = 'parameters'
from parameters import *


def microModel(phi, dt, key, extUpt, constrain):


  mesh = phi.mesh
  dist = distSYN
  rad = radSC
  dx = dx0

  xc, yc, zc = mesh.cellCenters()
  xfc, yfc, zfc = mesh.faceCenters()  
  VOL = mesh.cellVolumes

  
  # ________cells/faces where Michaelis-Menten-Reuptake takes place________________________________

  radXY =  np.sqrt(xfc**2 + yfc**2)

  if DATvar == "out":
    MMFaces = (zfc >= dist/2 + dx) * (zfc <= radSC + dx)    * (radXY <= radSC + dx) * mesh.exteriorFaces() * (radXY >= radSC - 2*dx)
  elif DATvar == "in":
    MMFaces = (zfc >= dist/2 - dx) * (zfc <= dist/2 + dx)   * (radXY <= radSC + dx) * mesh.exteriorFaces()
  elif DATvar == "rim":
    MMFaces = (zfc >= dist/2 - dx) * (zfc <= dist/2 + 2*dx) * (radXY <= radSC + dx) * mesh.exteriorFaces()
  elif DATvar == "all":
    MMFaces = (zfc >= dist/2 - dx) * (zfc <= radSC + dx)    * (radXY <= radSC + dx) * mesh.exteriorFaces()
  elif DATvar == "none":
    MMFaces = (zfc >= dist/2 + dx) * (zfc <= dist/2 - dx)     
  else: print 'ERROR, wrong DATvar Input'                  

 
  # ________calculations___________________________________________________________________________

  phi_old = np.array(phi())

  MMValue = np.array(Vmax*MMvar*phi.faceValue/(Km/MMvar+phi.faceValue)())
  IdValue = phi.faceValue()
  MMGradient = (MMValue <= IdValue)				
  IdGradient = (MMValue >  IdValue)


  phi.faceGrad.constrain(-MMValue*mesh.faceNormals, where = MMGradient * MMFaces)
  phi.faceGrad.constrain(-IdValue*mesh.faceNormals, where = IdGradient * MMFaces)

  if not constrain == False:
    phi.constrain(constrain, where = mesh.exteriorFaces*(np.sqrt(mesh.faceCenters()**2).sum(axis=0)>radSC+dist))
  if key == 'synaptic':    
    exSynCells = (np.sqrt(xc**2 + yc**2 + zc**2) >= rad)
  elif key == 'non-synaptic': 
    exSynCells = np.abs(zc) >= 0
  MichMen = Vmax*phi/(Km+phi) * np.array(exSynCells, 'int')
  
  eq = TransientTerm() == DiffusionTerm(coeff=(D))  - extUpt* MichMen/vf 
  eq.solve(var=phi, dt=dt)	

  if phi().min() < 0:
    print 'Warning!! negative concentration:'
    mask = (phi() < 0)
    phi.setValue(0, where = mask)
  massNEW   = np.dot(VOL, phi()) 
  uptake    = np.dot(VOL, phi()-phi_old)
  massEXTRA = np.dot(VOL, (phi()-phi_old) * exSynCells)
  massSYN   = massNEW - np.dot(VOL, phi() * exSynCells)
    
  return phi, massEXTRA, massSYN, uptake
    
