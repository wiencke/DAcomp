#!/usr/bin/env python

from fipy import *
import numpy as np
import random

cfgfile = 'parameters'
from parameters import *

dx = dx0
dist = distSYN
rad = radSC


def initialize(phiIN,N, vf):
  
  mesh = phiIN.mesh
  xc, yc, zc = mesh.cellCenters()
  xfc, yfc, zfc = mesh.faceCenters()  
  VOL = mesh.cellVolumes
  
  bottomFaces = (zfc >= dist/2 - dx) * mesh.exteriorFaces() * (zfc <= dist/2 + dx)
  releaseFaces = bottomFaces * (np.sqrt(xfc**2 + yfc**2) < radRelease)
  if not releaseFaces: print ('no release site')

  
  releaseCells = np.zeros(len(zc), bool)
  for id in mesh.faceCellIDs[0][releaseFaces].compressed():
    releaseCells[id] = True
  for id in mesh.faceCellIDs[1][releaseFaces].compressed():
    releaseCells[id] = True

  VOLreleaseCells = VOL*releaseCells

  MASSinCells = phiIN()*VOL*1e-15*vf
  releaseMass = random.gauss(N,N_var)/avoga # mass released during vesicle fusion
  releaseMASSinCells = releaseMass*VOLreleaseCells/VOLreleaseCells.sum()

  print releaseMASSinCells.sum(), ' mol released in ', VOLreleaseCells.sum(), 'qubic micro meter' 

  phiOUT = CellVariable(name = "concentration", mesh = mesh, value = (MASSinCells+releaseMASSinCells)/VOL/1e-15/vf) # transformation liter : cubic micro meter

  cleftCells = (np.sqrt(xc**2 + yc**2)<rad) * (np.abs(zc) < dist/2)
  VOLcleftCells = VOL*cleftCells
  initialMassSYN = np.dot(VOL,phiOUT()*cleftCells)
  initialMassEXT = np.dot(VOL, phiOUT()) - initialMassSYN

  return phiOUT, initialMassSYN, VOLcleftCells, initialMassEXT, VOL-VOLcleftCells
    
