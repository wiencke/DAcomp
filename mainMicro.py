#!/usr/bin/env python

import cPickle as pickle
import time
import numpy as np
cfgfile = 'parameters'
from parameters import *

from micro import microModel
from helpFunctions import conc_DTplane, cleanSpace1, cleanSpace2


def mainMicro(phi, startMacro, res, maxT, key, initialMassSYN, extUpt):
  mMesh = phi.mesh
  xfc, yfc, zfc = mMesh.faceCenters() 
  xc, yc, zc = mMesh.cellCenters()
  VOL = mMesh.cellVolumes
  T = 0.
  dt = dt_nano
  massSYN = [initialMassSYN]
  massEXT = [0.]
  massUPT = [0.]

  if res == False:
    conc = np.array((0,0))
  else:
    tmax = int(maxT/dt_macro)
    dmax = int(radOUT/res)*3/4
    conc = np.zeros((tmax, dmax))
    mass = np.zeros((tmax, dmax))
    vol  = np.zeros((dmax))
    conc[0,:], mass[0,:], vol[:] = conc_DTplane(phi, dmax, res)
 
  while T+dt < dt_macro:
    phi = microModel(phi, dt, key, extUpt, False)[0]
    T += dt
    if np.mod(int(np.round(T/dt_nano)),10) == 0: 
      target = open('%s/tmp/phi.pkl' %(startMacro), 'w')
      pickler = pickle.Pickler(target, -1)
      pickler.dump(phi)
      target.close()
      del phi
      target = open('%s/tmp/phi.pkl' %(startMacro), 'rb')
      phi = pickle.load(target)
      target.close()   

    dt = alpha*dt

  timestart = time.time()

  dt = dt_macro - T					
  phi, EXT, SYN, UPT = microModel(phi, dt, key, extUpt, False)
  massSYN.append(SYN)
  massEXT.append(EXT)
  massUPT.append(UPT)

  # attention: if maxT <2*dt_macro the script gives an error because it tries to index conc
  # where index 1 is out of bounds for axis 0 with size 1

  T += dt
  if not res == False:
    conc[1,:], mass[1,:], vol[:] = conc_DTplane(phi, dmax, res)
  
  print 'DA molecules diffuse into extra-cellular space...'

  while T < maxT:
    phi, EXT, SYN, UPT = microModel(phi, dt_macro, key, extUpt, False)
    massSYN.append(SYN)
    massEXT.append(EXT)
    massUPT.append(UPT)
    T += dt_macro
    index = np.round(T/dt_macro)
    if not res == False and index < tmax:
      conc[int(index),:], mass[int(index),:], vol[:] = conc_DTplane(phi, dmax, res)
    if np.mod(index,10) == 0:
      cleanSpace1(startMacro, T, phi)
      del phi
      phi = cleanSpace2(startMacro, T)
      print '... %f %s done in %f minutes'%(T*100/maxT, '%', (time.time() - timestart)/60.)
  return conc, phi, massEXT, massSYN, mass, vol



