#!/usr/bin/env python

import time, os, sys, shutil
import numpy as np
import cPickle as pickle

cfgfile = 'parameters'
from parameters import *
from micro_meshstr import microMesh, microMesh_nonsynaptic
from mainMicro import *
from spiketrain import spiketrain
from initialize import initialize



def mainSynapse(terminal, multiple, b_conc, train):

  res = 0.025

  mMeshes = {} 
  if terminal == 'synaptic':
    mMeshes['synaptic'] = microMesh(dx1, dx2, distSYN, radSC, radOUT)
  elif terminal == 'non-synaptic':
    mMeshes['non-synaptic'] = microMesh_nonsynaptic(dx1, dx2, distSYN, radSC, radOUT)
  elif terminal == 'both':
    mMeshes['synaptic'] = microMesh(dx1, dx2, distSYN, radSC, radOUT)
    mMeshes['non-synaptic'] = microMesh_nonsynaptic(dx1, dx2, distSYN, radSC, radOUT)
  else:
    print ('ERROR in function: mainSynapse(terminal, multiple, b_conc, train)')
    print ('variable "terminal" should be "synaptic", "non-synaptic" or "both"')
    sys.exit()
    
  startMacro = time.time()
  os.mkdir('%s'%startMacro)
  target = open('%s/parameters' %(startMacro), 'w')
  target.write('baseline conc. = %f\n' %(b_conc))
  target.write('time dt = %f and space resolution res = %f\n' %(dt_macro, res))
  target.write('vf = %f\n' %vf)
  target.write('DA molecules = %f\n' %N)
  target.write('DATvar = %s\n' %DATvar)
  target.write('MMvar = %f' %MMvar)
  target.close()

  if multiple == True:
    if train == False: 
      ST = spiketrain(30, 1.)	
      nMolecules = np.random.normal(N*0.1, N*0.1*N_var, len(ST))
      np.savetxt('%s/spiketrain.csv' %(startMacro), ST)
      np.savetxt('%s/nMolecules.csv' %(startMacro), nMolecules)
    else:
      ST = train
  elif multiple == False:
    nMolecules = np.random.normal(N, N*N_var, 1)
  else:
    print ('ERROR in function: mainSynapse(terminal, multiple, b_conc, train)')
    print ('variable "multiple" should be "True" or "False"')
    sys.exit()

  for key in mMeshes.keys():

    mMesh = mMeshes[key]
    os.mkdir('%s/data_%s'%(startMacro, key))
    os.mkdir('%s/tmp'%startMacro)

    target = open('%s/data_%s/mesh.pkl' %(startMacro, key), 'w')
    pickler = pickle.Pickler(target, -1)
    pickler.dump(mMesh)
    target.close()

    if multiple == False:
      for upt in list([0,1]):
        phi = CellVariable(name = "concentration", mesh = mMesh, value = b_conc)
        phi_initial, initialMassSYN, volSYN = initialize(phi, nMolecules[0], vf)[:3]
        extUpt = upt
        print 'initialize vesicular release from %s terminal...' %key, phi_initial.max(), phi_initial.min(), extUpt
        conc, tmp1, tmp2, tmp3, mass, vol = mainMicro(phi_initial, startMacro, res, microMaxT, key, initialMassSYN, extUpt)
        np.savetxt('%s/data_%s/conc%s.csv' %(startMacro, key, upt), conc)
        np.savetxt('%s/data_%s/vol%s.csv' %(startMacro, key, upt), vol)
    else:
      os.mkdir('%s/data_%s/multipleRelease'%(startMacro, key))
      for upt in list([0,1]):
        phi = CellVariable(name = "concentration", mesh = mMesh, value = b_conc)
        phi_micro, initialMassSYN, volSYN = initialize(phi, nMolecules[0], vf)[:3]
        extUpt = upt
        print 'initialize vesicular release from %s terminal...' %key, phi_micro.max(), phi_micro.min(), extUpt
        for st in xrange(len(ST)-1):
          dT = ST[st+1]-ST[st]
          conc, phi, massEXT, massSYN, mass, vol = mainMicro(phi_micro, startMacro, res, dT, key, initialMassSYN, extUpt)
          np.savetxt('%s/data_%s/multipleRelease/conc%s%s.csv' %(startMacro, key, upt, ST[st]), conc)
          phi_micro = initialize(phi, nMolecules[st+1], vf)[0]
        conc = mainMicro(phi_micro, startMacro, res, microMaxT, key, massSYN, extUpt)[0]
        np.savetxt('%s/data_%s/multipleRelease/conc%s%s.csv' %(startMacro, key, upt, ST[-1]), conc)

    shutil.rmtree('%s/tmp'%startMacro)
  return 



