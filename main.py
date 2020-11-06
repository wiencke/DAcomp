#!/usr/bin/env python

from fipy import *
import numpy as np
import cPickle as pickle
import time, os
from matplotlib import pyplot as plt
plt.ion()

cfgfile = 'parameters'
from parameters import *
from macro_meshstr import *
from helpFunctions import *

# **************** mesh ******************

startMacro = time.time()
os.mkdir('%s'%startMacro)
os.mkdir('%s/data'%startMacro)
os.mkdir('%s/tmp'%startMacro)
os.system('cp parameters.py %s/data/parameters'%startMacro) 

MMesh, MMeshstyle = macroMeshRect(dx3, radSPACE)
print 'macro mesh constructed'
VOL = MMesh.cellVolumes 

target = open('%s/data/MMesh.pkl'%startMacro, 'w')
pickler = pickle.Pickler(target, -1)
pickler.dump(MMesh)
target.close()

# **************** neuron population ******************

burstT = []				
burstM = []
dipT = 0		
start = 5.0
relax = 0.0
propTonic = 0.0
propMixed = 1.0

print 'configuring neuron population...'
nonNeur, trains = neuroPop(MMesh, rho, propTonic, propMixed, non_prop, burstT, burstM, dipT, start, relax, startMacro)

# **************** configure simulation ***************

maxIteration = int((start+(np.array(burstT)+relax).sum() + (np.array(dipT)+relax).sum())/dt_macro) +1
concCell = {}
concHood = {}
for neuron in trains.keys():
  concCell[neuron] = np.zeros((maxIteration))
  concHood[neuron] = np.zeros((maxIteration))
concCellSyn={}
concCellNon={}
concHoodSyn={}
concHoodNon={}
nonSample = random.sample(nonNeur, int(len(nonNeur)*0.05))
synSample = random.sample(set(trains.keys()).difference(set(nonNeur)), int((len(trains.keys())-len(nonNeur))*0.05))

monSpace = []
releaseMass = random.gauss(N,N_var)/avoga/1e-15
massLS = releaseMass*np.loadtxt('MassLeavingSynapse.csv')

lastSpike = {}
for neuron in trains.keys():
  lastSpike[neuron] = []

phi_macro = CellVariable(name = "concentration", mesh = MMesh, value = bc)
print 'start simulation:', ' conc:', np.dot(VOL, phi_macro())/VOL.sum()

j=0
while j < maxIteration:	

  for neuron in trains.keys():
    concCell[neuron][j] = phi_macro()[neuron]
    concHood[neuron][j] = (phi_macro.faceValue()[MMesh.cellFaceIDs.data[:,neuron]]).mean()
    if len(trains[neuron])>0:
      nextSpike = trains[neuron][0]
      mask = np.zeros(len(VOL), bool)
      mask[neuron] = True
      bc = phi_macro()[neuron]
      vol = MMesh.cellVolumes[neuron]
      if neuron in nonNeur:
        if nextSpike == j:
          print 'nonsynaptic release', len(trains[neuron])
          mass = releaseMass + bc*vol*vf
          phi_macro.setValue(mass/vol/vf, where = mask)    
          trains[neuron].remove(nextSpike)    
      else:
        if nextSpike == j:
          print 'synaptic terminal release'
          if not lastSpike[neuron]:
            lastSpike[neuron] = list(massLS)			
          else:
            ls = 0
            while ls < len(lastSpike[neuron]):
              lastSpike[neuron][ls]+= massLS[ls]; ls+=1
            while ls < len(massLS):
              lastSpike[neuron].append(massLS[ls]); ls+=1
          mass = bc*vol*vf + lastSpike[neuron][0]			
          phi_macro.setValue(mass/vol/vf, where = mask)
          del lastSpike[neuron][0]	
          trains[neuron].remove(nextSpike)
        else: 
          if lastSpike[neuron]:
            mass = bc*vol*vf + lastSpike[neuron][0]
            del lastSpike[neuron][0]	
            phi_macro.setValue(mass/vol/vf, where = mask)        


  MichMen = Vmax*phi_macro/(Km+phi_macro)
  eq = TransientTerm() == DiffusionTerm(coeff=D) - MichMen/vf
  eq.solve(var=phi_macro, dt=dt_macro)	
  monSpace.append(rOccSpace(phi_macro))		
  if np.mod(j,20) == 0:
    print 'iteration time [s]:', j*dt_macro, ' computation time [min]', (time.time() - startMacro)/60
    print 'conc:', np.dot(VOL, phi_macro())/VOL.sum()
  j+= 1
  target = open('%s/tmp/phimacrotmp.pkl'%startMacro, 'w')
  pickler = pickle.Pickler(target, -1)
  pickler.dump(phi_macro)
  target.close()
  del phi_macro, eq, MichMen
  target = open('%s/tmp/phimacrotmp.pkl'%startMacro, 'rb')
  phi_macro = pickle.load(target)
  target.close()

  if np.mod(j,1000) == 0:
    np.savetxt('%s/data/mean_concentration.csv'%(startMacro), np.array(monSpace))
    for neuron in nonSample:
      concCellNon[neuron] = concCell[neuron]
      concHoodNon[neuron] = concHood[neuron]
    for neuron in synSample:
      concCellSyn[neuron] = concCell[neuron]
      concHoodSyn[neuron] = concHood[neuron]
    target = open('%s/data/synNeur.pkl'%startMacro, 'w')
    pickler = pickle.Pickler(target, -1)
    pickler.dump(concCellSyn)
    target.close()
    target = open('%s/data/nonNeur.pkl'%startMacro, 'w')
    pickler = pickle.Pickler(target, -1)
    pickler.dump(concCellNon)
    target.close()
    target = open('%s/data/synHood.pkl'%startMacro, 'w')
    pickler = pickle.Pickler(target, -1)
    pickler.dump(concHoodSyn)
    target.close()
    target = open('%s/data/nonHood.pkl'%startMacro, 'w')
    pickler = pickle.Pickler(target, -1)
    pickler.dump(concHoodNon)
    target.close()
    'variables have been saved...'

print 'macro simulation duration:', (time.time()-startMacro)/60/60 , 'hours'


monSpaceArray = np.array(monSpace)
c = monSpaceArray[:,0]/vf
plt.plot(np.arange(len(c))*dt_macro, c)

crop = int(np.ceil(len(c)/4.))
mean = c[crop:].mean()
plt.plot(np.arange(len(c)-crop)*dt_macro+crop*dt_macro, np.ones(len(c)-crop)*mean, label = '%f'%(mean*1000) )
plt.legend()




