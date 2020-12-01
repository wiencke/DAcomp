import numpy as np
import cPickle as pickle
import random
from spiketrain import *
from scipy.stats import bernoulli

# ----------- creating neuron population --------------

def neuroPop(mesh, rho, propTonic, propMixed, non_prop, burstD, burstM, dipD, start, relax, startMacro): 	
  nNeur = int(rho * mesh.cellVolumes.sum())				# number of terminals
  if mesh.cellVolumes.shape[0]<nNeur:
    print 'there are more terminals than volume units'
  locNeur = random.sample(range(0, mesh.numberOfCells), nNeur)		# terminal location
  nonNeur = random.sample(locNeur, int(non_prop*nNeur))			# non-synaptic terminals
  np.savetxt('%s/data/nonNeur.txt'%startMacro, nonNeur)
  np.savetxt('%s/data/location.txt'%startMacro, locNeur)
  print 'number of terminals: ', len(locNeur), ' number of non-synaptic terminals: ', len(nonNeur)

  neuroMixed = set(random.sample(set(locNeur), int(nNeur * propMixed)))						# tonic and phasic neurons
  neuroTonic = set(random.sample(set(locNeur)-neuroMixed, int(nNeur * propTonic)))				# tonic only neurons
  neuroBurst = set(random.sample(set(locNeur)-neuroMixed-neuroTonic, int(nNeur * (1-propTonic-propMixed))))	# bursting only neurons
  print 'terminals in tonic mode: ', len(neuroTonic), ' terminals in burst mode: ', len(neuroBurst),  ' terminals in mixed mode: ', len(neuroMixed)

  tmax = start+(np.array(burstD)+relax).sum()*len(burstM) + (np.array(dipD)+relax).sum()
  trains = {}
  for i in neuroTonic:
    train = spiketrain(fRateTonic,tmax)
    trains[i] = np.array(np.round(np.array(train)/dt_macro), int).tolist()  
  for i in neuroBurst:
    trainburst = list()
    trial = start
    for bd in burstD:
      for bm in burstM:
        trainburst = np.hstack((trainburst,(spiketrain(bm,bd)+trial*np.ones((1)))))
        trial +=relax + bd
    trains[i] = np.array(np.round(np.array(trainburst)/dt_macro), int).tolist()
  for i in neuroMixed:
    trainburst = spiketrain(fRateTonic,start)
    trial = start
    for bd in burstD:
      for bm in burstM:
        trainburst = np.hstack((trainburst,(spiketrain(bm,bd)+trial*np.ones((1)))))
        trial += bd
        trainburst = np.hstack((trainburst,(spiketrain(fRateTonic,relax)+trial*np.ones((1)))))
        trial +=relax
      for dip in dipD:
        trial += dip
        trainburst = np.hstack((trainburst,(spiketrain(fRateTonic,relax)+trial*np.ones((1)))))
        trial +=relax      
    trains[i] = np.array(np.round(np.array(trainburst)/dt_macro), int).tolist()
  for neuron in trains.keys():
    removeN = len(trains[neuron]) - (bernoulli.rvs(rprob, size = len(trains[neuron]))).sum()
    for k in xrange(removeN):
      trains[neuron].remove(random.sample(trains[neuron],1)[0])
  target = open('%s/data/trains.pkl'%startMacro, 'w')
  pickler = pickle.Pickler(target, -1)
  pickler.dump(trains)
  target.close()

  return nonNeur, trains

# ----------- monitor macro model --------------

def rOccSpace(phi):
  MMesh = phi.mesh
  VOL = MMesh.cellVolumes
  interiorCells = np.ones(len(VOL), bool)
  for id in MMesh.faceCellIDs[0][MMesh.exteriorFaces()].compressed():
    interiorCells[id] = False
  for id in MMesh.faceCellIDs[1][MMesh.exteriorFaces()].compressed():
    interiorCells[id] = False
  phiBar = phi()*interiorCells
  VOLbar = VOL*interiorCells
  phiVolBar = np.dot(VOL,phiBar)
  phi_mean = np.dot(VOL,phiBar)/VOLbar.sum()
  phi_var = np.sqrt(np.dot(VOL,(phiBar-phi_mean)**2)/VOLbar.sum())
  return  phi_mean, phi_var, np.dot(VOL,phiBar/(phiBar+1.)).sum()/VOLbar.sum(), np.dot(VOL,phiBar/(phiBar+0.01)).sum() /VOLbar.sum()

    
# ----------- help functions for synaptic transmission analysis --------------

def conc_DTplane(phi_micro, dmax, res):
  mMesh = phi_micro.mesh
  xc, yc, zc = mMesh.cellCenters()
  VOL = mMesh.cellVolumes

  conc = np.zeros((dmax))
  mass = np.zeros((dmax))
  vol = np.zeros((dmax))
  for i in xrange(dmax):  
    Cells = (np.sqrt(xc**2 + yc**2 + zc**2)<(i+1)*res) * (np.sqrt(xc**2 + yc**2+ zc**2)>=i*res)
    conc[i] = np.dot(VOL,phi_micro()*Cells)/(VOL*Cells).sum()
    mass[i] = np.dot(VOL,phi_micro()*Cells)
    vol[i]  = (VOL*Cells).sum()
  print 'min conc: ', conc.min(), 'max conc: ', conc.max()
  return conc, mass, vol

# ----------- clear memory --------------

def cleanSpace1(startMacro, T, phi_micro):
      target = open('%s/tmp/phi%f.pkl' %(startMacro,T), 'w')
      pickler = pickle.Pickler(target, -1)
      pickler.dump(phi_micro)
      target.close()

def cleanSpace2(startMacro, T):
      target = open('%s/tmp/phi%f.pkl' %(startMacro,T), 'rb')
      phi_micro = pickle.load(target)
      target.close()   
      return phi_micro


