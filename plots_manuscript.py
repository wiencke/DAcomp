import numpy as np
from matplotlib import pyplot as plt
plt.ion()
from matplotlib.colors import LogNorm
from parameters import *
from matplotlib import ticker, cm




def concSN(dirName, res, dt, upt):
    concSyn = np.loadtxt('%s/data_synaptic/conc%s.csv'%(dirName, upt))
    concNon = np.loadtxt('%s/data_non-synaptic/conc%s.csv'%(dirName, upt))
    x=np.arange(0,concSyn.shape[1])*res
    y=np.arange(0,concSyn.shape[0])*dt
    levels = np.array((0.01, 0.1, 1.))

    # Fig1: concentration on log scale
    plt.figure()
    plt.pcolor(y, x, concSyn.T, vmin=0.0001, vmax=10., norm=LogNorm())
    plt.colorbar(); plt.title('conc syn. release')
    plt.contour(y,x,concSyn.T,levels, colors='white')
    
    # Fig2: concentration difference between synaptic and non-synaptic terminal 
    #plt.figure(); plt.pcolor(x, y, concNon-concSyn, vmin= 0.001, vmax= 1000., cmap= 'Blues', norm=LogNorm()); plt.colorbar()
    #plt.pcolor(x, y, concSyn-concNon, vmin= 0.001, vmax= 1000., cmap = 'Oranges', norm=LogNorm()); plt.colorbar()
    #plt.title('conc diff: synaptic(orange) vs. non-synaptic(blue)')

    # Fig3: dynamic fast binding
    dInd = list([4, 8, 20, 80, 160])
    color = list(['b', 'orange', 'g', 'r', 'k'])
    D2_tot = 100
    B0 = 0.
    D2Syn = np.ones((concSyn.shape))*B0
    D2Non = np.ones((concNon.shape))*B0
    plt.figure();
    for j in xrange(D2Syn.shape[1]):
        for t in xrange(D2Syn.shape[0]-1):
          D2Syn[t+1,j] = np.min((D2Syn[t,j]+dt*( 0.2772*1000*concSyn[t,j]*(D2_tot-D2Syn[t,j])-6.93*D2Syn[t,j]),D2_tot)) 
          D2Non[t+1,j] = np.min((D2Non[t,j]+dt*( 0.2772*1000*concNon[t,j]*(D2_tot-D2Non[t,j])-6.93*D2Non[t,j]),D2_tot)) 

    for i in xrange(len(dInd)):
        plt.plot(y, (D2Syn)[:,dInd[i]], color=color[i], label='distance = %s'%(res*dInd[i]))
    for i in xrange(len(dInd)):
        plt.plot(y, (D2Non)[:,dInd[i]], ':',color=color[i])
    plt.legend()
    plt.title('fast D2R binding')

    # Fig4: dynamic slow binding
    dInd = list([4, 8, 20, 80, 160])
    color = list(['b', 'orange', 'g', 'r', 'k'])
    D2_tot = 100
    B0 = 0.
    D2Syn = np.ones((concSyn.shape))*B0
    D2Non = np.ones((concNon.shape))*B0
    plt.figure();
    for j in xrange(D2Syn.shape[1]):
        for t in xrange(D2Syn.shape[0]-1):
          D2Syn[t+1,j] = np.min((D2Syn[t,j]+dt*( 0.000333*1000*concSyn[t,j]*(D2_tot-D2Syn[t,j])-0.00833*D2Syn[t,j]),D2_tot)) 
          D2Non[t+1,j] = np.min((D2Non[t,j]+dt*( 0.000333*1000*concNon[t,j]*(D2_tot-D2Non[t,j])-0.00833*D2Non[t,j]),D2_tot)) 

    for i in xrange(len(dInd)):
        plt.plot(y, (D2Syn)[:,dInd[i]], color=color[i], label='distance = %s'%(res*dInd[i]))
    for i in xrange(len(dInd)):
        plt.plot(y, (D2Non)[:,dInd[i]], ':',color=color[i])
    plt.legend()
    plt.title('slow D2R binding')

    return

def inhomUPT(dirName1, dirName2, dirNamie3, res, dt, upt):

    concSyn1 = np.loadtxt('%s/data_synaptic/conc%s.csv'%(dirName1, upt))
    concSyn2 = np.loadtxt('%s/data_synaptic/conc%s.csv'%(dirName2, upt))
    concSyn3 = np.loadtxt('%s/data_synaptic/conc%s.csv'%(dirName2, upt))
    x=np.arange(0,concSyn1.shape[1])*res
    y=np.arange(0,concSyn1.shape[0])*dt
    levels = np.array((0.01, 0.1, 1.))
    plt.figure()
    plt.pcolor(y, x, concSyn1.T, vmin=0.0001, vmax=10., norm=LogNorm())
    plt.colorbar(); plt.title('inhomogeneous uptake')
    plt.contour(y,x,concSyn1.T,levels, colors='white')
    plt.contour(y,x,concSyn2.T,levels, colors='gray')
    plt.contour(y,x,concSyn3.T,levels, colors='black')


def multipleSpiking(dirName, res, dt, upt):
    ST = np.loadtxt('%s/spiketrain.csv'%dirName)[:-1]
    concSyn = {}
    for st in ST:        
      concSyn[st] = np.loadtxt('%s/data_synaptic/multipleRelease/conc%s%s.csv'%(dirName, upt, st))
    concSynALL = concSyn[ST[0]]
    for s in xrange(len(concSyn.keys())-1):
      concSynALL = np.vstack((concSynALL[:int(np.ceil((ST[s+1]-ST[0])/ dt_macro)),:], concSyn[ST[s+1]]))
    y=np.arange(0,concSynALL.shape[1])*res
    x=np.arange(0,concSynALL.shape[0])*dt*1000
    plt.figure()
    
    # Fig1: concentrations with isolines conc = 0.01, 0.1 and 1.0
    plt.subplot(4,1,1); 
    levels = np.array((0.01, 0.1, 1.))
    plt.pcolor(x, y, concSynALL.T, vmin=0.000001, vmax=10., norm=LogNorm()); plt.colorbar()
    plt.contour(x,y,concSynALL.T,levels, colors='white')
    plt.ylabel('distance [micro meter]')
    plt.axis((-40,1040,0,10))

    # Fig2:  D2R equilibrium binding
    plt.subplot(4,1,2); 
    plt.pcolor(x, y, concSynALL.T/(concSynALL.T+0.01), vmin=0.0, vmax=1., cmap = 'Oranges'); plt.colorbar()
    plt.axis((-40,1040,0,10))

    # Fig3: slow dynamic D2R binding at specific distances
    D2_tot = 100
    dInd = list([4, 8, 20, 80, 160])
    D2Syn = np.ones((concSynALL.shape))*0.
    for j in xrange(D2Syn.shape[1]):
        for t in xrange(D2Syn.shape[0]-1):
          D2Syn[t+1,j] = D2Syn[t,j]+dt*(0.000333*1000*concSynALL[t,j]*(D2_tot-D2Syn[t,j])-0.00833*D2Syn[t,j]) # 0.000333 and 0.00833 (slow) vs 0.2772 and 6.93 (fast)
    plt.subplot(4,1,3); 
    for i in dInd:
       plt.plot(x, D2Syn[:,i], label='dist = %s'%(res*i)) 
    plt.legend()

    # Fig4: fast dynamic D2R binding at specific distances
    D2Syn = np.ones((concSynALL.shape))*0.
    for j in xrange(D2Syn.shape[1]):
        for t in xrange(D2Syn.shape[0]-1):
          D2Syn[t+1,j] = np.min((D2Syn[t,j]+dt*(0.2772*1000*concSynALL[t,j]*(D2_tot-D2Syn[t,j])-6.93*D2Syn[t,j]),D2_tot)) # 0.000333 and 0.00833 (slow) vs 0.2772 and 6.93 (fast)
    plt.subplot(4,1,4); 
    for i in dInd:
       plt.plot(x, D2Syn[:,i], label='dist = %s'%(res*i)) 

    return


