#!/usr/bin/env python

# medium/volume parameters

vf = 0.21 		# volume fraction = 0.21
tor = 1.54 		# tortuosity = 1.54
D = 763./tor/tor 	# free Diffusioncoefficient (= 763.)/tortuosity**2 in microm^2/s, see Nicholson 1995
avoga = 6.02214129e17 	# in microMol
rho = 0.104 		# density of striatal terminals in 1 microm^3 from Dreyer et al 2010
volAxon = 0.54e9 	# volume of axon in microM^^3

# mesh parameters

distSYN = 0.015 	# distance pre-post synapse (0.015 micrometer, Garris 1994)
radSC = 0.15 		# radius synapse (0.15 micrometer, Garris et al. 1994)
radOUT = 30. 		# 10. default radius for outer surface in micro model
radSPACE = 50.0 	# radius/edge length for whole simulation space (default 50.0)
radRelease = radSC/10
dx0 = distSYN/10
dx1 = distSYN/5
dx2 = radOUT/10 
dx3 = 2.

# time parameters

dt_macro = 0.00025	# dt_macro = 0.00025
dt_nano = 0.00000001	# dt_nano = 0.00000001
alpha = 1.6		# regulates how big timesteps are when initializing a synaptic DA release
microMaxT = 0.1

# parameters of interest for different scenarios of DA transmission

fRateTonic = 4.0 	# impulses per second in population firing (default = 4.0, from Dreyer et al 2010)
N = 3000. 		# number of dopamine molecules released during vesicle fusion (default = 3000, Pothos 1998)
N_var = 0.		# variance in DA molecules released
rprob = 0.06	 	# release probability (default = 0.06, from Dreyer et al 2010)
Vmax = 4.1	 	# default = 4.1, from Dreyer et al 2010 // in Wu et al 2001: Vmax=3.8 microM/s
Km = 0.4		# default = 0.21, from Dreyer et al 2010 // in Wu et al 2001: Km=0.22microM
DATvar = "none" 	# location of DAT: "none" (default), "out"=outside  (default for inhomogenous uptake), 
			# "in"=inside, "rim"= border of junction, "all"=whole terminal
MMvar = 1.		# to adjust surface uptake
bc = 0.01 		# initial baseline concentration macro model (microMolar = microMol/L), 10-50nM according to Grace(1994)
non_prop = 0.65 	# proportion of nonsynaptic terminals
PD_frac = 1.0   	# fraction of remained terminals in PD (default = 1.0)
rho = 0.104 * PD_frac

