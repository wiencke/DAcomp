# DAcomp

This repository includes code for simulating DA transmission. See the corresponding manuscript 'Dopamine release, diffusion and uptake: A computational model for synaptic and volume transmission' for more information.

## Set Parameters

To change the hard-coded parameters edit the `parameters.py` file. 

## Generate Figures

### Single release events: Fig. 3

This Figure in the above mentioned manuscript illustrates single DA release events from a (synaptic) terminal. To generate plots A-D you can e.g. type in your ipython shell:
```
%run plots_manuscript
%run main_terminal 
mainSynapse('both', False, bc, False)
```
The variable `bc` should be set to bc = 0.004 for plots A, C, E-H or bc = 0.05 for plots B & D. Alternatively, this variable can be set accordingly in the `parameters.py` file. By default `mainSynapse()` generates a folder that is called according to the unixtime when running the command. The folder contains two subfolders, `data_synaptic` and `data_non-synaptic`, and a file called `parameters`, which lists some of the parameters used for this simulation and that are important to remember, for example to run the following command with the correct values: 

```
concSN('unixtime', res, dt, upt)
```

For our plots `0.025` and `dt = 0.00025`. The value for the uptake parameter should either be `upt = 0` if you wish to generate plots A, B & F (simulations without uptake) or `upt = 1` for plots C & D (simulations with uptake). To create plot E you need to run `mainSynapse()` two more times whereby you need to adjust the DA transporter variable `DATvar` and the Michaelis-Menten variable `MMvar` in the `parameters.py` file. `DATvar` is used to describe where DA transportes are located. By default this variable is set to 'none', when only volumetric uptake is simulated. To simulate uptake occuring on the terminal membrane surface we set `DATvar = out`, such that uptake takes place only outside the synaptic cleft. `MMvar` is the factor to increase surface uptake. Thus, by default this value equals one. We ran a second simulation with `MMvar = ` which created a folder named `unixtime2` and a second simulation with `MMvar = ` which created a folder named `unixtime3`. Plot E can then be generated with the following command:

```
inhomUPT('unixtime', 'unixtime2', 'unixtime3', res, dt, 1):
```
### Multiple release events: Fig. 4

This Figure in our manuscript 'Dopamine release, diffusion and uptake: A computational model for synaptic and volume transmission' illustrates multiple DA release events from a single terminal. To generate plots A-D you can e.g. type in your ipython shell:
```
%run plots_manuscript
%run main_terminal 
mainSynapse('synaptic', True, bc, False)
```

### Volume transmission: Fig. 5

To illustrate issues regarding volume transmisson we simulated three scenarios. First, we refer to the normal or healthy condition as the "default setting simulation" (DSS). To simulate the presence of an agent stimulating DA release (SRS), in file `parameters.py` we increased the release probability (variable `rprob`) from 6% to 7.8%. This elevates DA tone to 129% of the DSS. For the uptake inhibition scenario (UIS), in file `parameters.py`, we changed the variables `Vmax` to 80% of the default and `Km` to 24μM, which elicits an increase of ~ 1273% overbaseline in extrasynaptic DA. With the respective parameters in `parameters.py` for each of the different scenarios we ran the following commands:

```
%run main
```
Again, by default this creates folders called according to the unixtime when running the command, here we refer to these folders as `unixtimeDSS`, `unixtimeSRS` and `unixtimeUIS`.
