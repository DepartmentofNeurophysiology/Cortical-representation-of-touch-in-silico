# Cortical-representation-of-touch-in-silico
Biologically inspired, computationally efficient network model of the somatosensory cortex after reconstructing the mouse barrel cortex in soma resolution and defining a mathematical model of cortical neurons whose action potential threshold adapts to the rate of ongoing network activity impinging onto the postsynaptic neuron. 
The model was used in the following two preprints:
* Huang, C., Englitz, B., Reznik, A., Zeldenrust, F., & Celikel, T. (2020). Information transfer and recovery for the sense of touch. BioRxiv, 2020.12.08.415729. https://doi.org/10.1101/2020.12.08.415729
* Huang, C., Zeldenrust, F., & Celikel, T. (2020). Cortical representation of touch in silico. Neuroinformatics. https://doi.org/10.1007/s12021-022-09576-5

Please cite the second paper when using this model.

A NetPyNE (http://www.netpyne.org/) implementation can be found in this repository: https://github.com/DepartmentofNeurophysiology/Cortical-representation-of-touch-in-silico-NetPyne

A few realizations of the connectivity calculations can be found in this data respository: https://doi.org/10.34973/tmf3-2m63


## About

Simulation of one or more barrel cortical columns. The model is a biologically inspired, computationally efficient network model of  somatosensory cortical columns. It was created by (1) reconstructing the barrel cortex in soma resolution using multi-channel mosaic scanning confocal microscopy, (2) defining a mathematical model (Izhikevich model [1]) of cortical neurons whose action potential threshold adapts to the rate of ongoing network activity impinging onto the postsynaptic neuron and (3) connecting each neuron in the network using statistical rules of pair-wise connectivity based on experimental observations. The input consists of whisker data (here: angle and curvature, but other metrics such as deviation from baseline are also possible). Based on this, thalamic spike trains are generated (thalamic neurons are considered simple filter-and-fire Poisson neurons), that then form the input to the cortical model. 

## Installation

Download files, preserving the folder structure. The model will run when this structure is added to the path (addpath(genpath([‘.’])), as is done in run_sim.m).

## General

A simple example simulation can be run by executing the script quick_example, which uses the function ‘run_sim’ (this function controls the simulation, and can be used as an example for doing user-defined simulations). A 3x1 grid of barrels in L4 and L23 will be constructed, with corresponding thalamic barreloids. Thalamic input spike trains will be generated in reaction to whisker input (in this case: whisker base angles and curvatures from the lab of Karel Svoboda [2], publicly available on Data Sets — CRCNS.org), and the cortical model will respond to these thalamic spike trains. 

The model includes Short Term Depression (STD) and Potentiation (STP) and (optional)
1. Direct whisker modulation by motor cortex [3].
2. Spike Timing-Dependent Plasticity (STDP) between L4 and L2/3 [4] between L4 neurons [5] and between L2/3 neurons [6]

## Input

As input, data files (whisker base angle and curvature) of the lab of Karel Svoboda [2] were used (publicly available on Data Sets — CRCNS.org) were used. The specified whisker data files are expected to be present in folder ‘Input data’, the example is written for file 'an171923_2012_06_04_data_struct.mat'. However, the file is too large to upload here. The function ‘make_thalamic_spike_trains_svoboda_recordings’ in folder ‘Make_New_Thinput’ makes spike trains from the whisker data. For the use of any other input, make a file similar to this example. 

## Simulation

A full simulation includes:
1.	generating new input spike trains (make_new_thalamic_input = 1), relevant functions in ‘Make_New_Thinput’
2.	making new thalamic filter neurons (make_new_thalamic_kernels = 1), relevant functions in ‘Make_New_Thinput’
3.	Make a new realisation of the network connectivity (make_new_connectivity=1), relevant functions in ‘Make_New_Connectivity’
4.  Choosing whether to turn on direct whisker modulution by motor cortex (includemodulationyn = 1)
5.  Choosing whether to turn on STDP (includeSTDPyn = 1)
6.  Choosing a name identifyer for the simulations (savename)
6.	Initialising and running the network, relevant functions in ‘Network Simulations’ (run_sim)

## Results

This will result in the following files in subfolder ‘Simulation Results -> (user defined name)

* cellinfo: information about all cells: (Number of Cells-by 6, note that L4 cells come first, followed by L23 cells) matrix, with
	* the first 3 columns are the location of each cell 
		* first column: position from the middle of the barrel (~rostral-caudal)
		* second column: position from the middle of the barrel (~dorsal-ventral) 
		* third column: depth, distance from pia
	* 4th column is the cell type 
		* 1 - L4 spiny stallet
		* 2 - L4 pyramidal
		* 3 - L4 Fast spike 
		* 4 - L4 low-threshold spike
		* 5 - L2/3 pyramidal 
		* 6 - L2/3 PV+ fast spike
		* 7 - L2/3 PV+ chandler
		* 8 - L2/3 PV+ bursting 
		* 9 - L2/3 SOM+ martinotti 
		* 10 - L2/3 SOM+ bitufted 
		* 11 - VIP+ double bouquet 
		* 12 - VIP+ bipolar
		* 13 - CR+ bipolar
		* 14 - CR+ multipolar/basket			
		* 15 - neurogliaform 
		Note: some L2/3 types have been merged with other types so in the matrix the numbers are 0); 
	* 5th column is the barrel identification (different number indication different barrels). 
	* 6th column is Exe/Inh indicator (1 for excitatory and -1 for inhibitory neurons)
* CMDMs: all connectivity data
* CMDMs_(…)_ConData
* CMDMs_(…)_ParaMat 
* CMDMs_(…)_ParaMat_reduced: a reduced version of the previous set, to speed up simulations
* CMDMs_(…)_WhiskerModulationModel: the model for the direct modulation by motor cortex [3]
* (…)_initialsettings: the initial settings to each simulation
* (…)_Thalamic_Kernels: thalamic filters for each thalamic neuron
* (…)_Thalamic_Spike_Trains: all input spike trains and whisker traces
* (…)_WhiskerModulation: the direct whisker modulation [3] for these input trials
* (…)_Simcolumn_(…)_(simulation number): spike times and membrane potentials of each simuation (what is saved depends on what is passed in variable 'whattosave' in run_sim)
* (…)_Simcolumn_(…)_(simulation number)_calciumdata: spike times convolved with an exponential kernel to mimic 2-photon calcium data

NB Note that the measured cell densities are in mice, but most axon/dendritic distribution patterns, connectivity and synaptic efficacy in the literature are from rats. To make these fit, the measured cell densities were scaled to fit rat data. So the resulting cell locations are appropriate for rat. Note that once the connectivity data are generated, these cell locations are not used in the network simulations (as neurons are simulated as point neurons). To get back to relevant mouse cell locations, use function cellinfo_mouse = rat_to_mouse_locations(cellinfo).

## References

[1] Izhikevich, E. M. (2003). Simple Model of Spiking Neurons. IEEE Transactions on Neural Networks, 14(6), 1572–1596. https://doi.org/10.1109/TNN.2003.820440

[2] Peron, S. P., Freeman, J., Iyer, V., Guo, C., & Svoboda, K. (2015). A Cellular Resolution Map of Barrel Cortex Activity during Tactile Behavior. Neuron, 86(3), 783–799. https://doi.org/10.1016/j.neuron.2015.03.027

[3] Crochet, S., Poulet, J. F. A., Kremer, Y., & Petersen, C. C. H. (2011). Synaptic Mechanisms Underlying Sparse Coding of Active Touch. Neuron, 69(6), 1160–1175. https://doi.org/10.1016/j.neuron.2011.02.022

[4] Celikel, T., Szostak, V. A., & Feldman, D. E. (2004). Modulation of spike timing by sensory deprivation during induction of cortical map plasticity. Nature Neuroscience, 7(5), 534–541. https://doi.org/10.1038/nn1222

[5] Egger, V., Nevian, T., & Bruno, R. M. (2008). Subcolumnar Dendritic and Axonal Organization of Spiny Stellate and Star Pyramid Neurons within a Barrel in Rat Somatosensory Cortex. Cerebral Cortex, 18(4), 876–889. https://doi.org/10.1093/cercor/bhm126

[6] Banerjee, A., González-Rueda, A., Sampaio-Baptista, C., Paulsen, O., & Rodríguez-Moreno, A. (2014). Distinct mechanisms of spike timing-dependent LTD at vertical and horizontal inputs onto L2/3 pyramidal neurons in mouse barrel cortex. Physiological Reports, 2(3), e00271. https://doi.org/10.1002/phy2.271


