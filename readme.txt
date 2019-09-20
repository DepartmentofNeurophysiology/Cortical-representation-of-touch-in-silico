———————————————————————— INSTALLATION————————————————————
Download all the files and folders, preserving the folder structure. The model will run when this structure is added to the path (addpath(genpath([‘.’])), as is done in run_sim.m).

————————————————————————— GENERAL————————————————————————
A simple example simulation can be run by executing the script quick_example, which uses the function ‘run_sim’ (this function controls the simulation, and can be used as an example for doing user-defined simulations). A 3x1 grid of barrels in L4 and L23 will be constructed, with corresponding thalamic barreloids. Thalamic input spike trains will be generated in reaction whisker input (in this case: whisker base angles and curvatures from the lab Karel Svoboda, publicly available on Data Sets — CRCNS.org), and the cortical model will respond to these thalamic spike trains. 

—————————————————————————— INPUT—————————————————————————
As input, data files (whisker base angle and curvature) of the lab Karel Svoboda were used (publicly available on Data Sets — CRCNS.org). The specified whisker data files are expected to be present in folder ‘Input data’. The function ‘make_thalamic_spike_trains_svoboda_recordings’ in folder ‘Make_New_Thinput’ makes spike trains from the whisker data. For the use of any other input, make a file similar to this example. 

————————————————————————— SIMULATION—————————————————————
A full simulation includes:
	1.	generating new input spike trains (make_new_thalamic_input = 1), relevant functions in ‘Make_New_Thinput’
	2.	making new thalamic filter neurons (make_new_thalamic_kernels = 1), relevant functions in ‘Make_New_Thinput’
	3.	Make a new realisation of the network connectivity (make_new_connectivity=1), relevant functions in ‘Make_New_Connectivity’
	4.	Initialising and running the network, relevant functions in ‘Network Simulations’

————————————————————————— RESULTS————————————————————————
This will result in the following files in folder ‘Simulation Results -> (user defined name)
	• cellinfo: information about all cells: (Number of Cells-by 6) matrix, with
		• the first 3 columns are the location of each cell 
			• first column: position from the middle of the barrel (~rostral-caudal)
			• second column: position from the middle of the barrel (~dorsal-ventral) 
			• third column: depth, distance from pia
		• 4th column is the cell type 
			• 1 - L4 spiny stallet
			• 2 - L4 pyramidal
			• 3 - L4 Fast spike 
			• 4 - L4 low-threshold spike
			• 5 - L2/3 pyramidal 
			• 6 - L2/3 PV+ fast spike
			• 7 - L2/3 PV+ chandler
			• 8 - L2/3 PV+ bursting 
			• 9 - L2/3 SOM+ martinotti 
			• 10 - L2/3 SOM+ bitufted 
			• 11 - VIP+ double bouquet 
			• 12 - VIP+ bipolar
			• 13 - CR+ bipolar
			• 14 - CR+ multipolar/basket			
			• 15 - neurogliaform 
			Note: some L2/3 types have been merged with other types so in the matrix the numbers are 0); 
		• 5th column is the barrel identification (different number indication different barrels). 
		• 6th column is Exe/Inh indicator (1 for excitatory and -1 for inhibitory neurons)
	• CMDMs: all connectivity data
	• CMDMs_(…)_ConData
	• CMDMs_(…)_ParaMat
	• CMDMs_(…)_ParaMat_reduced: a reduced version of the previous set, to speed up simulations
	• (…)_initialsettings: the initial settings to each simulation
	• (…)_Simcolumn_(…)_(simulation number): spike times and membrane potentials of each simuation
	• (…)_Thalamic_Kernels: thalamic filters for each thalamic neuron
	• (…)_Thalamic_Spike_Trains: all input spike trains and whisker traces