make_new_thalamic_input = 1;                                % 0 or 1: make new spike trains or use existing
make_new_thalamic_kernels = 1;                              % 0 or 1: make new kernels or use existing
make_new_connectivity = 1;                                  % make new (1) or use existing (0) connectivity

includemodulationyn = 0;                                    % direct whisker modulation of L23 neurons (0 = no, 1=yes), see Crochet et al. 2011
includeSTDPyn = 0;                                          % STDP (0 = no, 1=yes)

% How to save everything
savename = ['Test_sim_Svoboda'];

run_sim(make_new_thalamic_input,make_new_thalamic_kernels,make_new_connectivity,includemodulationyn,includeSTDPyn, savename)