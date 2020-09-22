close all
clear all

f = filesep;
addpath(genpath(['..']))

edges = [-5:0.1:5];
edgesav = -1:0.02:1;

%% Svoboda data params
InputDataStruct.loadfolder = '/Users/fleurzeldenrust/Documents/Data/Svoboda data/data/';
InputDataStruct.savefolder = '/Users/fleurzeldenrust/Documents/Data/Svoboda data/figures/compare for paper/';
InputDataStruct.animal = 'an171923';
InputDataStruct.session = {'2012_06_04'};
InputDataStruct.dataname = 'data_struct';

vol_l23 = 6;
vol_l4 = 8;

depth_l23 = [150 250]; 
depth_l4 = [250 400];

InputDataStruct.window.start = 'first touch';
InputDataStruct.window.window = [-2000,8000];

make_plot_Svoboda_data = 1;

%% Network params
% Input           
InputDataStruct.trialvec = [];                                % include all trials     
% Thalamic spike trains (filter neurons responding to whiser data above)
InputDataStruct.Nkernel_ba = 80;                              % # base angle kernels ('neurons')
InputDataStruct.Nkernel_c  = 80;                              % # curvature kernels ('neurons')
InputDataStruct.Nkernel_m  = 40;                              % # mixed kernels ('neurons')

includemodulationyn = 0;
includeSTDPyn = 0;
savename = 'Compare_data_sims';
InputDataStruct.savename = savename;

% for calcium traces
% NB The model is overparametrized, so the only data to fit are: mean and
% std of recordings, frame rate and decay time. The rest is arbitrary
mu_exp = 0.015;
sigma_exp = 0.15;
CalciumPara.frame_rate_c = 7/1000;  % : sampling rate calcium calculations (kHz)
CalciumPara.tau_c = 500;            % : calcium decay time constant (ms)
CalciumPara.Ca_b = 0;               % : baseline calcium concentration (microM) 
CalciumPara.A_c = 5;                % : calcium concentration 'jump' for each spike according to Vogelstein paper, baseline and jumps are relatively similar in size. 
CalciumPara.alpha = 1;              % : scale fluorescence signal
CalciumPara.beta =mu_exp;           % : offset fluorescence signal
CalciumPara.sigma_f = sigma_exp;    % : std noise for added for fluorescence signal
CalciumPara.sigma_c = 0;            % : std noise for calcium signal
CalciumPara.tmax = InputDataStruct.window.window(2)-InputDataStruct.window.window(1);

savefolder_sims = [InputDataStruct.savefolder 'Simulation_results' f savename '_data' f];

make_new_connectivity=0;
make_new_thalamic_kernels=0;
InputDataStruct.make_new_thalamic_input = 0;
do_new_simulations = 0;
CalciumPara.make_new_calcium_data = 0;
make_new_plots_sims = 0;

InputDataStruct_l23 = InputDataStruct;
InputDataStruct_l23.volume = vol_l23;                                   % trials for which recorded volume to use (see explanation Svoboda data)
InputDataStruct_l4 = InputDataStruct;
InputDataStruct_l4.volume = vol_l4;                                   % trials for which recorded volume to use (see explanation Svoboda data)


%% Load Svoboda data
if make_plot_Svoboda_data 
    [f_example_1vol_l23, f_average_1vol_l23,f_example_1vol_l4, f_average_1vol_l4, zrange_vol_l23, zrange_vol_l4] = load_plot_Svoboda_data_to_compare(InputDataStruct, vol_l23, vol_l4, depth_l23, depth_l4, edges, edgesav);
    InputDataStruct_l23.depth = zrange_vol_l23;
    InputDataStruct_l4.depth = zrange_vol_l4;
end

%% Simulate corresponding trials

% L23
disp('Simulating L23')
run_sim_compare_data(InputDataStruct_l23,make_new_thalamic_kernels,make_new_connectivity,do_new_simulations, includemodulationyn,includeSTDPyn, CalciumPara, savename)

% L4 
disp('Simulating L4')
run_sim_compare_data(InputDataStruct_l4,make_new_thalamic_kernels,make_new_connectivity,do_new_simulations, includemodulationyn,includeSTDPyn, CalciumPara, savename)

close all
%% Plot simulations
if make_new_plots_sims
    % L23
    disp('Plotting L23')
    savename_sims = [savename '_' InputDataStruct_l23.animal '_vol' num2str(vol_l23)];
    [f_trial1_l23_sims, f_average_l23_sims, f_trial1_l4_sims, f_average_l4_sims] = plot_activity_simulations(CalciumPara, InputDataStruct_l23, savefolder_sims, savename, savename_sims, 'mouse',edges, edgesav,  0);
    savefig(f_trial1_l23_sims,[InputDataStruct.savefolder 'single_trial_l23_vol' num2str(vol_l23) '_sims.fig'],'compact')
    savefig(f_average_l23_sims,[InputDataStruct.savefolder 'average_l23_vol' num2str(vol_l23) '_sims.fig'],'compact')

    % L4
    disp('Simulating L4')
    savename_sims = [savename '_' InputDataStruct_l4.animal '_vol' num2str(vol_l4)];
    [f_trial1_l23_sims, f_average_l23_sims, f_trial1_l4_sims, f_average_l4_sims] = plot_activity_simulations(CalciumPara, InputDataStruct_l4, savefolder_sims, savename, savename_sims, 'mouse',edges, edgesav,  0);
    savefig(f_trial1_l4_sims,[InputDataStruct.savefolder 'single_trial_l4_vol' num2str(vol_l4) '_sims.fig'],'compact')
    savefig(f_average_l4_sims,[InputDataStruct.savefolder 'average_l4_vol' num2str(vol_l4) '_sims.fig'],'compact')
end


%% Make combined figures
close all

fig_l23 = combine_figs(InputDataStruct, vol_l23, 'l23');
savefig(fig_l23,[InputDataStruct.savefolder fig_l23.Name],'compact')

fig_l4 = combine_figs(InputDataStruct, vol_l4, 'l4');
savefig(fig_l4,[InputDataStruct.savefolder fig_l4.Name],'compact')




