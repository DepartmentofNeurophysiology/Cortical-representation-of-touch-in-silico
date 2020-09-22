function run_sim_compare_data(InputDataStruct,make_new_thalamic_kernels,make_new_connectivity,do_new_simulations, includemodulationyn,includeSTDPyn, CalciumPara, savename)

%% General settings
f = filesep;
addpath(genpath(['.']))

savefolder = [InputDataStruct.savefolder 'Simulation_results' f savename '_data' f];

%% Specific settings thalamic input from Svoboda data
savename_input = [savename '_' InputDataStruct.animal '_vol' num2str(InputDataStruct.volume)];
if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end

%% Connectivity 
if make_new_connectivity
    disp('Generating connectivity')
    % Make barrel-grid
    Nbx = 1;                                                % # of barrels 'x-direction'
    Nby = 1;                                                % # of barrels 'y-direction'
    barrelstruct = cell(Nbx, Nby);
    for nbx = 1:Nbx
        for nby = 1:Nby
            barrelstruct{nbx, nby}.xpos         = (nbx-2)*300;      % barrel position
            barrelstruct{nbx, nby}.ypos         = (nby-1)*300;
            barrelstruct{nbx, nby}.Nthalamic    = 200;              % # filter neurons for this barrel
            barrelstruct{nbx, nby}.mainbarrel   = 3;
        end
    end
    % Choose main and secondary barrels
    barrelstruct{1,1}.mainbarrel    = 1; % main
    
    % Generate connectivity
    generate_connectivity(barrelstruct, savefolder, savename)
else
    load([savefolder 'CMDMs_' savename], 'barrelstruct');
end

% reorganize the connectivity into a single list in order to speed up simulations.
ConMatfilename = ['CMDMs_' savename]; % connectivity matrix file

if exist([savefolder ConMatfilename '_ConData.mat'],'file') == 2
    load([savefolder ConMatfilename '_ConData.mat'])
else 
    ConData = reorganize_conmat(savename ,savefolder, savefolder, ConMatfilename, includemodulationyn);    
end

%% Make/load thalamic spike trains 

if InputDataStruct.make_new_thalamic_input
    % Make new spike trains (from Svoboda recordings)
    SpikeTrainStruct = make_thalamic_spike_trains_svoboda_recordings(savefolder, savename_input, InputDataStruct, barrelstruct, make_new_thalamic_kernels);
    TSim = InputDataStruct.window.window(2)-InputDataStruct.window.window(1);
else
    % Load existing spike trains
    load([savefolder savename_input '_Thalamic_Spike_Trains']);
    load([savefolder savename '_Thalamic_Kernels']);
    TSim = length(SpikeTrainStruct{1}.PSTH{1,1})*(KernelStruct{1}.kerneltime(2)-KernelStruct{1}.kerneltime(1));
end
% Check
Input_spike_trains = check_reorganize_spike_trains(SpikeTrainStruct, barrelstruct);
clear SpikeTrainStruct KernelStruct SvobodaStruct

%% Simulation parameters
simdata.TSim = TSim;                    % (ms) Length of the simulation: can be single number (so same length for all trials) or an array of numbers
simdata.timestep = 0.1;                 % ms
simdata.Trials = 1;                   % This is the amount of times the simulation is run for the same initial conditions/parameters/trace etc
simdata.inputvec = [];               % Which input spike traces (trials) from SpikeTrainStruct to use (optional, if empty all will be used).
if includeSTDPyn
    simdata.STDP = true;                    % Turn on STDP          
else
    simdata.STDP = false;                    % Turn off STDP        
end

% Parameters and variables to save for each trial:
simdata.whattosave = {'Neuron_Para', 'Tau_plas', 'cellinfo_all', 'cellinfo_input', ...
        'modelsc', 'modelspt', 'inputsc', 'inputspikes', 'simLen', 'step', ...
        'V0','U0','VT0','vr', 'V', 'seeds', 'timestepseed_input', ...
        'timestepseed_model','simdata', 'final_connectivity', 'initial_connectivity'};
% Now left out (too large), but can be added: 
%   * initdata 
%   * MIn (whisker modulation for each cell, can easily be reconstructed)
%   * U (recovery variable U for each cell) 
%   * VT (dynamic threshold trace for each cell)    

%% Make whisker angles for direct modulation (optional)  
if includemodulationyn
    simdata.whiskermodulation = true;   % Turn on direct whisker modulation
    if make_new_thalamic_input
        load([savefolder savename_input '_Thalamic_Spike_Trains'], 'WhiskerTrace')
        WhPara = WhiskerPara_direct_modulation(WhiskerTrace, barrelstruct, simdata.timestep, savefolder, savename_input);
    else
        load([savefolder savename_input '_WhiskerModulation']);
    end
else
    simdata.whiskermodulation = false;   % Turn off direct whisker modulation
    WhPara = [];
end
clear WhiskerTrace barrelstruct 

%% Initial conditions
% Initial parameters
initdata.Vrest      = [-65];            % mean resting membrane potential simulated cells; can be list (multiple simulations)
initdata.stdVrest   = [0];              % should have same size as Vrest   
initdata.V0         = -65;              % mean starting membrane potential simulated cells; can be list (multiple simulations)
initdata.stdV0      = 5;                % should have same size as V0
initdata.u0         = 0;                % mean starting u simulated cells; can be list (multiple simulations)
initdata.stdu0      = 0;                % should have same size as u0

% Thresholds
initdata.Vthresdyn  = 1;                % whether to use a dynamic threshold (1) or not (0) for excitatory neurons 
initdata.setVthres.type = 'pertype';    % how to set the threshold: 'distribution', 'pertype' or 'individual'
if strcmp(initdata.setVthres.type, 'distribution')
    initdata.Vthres     = -55;          % mean threshold potential simulated cells; ; can be list (multiple simulations)
    initdata.stdVthres  = 0;            % should have same size as Vthres 
    initdata.setVthres.nsim = length(initdata.Vthres);
elseif strcmp(initdata.setVthres.type, 'pertype')
    % NB mean and std per celltype are set in ConData.Neuron_Vt; this is used
    initdata.setVthres.nsim = 1;        % multiple realizations can be done with Trials
    initdata.Vthresvar  = 0;            % whether to use an initial distribution around the mean (1) or not (0)
elseif strcmp(initdata.setVthres.type, 'individual')
    % set inidividual values, like shown in next section
end

% total # simulations
if isempty(simdata.inputvec)    
    inputvec = 1:length(Input_spike_trains);
else
    inputvec = simdata.inputvec;
end
simdata.Nsim =length(inputvec)*length(initdata.Vrest)*initdata.setVthres.nsim*length(initdata.V0)*length(initdata.u0)*length(simdata.Trials); 
% NB simulations are done in the order of for loops above (so Trials is innermost loop) 

%% Set initial conditions for individual cells (optional)
% If you want to set parameters / initial conditions for cells specific, set this in
% setindcell, in a matrix of (Number of neurons NAll x total # simulations) 
% Here for example: Vrest L4 is more depolarized
% Note: this overwrites anything that would be set with the mean and std of initdata (previous section)!

% NsimperVr = length(initdata.setVthres.nsim)*length(initdata.V0)*length(initdata.u0)*length(simdata.Trials);
% initdata.setindcell.Vrest = zeros(ConData.NAll, Nsim);
% % set all Vrest to initdata.Vrest
% for tt = 1:length(initdata.Vrest)
%     initdata.setindcell.Vrest(:, (tt-1)*NsimperVr+1:tt*NsimperVr) = initdata.Vrest(tt)*ones(ConData.NAll,NsimperVr);
% end
% % make L4 cells more depolarized
% initdata.setindcell.Vrest(ConData.Cellinfo_All(:,4)<3,:) = initdata.setindcell.Vrest(ConData.Cellinfo_All(:,4)<3,:)+4.1*ones(sum(ConData.Cellinfo_All(:,4)<3),Nsim);
% initdata.setindcell.V0 = initdata.setindcell.Vrest;

%% External inputs (for instance random, optional)
% initdata.ExternalInput = zeros(simdata.Nsim, ConData.NAll,simdata.TSim/simdata.timestep); % (optional) random inputs to each cell 

%% Seeding (optional)
% NB should be as long as nr of trials, or non-existent
% seeds.initseed = [3,3,3,3,3,3,2,2];  % for seeding initial values (if stdV>0 or stdu0>0)
% seeds.runseed  = [5,5,5,5,2,2,2,2];  % for seeding synaptic failures and amplitudes
%% Run simulation
savename_sims = savename_input;
if do_new_simulations
    ConData.FnametoSave = savename_sims;
    run_trials(ConData, initdata, simdata, Input_spike_trains, [], WhPara)
end

%% Make 'luminescence' data: convolve and downsample, see Vogelstein et al 2009
if CalciumPara.make_new_calcium_data
    disp('Calculating calcium data')
    simulation_to_calcium( CalciumPara, savefolder,savename_sims, initdata, simdata.Nsim);
end

%% Plot
% plot_simulation(simdata.Nsim, initdata, ConData, savename_input)

