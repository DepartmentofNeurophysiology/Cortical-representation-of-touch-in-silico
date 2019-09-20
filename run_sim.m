function run_sim(make_new_thalamic_input,make_new_thalamic_kernels,make_new_connectivity,includemodulationyn,savename)

%% General settings
f = filesep;
addpath(genpath(['.']))

savefolder = ['Simulation_results' f savename '_data' f];

%% Specific settings thalamic input fromSvoboda data
if make_new_thalamic_input
    % What whisker data to use; here specific for Svoboda data
    SvobodaStruct.loadfolder = ['..' f 'Input data' f];         % folder where input data are stored
    SvobodaStruct.animal = 'an171923';                          % animal ID
    SvobodaStruct.sessionvec = {'2012_06_04'};                  % which sessions to load
    SvobodaStruct.dataname = 'data_struct';                     % addition on file names
    SvobodaStruct.volume = 2;                                   % trials for which recorded volume to use (see explanation Svoboda data)
    SvobodaStruct.window.start = 'first touch';                 % How to align trials ('first touch', 'first' or 'pole in reach')
    SvobodaStruct.window.window = [-2000,4000];                   % window (ms) around start time (above)
%     SvobodaStruct.trialvec = [8,8,9];                           % (optional) which of the selected trials to use
    
    savename = [savename '_' SvobodaStruct.animal '_' date];
    savefolder = ['Simulation_results' f savename '_data' f];

    % Thalamic spike trains (filter neurons responding to whiser data above)
    SvobodaStruct.Nkernel_ba = 80;                              % # base angle kernels ('neurons')
    SvobodaStruct.Nkernel_c  = 80;                              % # curvature kernels ('neurons')
    SvobodaStruct.Nkernel_m  = 40;                              % # mixed kernels ('neurons')
end

if ~exist(savefolder, 'dir')
    mkdir(savefolder)
end

%% Connectivity 
if make_new_connectivity
    disp('Generating connectivity')
    % Make barrel-grid
    Nbx = 3;                                                % # of barrels 'x-direction'
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
    barrelstruct{2,1}.mainbarrel    = 1; % main
    barrelstruct{1,1}.mainbarrel    = 2; % secondary
    barrelstruct{3,1}.mainbarrel    = 2; % tertiary
    
    % Generate connectivity
    generate_connectivity(barrelstruct, savefolder, savename)
else
    load([savefolder 'CMDMs_' savename], 'barrelstruct');
    [Nbx, Nby] = size(barrelstruct);
end

% reorganize the connectivity into a single list in order to speed up simulations.
ConMatfilename = ['CMDMs_' savename]; % connectivity matrix file
Nbarrel = Nbx*Nby;
ConData = reorganize_conmat(savename ,savefolder, savefolder, ConMatfilename, includemodulationyn); 

%% Make/load thalamic spike trains 

if make_new_thalamic_input
    % Make new spike trains (from Svoboda recordings)
    SpikeTrainStruct = make_thalamic_spike_trains_svoboda_recordings(savefolder, savename, SvobodaStruct, barrelstruct, make_new_thalamic_kernels);
    TSim = SvobodaStruct.window.window(2)-SvobodaStruct.window.window(1);
else
    % Load existing spike trains
    load([savefolder savename '_Thalamic_Spike_Trains']);
    load([savefolder savename '_Thalamic_Kernels']);
    TSim = length(SpikeTrainStruct{1}.PSTH{1,1})*(KernelStruct{1}.kerneltime(2)-KernelStruct{1}.kerneltime(1));
end
% Check
Input_spike_trains = check_reorganize_spike_trains(SpikeTrainStruct, barrelstruct);
clear SpikeTrainStruct WhiskerTrace barrelstruct KernelStruct

%% Load whisker angles for direct modulation (optional)
% Needs to be made

%% Simulation parameters
simdata.TSim = TSim;                                                                % (ms) Length of the simulation: can be single number (so same length for all trials) or an array of numbers
simdata.timestep = 0.1;                                                             % ms
simdata.Trials = 1:2;                                                               % This is the amount of times the simulation is run for the same initial conditions/parameters/trace etc
simdata.inputvec = [1,2];                                                           % Which input spike traces (trials) from SpikeTrainStruct to use (optional, if empty all will be used).

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
initdata.ExternalInput = zeros(simdata.Nsim, ConData.NAll,simdata.TSim/simdata.timestep); % (optional) random inputs to each cell 

%% Seeding (optional)
% NB should be as long as nr of trials, or non-existent
seeds.initseed = [3,3,3,3,3,3,2,2]; % for seeding initial values (if stdV>0 or stdu0>0)
seeds.runseed  = [5,5,5,5,2,2,2,2];  % for seeding synaptic failures and amplitudes
%% Run simulation
run_trials(ConData, initdata, simdata, Input_spike_trains, seeds)

%% Make 'luminescence' data: convolve and downsample, see Vogelstein et al 2009
disp('Calculating calcium data')
Para.frame_rate_c = 7/1000; % : sampling rate calcium calculations (kHz)
Para.tau_c = 500;           % : calcium decay time constant (ms)
Para.Ca_b = 0.1;            % : baseline calcium concentration (microM)
Para.A_c = 5;               % : calcium concentration 'jump' for each spike
Para.sigma_c = 1;           % : std noise for calcium signal
Para.alpha = 1;             % : scale fluorescence signal
Para.beta =0;               % : offset fluorescence signal
Para.sigma_f =1;            % : std noise for luminescence
Para.tmax = TSim;
simulation_to_calcium( Para, savefolder,savename, initdata, simdata.Nsim);

%% Plot
plot_simulation(simdata.Nsim, initdata, ConData)

