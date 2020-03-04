function ConData = simcolumn_indtrial(ConData, initdata, simdata, inputspikes, seeds, WhPara, varargin)
% Simcolumn_indtrial 
% * runs single trial simulation (needed: 
% * sets for this trial initial Vr, Vt dynamic threshold, V and u
% * saves each trial into a mat file

% INPUTS:
% * Data: structure with connectivity data and parameter values from reorganize_conmat
% * initdata: structure with initial conditions: Vrest(mean resting
% membrane potential; can be list); stdVrest; Vthres (mean threshold
% potential; stdVthres; V0 (mean starting membrane potential); stdV0; u0;
% (mean starting u); stdu0; Vthresdyn (dynamic threshold (1) or not (0));
% setVthres.type (how to set the threshold: 'distribution', 'pertype' or
% 'individual') and optional setindcell (.Vrest; .Vthres, .V0, .u0) if one wants to set
% them for each cell individually. In this case, set the desired value in
% setindcell, in a matrix of (Number of neurons NAll x total # simulations) ,  and (optional) ExternalInput = NAll x SimLen
% array with random inputs to each cell 
% * simdata: structure with simulation data: Tsim (total length simulation)
% and timestep, both in ms
% * inputspikes: cell structure with input spikes inputspikes{# iterations per condition}(Ncell x Nspikes_max) 
% * seeds (optional): structure with seeds to re-do simulations: initseed = integer (seed for initial conditions); runseed = integer (seed for running: synaptic failures and
% amplitude)
% * WhPara (optional): structure with base-angles and amplitudes whiskers for whisker
% modulation l2/3 (Crochet 2011)

% GPU version of simulation
% the strategy is to put everything on GPU and run from there
%run the simulation with full matrix operations

% seeds are optional: seeds.initseed for initial values membrane potential
% and threshold, seeds.runseed for running the simulations (synaptic
% failures and random amplitude synapse)




%% neuronal model parameters
a=ConData.Neuron_Para(:,1);
b=ConData.Neuron_Para(:,2);
% b(Data.Cellinfo_All(:,4)>5) = b(Data.Cellinfo_All(:,4)>5) + 0.3;
b(ConData.Cellinfo_All(:,4)==5) = b(ConData.Cellinfo_All(:,4)==5) - 0.15;
c=ConData.Neuron_Para(:,3);
d=ConData.Neuron_Para(:,4);
d(ConData.Cellinfo_All(:,4)>5) = d(ConData.Cellinfo_All(:,4)>5) + 2;

Tau_plas = 120; % time constant for short term synaptic dynamics in ms; 

% parameters for dynamic threshold model
DynThresMat = nan*ones(ConData.NAll, 6);
if initdata.Vthresdyn == 1
    for nt = 1:length(ConData.NtAll)
        cells = find(ConData.Cellinfo_All(:,4)==nt);
        ncells = length(cells);
        DynThresMat(cells,1) = ConData.Neuron_Vt.alpha(nt) *ones(ncells,1);
        DynThresMat(cells,2) = ConData.Neuron_Vt.Vi(nt)    *ones(ncells,1);
        DynThresMat(cells,3) = ConData.Neuron_Vt.Vmin(nt)  *ones(ncells,1);
        DynThresMat(cells,4) = ConData.Neuron_Vt.Ka(nt)    *ones(ncells,1);
        DynThresMat(cells,5) = ConData.Neuron_Vt.Ki(nt)    *ones(ncells,1);
        DynThresMat(cells,6) = ConData.Neuron_Vt.tau(nt)   *ones(ncells,1);
    end
end

%% simulation time
step = simdata.timestep;
simLen = simdata.TSim;
Nsim = simLen/step;

%% Crochet Whisker modulation
% the membrane potential of L2/3 neurons is directly modulated by whisking
% (Crochet et al, Neuron 2011)
% * model: Data.WhiskModel made in 'reorganize_conmat'
% * data: structure WhPara with base-angles (degrees), amplitudes and phase of the whiskers
% Assumptions:
% * Whisker angle correlated inputs to L2/3 population; 
% * not all neurons are modulated by whisking 
% * the distribution of whisking modulation is depth dependent (90% in L3 and 20% in L2); 
% * whisker modulated Neurons have evenly distributed phase preference; 
% * modulation strength is normally distributed with 0.08 +- 0.02 mV/deg; 
% * additional white noise inputs to reduce the r2 (in petersen paper the r2 is low; not
% implemented now)
% * only neurons in the principal barrel are modulated
MIn = zeros(ConData.NAll, Nsim);
if simdata.whiskermodulation
    if (exist('WhPara', 'var') && (~isempty(WhPara) && isfield(ConData, 'WhiskModel')))
        disp('Calculating parameters whisker modulation')
        % needed: whisker base angle trace for this trial in degrees,
        % Amplitude and Phase, it is calculated with
        % WhiskerPara_direct_modulation
        WhPara.MeanAm = nanmean(WhPara.Baseangles_degrees);
        WhPara.MinAm = nanmin(WhPara.Baseangles_degrees);
        WhPara.MeanProtraAm = nanmean(WhPara.Amplitude);
        if ~(length(WhPara.Phase) == Nsim)
            error('Please use a whisker phase-trace that has the same length as the curvature trace')
        else
            MIn = 0.5*WhPara.MeanProtraAm*repmat(ConData.WhiskModel.k,[1,Nsim]).*(sin(repmat(WhPara.Phase,[ConData.NAll,1]) - repmat(ConData.WhiskModel.Pha0, [1,Nsim])));
        end
        if ~(size(MIn,1) == ConData.NAll)
            error('Size matrix direct whisker modulation does not have correct number of neurons')
        end
        if ~(size(MIn,2) == Nsim)
            error('Size matrix direct whisker modulation does not have correct number of time steps')
        end
    elseif ~exist('WhPara', 'var') &&  isfield(ConData, 'WhiskModel')
        disp('Whisker modulation model (Crochet 2011) defined, but no data given; ignoring model.')

    elseif  exist('WhPara', 'var') && ~isfield(ConData, 'WhiskModel')
        disp('Whisker modulation data given, but no model (Crochet 2011) defined; ignoring data.')
    else
        disp('No direct whisker modulation')
    end
end

%% Random inputs to network
if isfield(initdata, 'ExternalInput')
    disp('Using given external input')
    RIn = initdata.ExternalInput;
else
    disp('No external input given')
    RIn = zeros(ConData.NAll,Nsim);
end

%% assign initial parameters for vr, vt, v, u
if (exist('seeds','var') && (~isempty(seeds) && isfield(seeds, 'initseed')))
    rng(seeds.initseed)
else
    rng('shuffle')
    scurr = rng;
    seeds.initseed = scurr.Seed;
end


vr=initdata.Vrest*ones(ConData.NAll,1) + initdata.stdVrest*randn(ConData.NAll,1);     % resting membrane potential

% threshold for spiking
if strcmp(initdata.setVthres.type, 'distribution')
    vt=initdata.Vthres*ones(ConData.NAll,1) + initdata.stdVthres*randn(ConData.NAll,1);   
else
    if strcmp(initdata.setVthres.type, 'pertype')
        % NB Neurons are in order of type, so one can loop over types
        Nt_temp = [0, cumsum(ConData.NtAll)];
        for i = 1:length(Nt_temp) - 1
            vt(Nt_temp(i)+1:Nt_temp(i+1),1) = ConData.Neuron_Vt.avg(i) + initdata.Vthresvar*ConData.Neuron_Vt.std(i)*randn(ConData.NtAll(i), 1);
            % The next parameters are only for the dynamic threshold model (fixed
            % threshold does not use them)
            % currently spike threshold of excitatory neurons are dynamic, while those
            % of inhibitory neurons are static
            VtModel.sl(Nt_temp(i)+1:Nt_temp(i+1), 1) = ConData.Neuron_Vt.sl(i);
            VtModel.a(Nt_temp(i)+1:Nt_temp(i+1), 1) = ConData.Neuron_Vt.a(i);
        end
    elseif strcmp(initdata.setVthres.type, 'individual')
        vt = initdata.setindcell.Vthres;
    else
        error('Choose method to set threshold: initdata.setVthres.type should be distribution, pertype or individual')
    end
end
VT0 = vt;
    
% initial values for neuron model
v=initdata.V0.*ones(ConData.NAll,1)+initdata.stdV0.*randn(ConData.NAll,1); % Initial values of v
u=initdata.u0.*ones(ConData.NAll,1)+initdata.stdu0.*randn(ConData.NAll,1); % Initial values of u
% required initial value for computation

if isfield(initdata, 'setindcell')
    namevec = fieldnames(initdata.setindcell);
    for ii = 1:length(namevec)
        if strcmp(namevec{ii}, 'Vrest')
            vr = initdata.setindcell.Vrest;
        elseif strcmp(namevec{ii}, 'Vthres')
            % do nothing; already done
        elseif strcmp(namevec{ii}, 'V0')
            v = initdata.setindcell.V0;
        elseif strcmp(namevec{ii}, 'u0')
            u = initdata.setindcell.u0;
        else
            error('Non-existing field for setindcell')
        end
    end
end

V0 = v;
U0 = u;

%% Preallocate
inputsc = single(zeros(ConData.NIn,1));
modelsc = single(zeros(ConData.NAll,1));
inputin = single(zeros(ConData.NAll,1));
modelin = single(zeros(ConData.NAll,1));
Finput = [];
Fmodel = [];


Vminus = v;
V = nan*ones(ConData.NAll, Nsim);
U = nan*ones(ConData.NAll, Nsim);
if initdata.Vthresdyn == 1
    VT = nan*ones(ConData.NAll, Nsim);
end
dS = zeros(ConData.NAll, 1);
vt_sp=zeros(ConData.NAll,1);
vtct=zeros(ConData.NAll,3);
modelvt=zeros(ConData.NAll, 1);
modelspt=zeros(ConData.NAll, 1);

inputspikes(size(inputspikes,1)+1:ConData.NIn,:) = 0;
inputspikes(ConData.NIn+1:end, :) = [];

%% matrix for running variables
% these varibles run in GPU
STD_input = (ConData.PMat_IntoAll.STD);
STD_model = (ConData.PMat_AlltoAll.STD);

g_input = (ConData.PMat_IntoAll.g);
s_input = (ConData.PMat_IntoAll.s);
g_model = (ConData.PMat_AlltoAll.g);
s_model = (ConData.PMat_AlltoAll.s);

CN_input = (ConData.PMat_IntoAll.CN);
CN_model = (ConData.PMat_AlltoAll.CN);

Trise_input = (ConData.PMat_IntoAll.Trise);
Tfall_input = (ConData.PMat_IntoAll.Tfall);
Trise_model = (ConData.PMat_AlltoAll.Trise);
Tfall_model = (ConData.PMat_AlltoAll.Tfall);

inputspt_temp = inputspikes;

%% variables needed to run STDP
spikepairs_IntoAll = zeros(length(ConData.PMat_IntoAll.Am), 4);
spikepairs_AlltoAll = zeros(length(ConData.PMat_AlltoAll.Am), 4);
Ampost_IntoAll = ConData.PMat_IntoAll.Am;
Ampost_AlltoAll = ConData.PMat_AlltoAll.Am;
initial_connectivity = Ampost_AlltoAll;
% Ampre_IntoAll = Data.PMat_IntoAll.Am;
% Ampre_AlltoAll = Data.PMat_AlltoAll.Am;
% celltype is organized as postType, postEI, preType, preEI
Celltype_IntoAll = [ConData.Cellinfo_All(ConData.PMat_IntoAll.preCell(:,1), 4), ConData.Cellinfo_All(ConData.PMat_IntoAll.preCell(:,1), 6), ...
    ConData.Cellinfo_In(ConData.PMat_IntoAll.preCell(:,2), 4), ConData.Cellinfo_In(ConData.PMat_IntoAll.preCell(:,2), 6)];
Celltype_AlltoAll = [ConData.Cellinfo_All(ConData.PMat_AlltoAll.preCell(:,1), 4), ConData.Cellinfo_All(ConData.PMat_AlltoAll.preCell(:,1), 6), ...
    ConData.Cellinfo_All(ConData.PMat_AlltoAll.preCell(:,2), 4), ConData.Cellinfo_All(ConData.PMat_AlltoAll.preCell(:,2), 6)];

%% simulation starts
if exist('seeds','var') && isfield(seeds, 'runseed')
    rng(seeds.runseed)
else
    rng('shuffle')
    scurr = rng;
    seeds.runseed = scurr.Seed;
end
% seeds for synaptic failures and Am
timestepseed_input = randi(2^32, 1, Nsim);
timestepseed_model = randi(2^32, 1, Nsim);

for t=1:Nsim % simulation of Ni ms

    tm = t*step;

    

        


    [ind_input,ind_sc]=find(inputspt_temp <= tm & inputspt_temp > 0);
    inputspt_temp(sub2ind(size(inputspikes), ind_input, ind_sc)) = 0;
    
    inputpre = [];
    
    if isempty(ind_input) == 0
        for i = 1:length(ind_input)
            inputpre = [inputpre; ConData.PMat_IntoAll.postIdx{ind_input(i)}];
        end
        inputsc(ind_input) = inputsc(ind_input)+1;
        clear ind_input ind_sc
        
        %now we can directly get all the parameters (for now just F4)
        %we need Am, Trise, Tfall, Plas and synaptic delay
        %Am will be generated from CV and Am matrix
        
        %  check for input spikes from input layer 
        [Finput, preSpikes_input, timestepseed_input(t)] = simcolumn_synapse_onestep(inputpre, Ampost_IntoAll, STD_input, ConData.PMat_IntoAll.CV, ...
            ConData.PMat_IntoAll.Fail, ConData.PMat_IntoAll.Delay, ConData.PMat_IntoAll.preID, ConData.NAll, Finput, tm, timestepseed_input(t));

        
    end
    
    % update running matrix for the synaptic current model
    if isempty(Finput) == 0
        [s_input, STD_input, Finput] = simcolumn_updatesynap_ver2(s_input, ...
            STD_input, ConData.PMat_IntoAll.Plas, Finput, tm);
    end
    
    idx = find(v>vt);%neuron reaches spike threshold
    idx_spingen = idx(vtct(idx, 3) == modelsc(idx));
    vt_sp(idx_spingen) = vt(idx_spingen);
    vtct(idx_spingen, 1) = 1; % keep track of neurons with spike generation in progress
    vtct(idx_spingen, 2) = tm; % time of first spike threshold cross
    vtct(idx_spingen, 3) = vtct(idx_spingen, 3) + 1;
    
    modelfired = find(v>0);    % find cells those are firing in layer2/3
    
    % check failed spikes
    idx = tm - vtct(:, 2) > 15 & vtct(:, 1) == 1;
    vtct(idx, 1) = 0;
    vtct(idx, 3) = vtct(idx, 3) - 1;
    
    modelpre= [];
    modelpost = [];
    
    if isempty(modelfired) == 0
        for i = 1:length(modelfired)
            % modelpre uses pre-post direction, need to use postIdx (which
            % contains all postsynaptic partner to one presynaptic neuron)
            % modelpost uses post-pre direction to calculate STDP, so use
            %  preIdx instead 
            modelpre  = [modelpre; ConData.PMat_AlltoAll.postIdx{modelfired(i)}];
            modelpost = [modelpost; ConData.PMat_AlltoAll.preIdx{modelfired(i)}];
            modelsc(modelfired(i)) = modelsc(modelfired(i)) + 1;
            modelspt(modelfired(i), modelsc(modelfired(i))) = tm - step;
            modelvt(modelfired(i), modelsc(modelfired(i))) = vt_sp(modelfired(i));
        end
        % reset spike generation tracker
        vtct(modelfired, 1) = 0;
        
        % check for input spikes from simulated layer
        [Fmodel, preSpikes_model, timestepseed_model(t)] = simcolumn_synapse_onestep(modelpre, Ampost_AlltoAll, STD_model, ...
        ConData.PMat_AlltoAll.CV, ConData.PMat_AlltoAll.Fail, ConData.PMat_AlltoAll.Delay, ...
        ConData.PMat_AlltoAll.preID, ConData.NAll, Fmodel, tm, timestepseed_model(t));
        
                

        % update spike pairs, which is used to run STDP 
        spikepairs_AlltoAll(preSpikes_model(:,1), 1) = preSpikes_model(:,2);
        spikepairs_AlltoAll(preSpikes_model(:,1), 3) = 1;
        % post synaptic spike timing
        spikepairs_AlltoAll(modelpost, 2) = tm;
        spikepairs_AlltoAll(modelpost, 4) = 1;
    end
    
    %% STDP rule
    % right now only run it for L4-L23 network; i.e. assume no plasticity
    % in thalamo-cortical synapses
    if simdata.STDP
        spikeidx_AlltoAll = [modelpre; modelpost];
        if ~isempty(spikeidx_AlltoAll)
            [spikepairs_AlltoAll, Ampost_AlltoAll] = simcolumn_Plasticity_STDP...
                (spikepairs_AlltoAll, spikeidx_AlltoAll, Celltype_AlltoAll, Ampost_AlltoAll, 3);
        end
    end
    
    
    %% update running matrix for the synaptic current model
    if isempty(Fmodel) == 0
        [s_model, STD_model, Fmodel] = simcolumn_updatesynap_ver2(s_model, ...
            STD_model, ConData.PMat_AlltoAll.Plas, Fmodel, tm);
    end
    
    
    
    % reset the neuron model after firing
    inter=find(v>0);
    v(inter)=c(inter);
    u(inter)=u(inter)+d(inter);

    % running synaptic current model with Euler method
    % the differential function is
    %dg/dt = ((tau2 / tau1) ** (tau1 / (tau2 - tau1))*s-g)/tau1 
    %ds/dt = -s/tau2 
    g_input = g_input + step*((CN_input.*s_input-g_input)./Tfall_input);
    s_input = s_input + step*( -s_input ./ Trise_input);
    
    g_model = g_model + step*((CN_model.*s_model-g_model)./Tfall_model);
    s_model = s_model + step*( -s_model ./ Trise_model);
    
    % running STD model with Euler method
    % the differential function is
    % dW/dt = -(W-1)/tau_plas
    STD_input = STD_input + step * (-(STD_input-1)./Tau_plas);
    STD_model = STD_model + step * (-(STD_model-1)./Tau_plas);
    
    % calculate synaptic current
    inputin = simcolumn_calcualteI_Mat_ver2(g_input, ConData.PMat_IntoAll.preEI, v);
    modelin = simcolumn_calcualteI_Mat_ver2(g_model, ConData.PMat_AlltoAll.preEI, v);
    

    I = inputin+modelin;
    
    % running neuron model with Euler method
    v_temp = v;
    v=v+step*(0.04*(v-vr).*(v-vt)-u + I + RIn(:,t) + MIn(:,t));
    u=u+step*a.*(b.*(v_temp-vr)-u); 
    
    V(:,t) = v;
    U(:,t) = u;
    
    %% spike threshold model from Fontaine et al Plos Comp Bio 2014
    if initdata.Vthresdyn == 1    
        Vt_thetainf = DynThresMat(:,1).*(v - DynThresMat(:,2))+DynThresMat(:,3)+DynThresMat(:,4).*log(1+exp((v - DynThresMat(:,2))./DynThresMat(:,5)));
        vt = vt + (Vt_thetainf - vt).*step./DynThresMat(:,6);
        
        VT(:,t) = vt;
    end
    % fix spike threshold for neurons under spike generation: 'freeze' the vt value once v crosses the vt. Otherwise the vt keeps changing and sometimes results in false negative spikes 
    vt(vtct(:, 1) == 1) = vt_sp(vtct(:,1) == 1);
    
end


%% Save
final_connectivity = Ampost_AlltoAll;

V = single(V);
U = single(U);
if initdata.Vthresdyn == 1
    VT = single(VT);
end




cellinfo_all = ConData.Cellinfo_All;
cellinfo_input = ConData.Cellinfo_In;
Neuron_Para = ConData.Neuron_Para;

a = exist(ConData.savefolder, 'dir');
if ~(a==7) 
    mkdir(ConData.savefolder)
end

if isfield(simdata, 'whattosave')
    whattosave = simdata.whattosave;
else
    whattosave = {'Neuron_Para', 'Tau_plas', 'cellinfo_all', 'cellinfo_input', ...
        'modelsc', 'modelspt', 'inputsc', 'inputspikes', 'simLen', 'step', ...
        'V0','U0','VT0', 'vr', 'seeds', 'timestepseed_input', ...
        'timestepseed_model', 'initdata','simdata', 'final_connectivity', 'initial_connectivity'};
end

if initdata.Vthresdyn == 1
    thresholdname = ['dynthreshold_set'];
else
    thresholdname = ['fixthreshold_set'];
end



savename = [ConData.savefolder,ConData.FnametoSave '_Simcolumn_' thresholdname '_' initdata.setVthres.type '_simulation_' num2str(initdata.nsim)];
save(savename, whattosave{:}, '-v7.3')

%% Overwrite Am in ConData (this changes in case of STDP, otherwise it stays the same)
if simdata.STDP
    ConData.PMat_AlltoAll.Am = final_connectivity;
end

