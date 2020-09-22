function ConData = reorganize_conmat(fname, savefolder, CMDM_folder, CMDMs_file, includemodulationyn, varargin)
% reorganize_conmat 
% * reorganizes the connectivity into a single list in order to speed up
% simulations. This all gets put into the structure 'ConData' that is needed
% to run each trial (function: simcolumn_indtrial)
% * If the files '_ParaMat', '_ParaMat_reduced.mat' and
% '_WhiskerModulationModel.mat' do not exist in the folder of the
% connectivity file, it will create them and save them there (based on
% random numbers, so different each time you create one). 

% INPUTS:
% * fname: name to save data file
% * savefolder: foldername to save data files (should end with '/')
% * CMDMs_folder: connectivity data folder (should end with '/')
% * CMDMs_file: connectivity data file name
% * includemodulationyn: 1 = include direct whisker modulation (Crochet),
% otherwise 0.

% OUTPUT:
% * Data structure to be used by run_save_trials



%% check optional parameters
P = parsePairs(varargin);
checkField(P, 'InputLabel', 'Mat');
checkField(P, 'WhPara', {});


%% load the connectivity matrix file
% the file should contain CM, DM, cellinfo, cellpara (for Izhikevich model)
ConData.filename = CMDMs_file;
try 
%     load([CMDM_folder CMDMs_file]);
    load([CMDM_folder CMDMs_file], '-regexp',  '^(?!savefolder)\w');
catch
    disp('Please specify a connectivity file. It can be constructed after the example of cellinfo_newbarrel')
    keyboard
end

%% generate synaptic parameter matrix
% load existing matrix if it is present

% if isempty(dir([CMDM_folder CMDMs_file(1:end-4) '_ParaMat.mat']))==0
%     load([CMDM_folder CMDMs_file(1:end-4) '_ParaMat.mat']);
if isempty(dir([CMDM_folder CMDMs_file '_ParaMat.mat']))==0
    load([CMDM_folder CMDMs_file '_ParaMat']);
else
    disp('Generating new synapse parameters Thalamus - L4 - L23')
    % use function simcolumn_ParaMat_ver2
    % thalamus to L4   
    Pthto4 = Synapse_parameter_ThtoL4_func;
    ParaMat_thto4 = simcolumn_ParaMat_ver2(Pthto4, thainfo, l4info, CMThtol4, DMThtol4);
    % thalamus to L23
    Pthto2 = Synapse_parameter_ThtoL23_func;
    ParaMat_thto2 = simcolumn_ParaMat_ver2(Pthto2, thainfo, l23info, CMThtol23, DMThtol23);
    % L4 to L4
    P4to4 = Synapse_parameter_L4toL4_func;
    ParaMat_4to4 = simcolumn_ParaMat_ver2(P4to4, l4info, l4info, CMl4tol4, DMl4tol4);
    % L4 to L23
    P4to2 = Synapse_parameter_L4toL23_func;
    ParaMat_4to2 = simcolumn_ParaMat_ver2(P4to2, l4info, l23info, CMl4tol23, DMl4tol23);
    % L23 to L23
    P2to2 = Synapse_parameter_L23toL23_func;
    ParaMat_2to2 = simcolumn_ParaMat_ver2(P2to2, l23info, l23info, CMl23tol23, DMl23tol23);
    % L23 to L4
    
    
    %% assemble CMs and DMs used for running the simulation
    % generate two pairs of matrix, one is input layer (neurons used to
    % provide input but not simulated) to simulated layers, and the other is
    % simulated layers to itself
    CM_IntoAll = sparse([CMThtol4; CMThtol23]);
    CM_AlltoAll = simcolumn_assemble_ParaMat({l4info, l23info}, {CMl4tol4, CMl4tol23, ...
        [], CMl23tol23});
    Am_IntoAll = sparse([ParaMat_thto4.Am; ParaMat_thto2.Am]);
    Am_AlltoAll = simcolumn_assemble_ParaMat({l4info, l23info}, {ParaMat_4to4.Am, ...
        ParaMat_4to2.Am, [], ParaMat_2to2.Am});
    Trise_IntoAll = sparse([ParaMat_thto4.Trise; ParaMat_thto2.Trise]);
    Trise_AlltoAll = simcolumn_assemble_ParaMat({l4info, l23info}, {ParaMat_4to4.Trise, ...
        ParaMat_4to2.Trise, [], ParaMat_2to2.Trise});
    Tfall_IntoAll = sparse([ParaMat_thto4.Tfall; ParaMat_thto2.Tfall]);
    Tfall_AlltoAll = simcolumn_assemble_ParaMat({l4info, l23info}, {ParaMat_4to4.Tfall, ...
        ParaMat_4to2.Tfall, [], ParaMat_2to2.Tfall});
    Fail_IntoAll = sparse([ParaMat_thto4.Fail; ParaMat_thto2.Fail]);
    Fail_AlltoAll = simcolumn_assemble_ParaMat({l4info, l23info}, {ParaMat_4to4.Fail, ...
        ParaMat_4to2.Fail, [], ParaMat_2to2.Fail});
    Plas_IntoAll = sparse([ParaMat_thto4.Plas; ParaMat_thto2.Plas]);
    Plas_AlltoAll = simcolumn_assemble_ParaMat({l4info, l23info}, {ParaMat_4to4.Plas, ...
        ParaMat_4to2.Plas, [], ParaMat_2to2.Plas});
    CV_IntoAll = sparse([ParaMat_thto4.CV; ParaMat_thto2.CV]);
    CV_AlltoAll = simcolumn_assemble_ParaMat({l4info, l23info}, {ParaMat_4to4.CV, ...
        ParaMat_4to2.CV, [], ParaMat_2to2.CV});
    Delay_IntoAll = sparse([ParaMat_thto4.Delay; ParaMat_thto2.Delay]);
    Delay_AlltoAll = simcolumn_assemble_ParaMat({l4info, l23info}, {ParaMat_4to4.Delay, ...
        ParaMat_4to2.Delay, [], ParaMat_2to2.Delay});
    
    
%     % only simulate L4 to L23
%     CM_IntoAll = sparse([CMl4tol23]);
%     CM_AlltoAll = simcolumn_assemble_ParaMat({l23info}, {CMl23tol23});
%     Am_IntoAll = sparse([ParaMat_4to2.Am]);
%     Am_AlltoAll = simcolumn_assemble_ParaMat({l23info}, {ParaMat_2to2.Am});
%     Trise_IntoAll = sparse([ParaMat_4to2.Trise]);
%     Trise_AlltoAll = simcolumn_assemble_ParaMat({l23info}, {ParaMat_2to2.Trise});
%     Tfall_IntoAll = sparse([ParaMat_4to2.Tfall]);
%     Tfall_AlltoAll = simcolumn_assemble_ParaMat({l23info}, {ParaMat_2to2.Tfall});
%     Fail_IntoAll = sparse([ParaMat_4to2.Fail]);
%     Fail_AlltoAll = simcolumn_assemble_ParaMat({l23info}, {ParaMat_2to2.Fail});
%     Plas_IntoAll = sparse([ParaMat_4to2.Plas]);
%     Plas_AlltoAll = simcolumn_assemble_ParaMat({l23info}, {ParaMat_2to2.Plas});
%     CV_IntoAll = sparse([ParaMat_4to2.CV]);
%     CV_AlltoAll = simcolumn_assemble_ParaMat({l23info}, {ParaMat_2to2.CV});
%     Delay_IntoAll = sparse([ParaMat_4to2.Delay]);
%     Delay_AlltoAll = simcolumn_assemble_ParaMat({l23info}, {ParaMat_2to2.Delay});
    
    
    % put these data input structure
    PMat_IntoAll.CM = CM_IntoAll;
    PMat_IntoAll.Am = Am_IntoAll;
    PMat_IntoAll.Trise = Trise_IntoAll;
    PMat_IntoAll.Tfall = Tfall_IntoAll;
    PMat_IntoAll.Plas = Plas_IntoAll;
    PMat_IntoAll.Fail = Fail_IntoAll;
    PMat_IntoAll.CV = CV_IntoAll;
    PMat_IntoAll.Delay = Delay_IntoAll;
    
    PMat_AlltoAll.CM = CM_AlltoAll;
    PMat_AlltoAll.Am = Am_AlltoAll;
    PMat_AlltoAll.Trise = Trise_AlltoAll;
    PMat_AlltoAll.Tfall = Tfall_AlltoAll;
    PMat_AlltoAll.Plas = Plas_AlltoAll;
    PMat_AlltoAll.Fail = Fail_AlltoAll;
    PMat_AlltoAll.CV = CV_AlltoAll;
    PMat_AlltoAll.Delay = Delay_AlltoAll;
    
    % parameter for neuronal model
    Neuron_Para = [l4para; l23para];
%     Neuron_Para = l23para;
    
    % save the parameter matrix
%     save([CMDM_folder CMDMs_file(1:end-4) '_ParaMat.mat'], 'PMat_IntoAll', 'PMat_AlltoAll', 'Neuron_Para');
    save([CMDM_folder CMDMs_file '_ParaMat.mat'], 'PMat_IntoAll', 'PMat_AlltoAll', 'Neuron_Para');

end

%% some testing routines
% set th to l23 connection to 0
% PMat_IntoAll.CM(size(l4info,1)+1:end, :) = 0;
% set l4 inhibitory to l23 connection to 0
% temp = PMat_AlltoAll.CM(size(l4info, 1) + 1:end, 1:size(l4info, 1));
% temp(temp < 0) = 0;
% PMat_AlltoAll.CM(size(l4info, 1) + 1:end, 1:size(l4info, 1)) = temp;
% set L4 to L23 connection to 0
% PMat_AlltoAll.CM(size(l4info,1)+1:end, 1:size(l4info, 1)) = 0;
% set L23 to L23 connection to 0
% PMat_AlltoAll.CM(size(l4info,1)+1:end, size(l4info, 1)+1:end) = 0;


%% assign spike threshold value for different type of neuron
% the values are mean+-std; individual incidents are generated within
% individual simulation runs
th_l4 = [-41.2, -41.2, -43.7, -43.7];
th_l4_std = [3,3,3,3];
th_l4_sl = [1, 1, 1, 1];
th_l4_a = [-38, -38, -41, -40];
th_l23 = [-43, -42, -39, -43, -46.27, -46, -43.8, -41, -43.05, -43, -46];
th_l23_std = [3, 3, 3, 2.7, 2.92, 3, 2.8, 3, 2.89, 2.89, 3];
th_l23_sl = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1];
th_l23_a = [-40, -38, -37, -41, -43, -41, -39, -39, -38, -39, -41];
thresh.avg = [th_l4, th_l23];
thresh.std = [th_l4_std, th_l23_std];
thresh.sl = [th_l4_sl, th_l23_sl];
thresh.a = [th_l4_a, th_l23_a];

% Fontaine 2014 model for spike threshold 
th_l4_Vmin = [-37, -37, -39, -39];
th_l4_Vi = [-41, -41, -43, -43];
th_l4_alpha = [0.3, 0.3, 0.3, 0.3];
th_l4_Ka = [7, 7, 7, 7];
th_l4_beta = [1.1, 1.1, 1.1, 1.1];
th_l4_Ki = th_l4_Ka./(th_l4_beta - th_l4_alpha);
th_l4_tau = [6.5, 6.5, 6.5, 6.5];

th_l23_Vmin = [-39.5, -38.5, -36, -39, -41, -41, -39.5, -37, -39, -39, -41];
th_l23_Vi = [-43, -42.5, -41, -43, -45, -45, -43, -41, -43, -43, -45];
th_l23_alpha = 0.3*ones(size(th_l23_Vmin));
th_l23_Ka = 7*ones(size(th_l23_Vmin));
th_l23_beta = 1.1*ones(size(th_l23_Vmin));
th_l23_Ki = th_l23_Ka./(th_l23_beta - th_l23_alpha);
th_l23_tau = 6.5*ones(size(th_l23_Vmin));

thresh.Vmin = [th_l4_Vmin, th_l23_Vmin];
thresh.Vi = [th_l4_Vi, th_l23_Vi];
thresh.alpha = [th_l4_alpha, th_l23_alpha];
thresh.Ka = [th_l4_Ka, th_l23_Ka];
thresh.beta = [th_l4_beta, th_l23_beta];
thresh.Ki = [th_l4_Ki, th_l23_Ki];
thresh.tau = [th_l4_tau, th_l23_tau];

% %only simulate L4 to L23
% thresh.avg = [th_l23];
% thresh.std = [th_l23_std];

%% reduce matrix size to speed up simulation
% if isempty(dir([CMDM_folder CMDMs_file(1:end-4) '_ParaMat_reduced.mat']))==0
%     load([CMDM_folder CMDMs_file(1:end-4) '_ParaMat_reduced.mat']);
if isempty(dir([CMDM_folder CMDMs_file '_ParaMat_reduced.mat']))==0
    load([CMDM_folder CMDMs_file '_ParaMat_reduced.mat']);    
else
    disp('Generating new reduced parameter matrix')
    %% adjust cell index in cellinfo
    cellinfo_In = thainfo;
    l23info(:,4) = l23info(:,4) + max(l4info(:,4));
    cellinfo_All = [l4info; l23info];
    
    % %only simulate L4 to L23
    % l23info(:, 4) = l23info(:,4) - max(l4info(:,4));
    % cellinfo_In = l4info;
    % cellinfo_All = l23info;
    
    NIn = size(cellinfo_In, 1);
    NtIn = [];
    Nt = unique(cellinfo_In(:,4));
    for i = 1:length(Nt)
        NtIn(Nt(i)) = length(find(cellinfo_In(:,4) == Nt(i)));
    end
    NAll = size(cellinfo_All, 1);
    NtAll = [];
    Nt = unique(cellinfo_All(:,4));
    for i = 1:length(Nt)
        NtAll(Nt(i)) = length(find(cellinfo_All(:,4) == Nt(i)));
    end
    
    %% convert to homogeneous time constant condition
    % to avoid using isnan function, replace 0 entries in Trise and Tfall
    % matrix
    % use homogenerous time constant for each type of connection
    % convert time constant matrix into homogeneous form, which should have
    % size of (Npost, Ntype_pre)
    temp = simcolumn_convert2HomoTau(PMat_IntoAll.Trise, NtIn, NtAll);
    temp(temp == 0) = 1;
    PMat_IntoAll.Trise = temp;
    temp = simcolumn_convert2HomoTau(PMat_IntoAll.Tfall, NtIn, NtAll);
    temp(temp == 0) = 2;
    PMat_IntoAll.Tfall = temp;
    temp = simcolumn_convert2HomoTau(PMat_AlltoAll.Trise, NtAll, NtAll);
    temp(temp == 0) = 1;
    PMat_AlltoAll.Trise = temp;
    temp = simcolumn_convert2HomoTau(PMat_AlltoAll.Tfall, NtAll, NtAll);
    temp(temp == 0) = 2;
    PMat_AlltoAll.Tfall = temp;
    
    % keyboard
    %% matrix size reduction
    [PMat_IntoAll.CM, PMat_IntoAll.preCell, PMat_IntoAll.preIdx, PMat_IntoAll.postIdx]...
        = simcolumn_CMreduce_all(full(PMat_IntoAll.CM), cellinfo_In);
    [PMat_AlltoAll.CM, PMat_AlltoAll.preCell, PMat_AlltoAll.preIdx, PMat_AlltoAll.postIdx]...
        = simcolumn_CMreduce_all(full(PMat_AlltoAll.CM), cellinfo_All);
    
    PMat_IntoAll.STD = ones(size(PMat_IntoAll.CM));  % STD is used to model short term dynamics
    PMat_AlltoAll.STD = ones(size(PMat_AlltoAll.CM));
    
    PMat_IntoAll.Am = simcolumn_PMatreduce_all(full(PMat_IntoAll.Am), PMat_IntoAll.CM);
    PMat_IntoAll.CV = simcolumn_PMatreduce_all(full(PMat_IntoAll.CV), PMat_IntoAll.CM);
    PMat_IntoAll.Fail = simcolumn_PMatreduce_all(full(PMat_IntoAll.Fail), PMat_IntoAll.CM);
    % PMat_IntoAll.Trise = simcolumn_PMatreduce_all(full(PMat_IntoAll.Trise), PMat_IntoAll.CM);
    % PMat_IntoAll.Tfall = simcolumn_PMatreduce_all(full(PMat_IntoAll.Tfall), PMat_IntoAll.CM);
    PMat_IntoAll.Delay = simcolumn_PMatreduce_all(full(PMat_IntoAll.Delay), PMat_IntoAll.CM);
    PMat_IntoAll.Plas = simcolumn_PMatreduce_all(full(PMat_IntoAll.Plas), PMat_IntoAll.CM);
    PMat_IntoAll.preID = simcolumn_preID_ver2(PMat_IntoAll.CM, cellinfo_In, cellinfo_All);
    
    PMat_AlltoAll.Am = simcolumn_PMatreduce_all(full(PMat_AlltoAll.Am), PMat_AlltoAll.CM);
    PMat_AlltoAll.CV = simcolumn_PMatreduce_all(full(PMat_AlltoAll.CV), PMat_AlltoAll.CM);
    PMat_AlltoAll.Fail = simcolumn_PMatreduce_all(full(PMat_AlltoAll.Fail), PMat_AlltoAll.CM);
    % PMat_AlltoAll.Trise = simcolumn_PMatreduce_all(full(PMat_AlltoAll.Trise), PMat_AlltoAll.CM);
    % PMat_AlltoAll.Tfall = simcolumn_PMatreduce_all(full(PMat_AlltoAll.Tfall), PMat_AlltoAll.CM);
    PMat_AlltoAll.Delay = simcolumn_PMatreduce_all(full(PMat_AlltoAll.Delay), PMat_AlltoAll.CM);
    PMat_AlltoAll.Plas = simcolumn_PMatreduce_all(full(PMat_AlltoAll.Plas), PMat_AlltoAll.CM);
    PMat_AlltoAll.preID = simcolumn_preID_ver2(PMat_AlltoAll.CM, cellinfo_All, cellinfo_All);
    
    % keyboard
    % to avoid using isnan function, replace 0 entries in Trise and Tfall
    % matrix
    PMat_IntoAll.Trise(PMat_IntoAll.Trise == 0) = 1;
    PMat_IntoAll.Tfall(PMat_IntoAll.Tfall == 0) = 2;
    PMat_AlltoAll.Trise(PMat_AlltoAll.Trise == 0) = 1;
    PMat_AlltoAll.Tfall(PMat_AlltoAll.Tfall == 0) = 2;
    
    % now we can have normalization factor in double-exponential function as
    % constant
    PMat_IntoAll.CN = (PMat_IntoAll.Trise./PMat_IntoAll.Tfall).^...
        (PMat_IntoAll.Tfall./(PMat_IntoAll.Trise - PMat_IntoAll.Tfall));
    PMat_AlltoAll.CN = (PMat_AlltoAll.Trise./PMat_AlltoAll.Tfall).^...
        (PMat_AlltoAll.Tfall./(PMat_AlltoAll.Trise - PMat_AlltoAll.Tfall));
    % keyboard
    %% running variable for the dynamic system to calculate synaptic current
    % size of these matrix will depend on if using homogeneous tau or hetero
    % form
    % in case of homogeneous, will be size of Npost-by-Ntypepre; for
    % heterogeneouse, will be Npost-by-Npre
    PMat_IntoAll.g = zeros(size(PMat_IntoAll.Trise));   % g, s are needed to model double exponential function
    PMat_IntoAll.s = zeros(size(PMat_IntoAll.Trise));
    PMat_AlltoAll.g = zeros(size(PMat_AlltoAll.Trise));
    PMat_AlltoAll.s = zeros(size(PMat_AlltoAll.Trise));
    
    %% check which pre-cell indexing method to use in the simulation
    % use preID when homogeneouse tau is used; otherwise use preCell
    % also need to get an array to indicate whether the presynaptic cell is
    % excitatory or inhibitory
    if size(PMat_AlltoAll.g, 2) == size(cellinfo_All, 1)
        % heterogeneouse tau
        PMat_IntoAll.preID = PMat_IntoAll.preCell;
        PMat_AlltoAll.preID = PMat_AlltoAll.preCell;
        PMat_IntoAll.preEI = PMat_IntoAll.preCell(:,3);
        PMat_AlltoAll.preEI = PMat_AlltoAll.preCell(:,3);
        
    elseif size(PMat_AlltoAll.g, 2) == length(NtAll)
        % homogeneouse tau
        cNt = cumsum(NtIn);
        for i = 1:length(NtIn)
            PMat_IntoAll.preEI(i) = cellinfo_In(cNt(i), 6);
        end
        cNt = cumsum(NtAll);
        for i = 1:length(NtAll)
            PMat_AlltoAll.preEI(i) = cellinfo_All(cNt(i), 6);
        end
        
    else
        error('size of matrix g,s is not correct')
    end
    
%     save([CMDM_folder CMDMs_file(1:end-4) '_ParaMat_reduced.mat'], 'PMat_IntoAll', 'PMat_AlltoAll', 'Neuron_Para', 'cellinfo_In', 'cellinfo_All');
    save([CMDM_folder CMDMs_file '_ParaMat_reduced.mat'], 'PMat_IntoAll', 'PMat_AlltoAll', 'Neuron_Para', 'cellinfo_In', 'cellinfo_All');
end

%% Whisker modulation on S1 neurons
% the model is WI(t) = k*A(t)*sin(pha(t)-pha0) + RetraSet - mean(A) -min(A)
%   A(t): whisking amplitude at time t; 
%   k: modulation depth (see Crochet et al 2011)
%   pha(t): whisking phase at time t
%   pha0: prefered phase of the given neuron; evenly distributed at [-pi, pi]
%   RetraSet: set point (most retracted point) in each circle
% absolute positional information is provided when RetraSet is taking into
% account. otherwise absolute whisker position is lost
%
% assuming: Whisker angle correlated inputs to L2/3 populaion; not all 
% neurons are modulated by whisking and the distribution of whisking 
% modulation is depth dependent (90% in L3 and 20% in L2); whisker 
% modulated Neurons have evenly distributed phase preference; modulation
% strength is log-normally distributed with 0.08 +- 0.02 mV/deg; additional
% white noise inputs to reduce the r2 (in petersen paper the r2 is low; not
% implemented now)
% SST neurons (type 9 and 10) are not modulated
%
% modulation depth
% significantly modulated neurons
if includemodulationyn
%     if isempty(dir([CMDM_folder CMDMs_file(1:end-4) '_WhiskerModulationModel.mat']))==0
%         load([CMDM_folder CMDMs_file(1:end-4) '_WhiskerModulationModel.mat']);
%     if isempty(dir([CMDM_folder CMDMs_file '_WhiskerModulationModel.mat']))==0
%         load([CMDM_folder CMDMs_file '_WhiskerModulationModel.mat'], '-regexp',  '^(?!savefolder)\w');
    if exist([CMDM_folder CMDMs_file '_WhiskerModulationModel.mat'], 'file') == 2
        disp('Found whisker modulation model');
        load([CMDM_folder CMDMs_file '_WhiskerModulationModel.mat']);
    else
        disp('Connectivity file does not contain direct L23 modulation by motor cortex')
        s = input('Do you want to include direct L23 modulation? (y/n)','s');
        if strcmp(s,'y')
            disp('Generating modulation data')
            s = input('Do you want to seed the rng? (y/n)', 's');
            if strcmp(s, 'y')
                ConData.modulationseed = input('Please give a seed');
                rng(ConData.modulationseed)
            else
                rng('shuffle')
                scurr = rng;
                ConData.modulationseed = scurr.Seed;
            end

            
            % find main barrel: only cells in principal barrel are
            % modulated
            [Nbx, Nby] = size(barrelstruct);
            mainbarrelfound = 0;
            nb = 0;
            for nbx = 1:Nbx
                for nby = 1:Nby
                    nb = nb+1;
                    if barrelstruct{nbx,nby}.mainbarrel == 1
                        mainbarrelfound = mainbarrelfound+1;
                        mainbarrelnr = nb;
                    end
                end
            end
            if mainbarrelfound>1
                error('Multiple main barrels')
            elseif mainbarrelfound==0
                error('No main barrel defined')
            end
            
            Midx = [];
            idx = find((cellinfo_All(:,5) == mainbarrelnr) & ((cellinfo_All(:,4)>4) & (cellinfo_All(:,3)<430/3))); %L2 cells
            tidx = randperm(length(idx));
            tidx = idx(tidx(1:round(0.2*length(idx))));
            Midx = tidx;
            idx = find((cellinfo_All(:,5) == mainbarrelnr) & ((cellinfo_All(:,4)>4) & (cellinfo_All(:,3)>=430/3))); %L3 cells
            tidx = randperm(length(idx));
            tidx = idx(tidx(1:round(0.9*length(idx))));
            Midx = [Midx; tidx];
            % generate modulation depth k
            k = zeros(size(cellinfo_All, 1), 1);
            k(Midx) = simcolumn_generatePara_lognorm([0.08, 0.02], size(Midx));
            % prefered phase, evenly distributed between [-pi, pi]
            Pha0 = 2*pi*rand(size(k)) - pi;
            WhMod.k = k;
            WhMod.k(cellinfo_All(:,4) == 9) = 0;
            WhMod.k(cellinfo_All(:,4) == 10) = 0;
            WhMod.Pha0 = Pha0;
%             save([CMDM_folder CMDMs_file(1:end-4) '_WhiskerModulationModel.mat'], 'WhMod','modulationseed');
            modulationseed = ConData.modulationseed;
            save([CMDM_folder CMDMs_file '_WhiskerModulationModel.mat'], 'WhMod','modulationseed');
        end
    end
end
%%
ConData.FnametoSave = fname;
ConData.savefolder = savefolder;
ConData.PMat_IntoAll = PMat_IntoAll;
ConData.PMat_AlltoAll = PMat_AlltoAll;
ConData.Neuron_Para = Neuron_Para;
ConData.Neuron_Vt = thresh;
ConData.Cellinfo_In = cellinfo_In;
ConData.Cellinfo_All = cellinfo_All;
if isnumeric(P.InputLabel)
    ConData.InputLabel = num2str(P.InputLabel);
elseif ischar(P.InputLabel)
    ConData.InputLabel = P.InputLabel;
end

if exist('WhMod', 'var')
    ConData.WhiskModel = WhMod;
end

%% number of cells in each group, i.e. input and simulated
NIn = size(ConData.Cellinfo_In, 1);
NtIn = [];
Nt = unique(ConData.Cellinfo_In(:,4));
for i = 1:length(Nt)
    NtIn(Nt(i)) = length(find(ConData.Cellinfo_In(:,4) == Nt(i)));
end
NAll = size(ConData.Cellinfo_All, 1);
NtAll = [];
Nt = unique(ConData.Cellinfo_All(:,4));
for i = 1:length(Nt)
    NtAll(Nt(i)) = length(find(ConData.Cellinfo_All(:,4) == Nt(i)));
end

ConData.NIn = NIn;
ConData.NtIn = NtIn;
ConData.NAll = NAll;
ConData.NtAll = NtAll;
ConData.Nt = Nt;

%% Save
a = exist(ConData.savefolder, 'dir');
if ~(a==7) 
    mkdir(ConData.savefolder)
end
% save([savefolder CMDMs_file(1:end-4) '_ConData.mat'], 'Data');
save([savefolder CMDMs_file '_ConData.mat'], 'ConData');


