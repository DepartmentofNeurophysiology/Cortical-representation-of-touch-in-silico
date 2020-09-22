function [ lummat, calcmat ] = simulation_to_calcium(CalciumPara, savefolder, savename, initdata, Nsim, varargin)
% Loads the results of simulations, and generates calcium and luminescence
% data
% INPUT:
% * Para struct with with fields (see Vogelstein et al 2009)
%     * frame_rate_c        % : sampling rate calcium calculations (kHz)
%     * tau_c               % : calcium decay time constant (ms)
%     * Ca_b                % : baseline calcium concentration (microM)
%     * A_c                 % : calcium concentration 'jump' for each spike
%     * sigma_c             % : std noise for calcium signal
%     * alpha               % : scale fluorescence signal
%     * beta                % : offset fluorescence signal
%     * sigma_f             % : std noise for luminescence
%     * tmax                % (optional) maximal time for luminescience calculation
% * savefolder: folder to save result
% * savename: name to save result
% * initdata (see simulations)
% * Nsim: # simulations
% * 'simvec' (optional): vector with which simulations to use

%% Initialize
simvec = 1:Nsim;
% Look for 'varargin' inputs
len = length(varargin);
% check "len" for even number
if mod(len,2) > 0
    error('Wrong arguments: must be name-value pairs.');
end
for i = 1:2:len
    switch lower(varargin{i})
        case 'simvec'
            simvec=varargin{i+1};
            Nsim = length(simvec);
        otherwise
            % neglect invalid option
            disp(['Ignoring invalid input ' varargin{i}])
    end
end

if initdata.Vthresdyn
    thresholdname = ['dynthreshold_set'];
else
    thresholdname = ['fixthreshold_set'];
end

%% Preallocate
Nlumtime = ceil(CalciumPara.tmax*CalciumPara.frame_rate_c);

%% Loop over sims and cells to calculate calcium traces
for ns = 1:Nsim
    disp(['Calculating calcium traces for simulation ', num2str(simvec(ns))])
    load([savefolder savename '_Simcolumn_' thresholdname '_' initdata.setVthres.type '_simulation_' num2str(simvec(ns))], 'modelspt');
    [NAll, ~] = size(modelspt);
    lummat = zeros(NAll, Nlumtime+1);
    calcmat = zeros(NAll, Nlumtime+1);
    for nn = 1:NAll
        if nn==1
            [calcmat_temp, lummat_temp, ~] = spike_train_to_calc_lum(modelspt(nn,:), CalciumPara);
            Nlumtime = length(lummat_temp);
            calcmat(nn,1:Nlumtime) = calcmat_temp;
            lummat(nn,1:Nlumtime) = lummat_temp;  
            lummat = lummat(:, 1:Nlumtime);
            calcmat = calcmat(:, 1:Nlumtime);
        else
            [calcmat(nn,:), lummat(nn,:), ~] = spike_train_to_calc_lum(modelspt(nn,:), CalciumPara);
        end
    end
    save([savefolder savename '_Simcolumn_' thresholdname '_' initdata.setVthres.type '_simulation_' num2str(simvec(ns)) '_calciumdata'],'lummat', 'calcmat', 'CalciumPara', 'simvec','-v7.3')
end


