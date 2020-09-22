function [calc, lum, usedseed] = spike_train_to_calc_lum(spiketimes, Para, varargin)
% Calculate calcium and luminescence trace on the basis of neural spikes
% (see Vogelstein et al 2009, eq 1 and 2)

% INPUT
% * spiketimes (array with spike times)
% * Para struct with fields
%       * frame_rate_c 
%       * tau_c : calcium decay time constant
%       * Ca_b: baseline calcium concentration
%       * A_c: calcium concentration 'jump' for each spike
%       * sigma_c: std noise for calcium signal
%       * alpha: scale fluorescence signal
%       * beta: offset fluorescence signal
%       * sigma_f: std noise for luminescence
%       * tmax (optional): maximum time for calculation
% * seed (optional) for random number generator (optional)

% NB Make sure everything is in the same units (i.e. ms and kHz or s and Hz for spike
% times, tau, frame_rate_c)!


%% Initialize
% Default values
rng('shuffle');
scurr = rng;
usedseed = scurr.Seed;
plotyn = 0;

% Look for 'varargin' inputs
len = length(varargin);
% check "len" for even number
if mod(len,2) > 0
    error('Wrong arguments: must be name-value pairs.');
end
for i = 1:2:len
    switch lower(varargin{i})
        case 'seed'
            usedseed=varargin{i+1};
            rng(usedseed)
        case 'plotyn'
            plotyn=varargin{i+1};
        otherwise
            % neglect invalid option
            disp(['Ignoring invalid input ' varargin{i}])
    end
end

if ~isfield(Para, 'tmax')
    Para.tmax = max(spiketimes);
end
bin_c = 1/Para.frame_rate_c;
NT = round(Para.tmax/bin_c);

%% Preallocate
calc = zeros(NT+1,1);
calc(1) = Para.Ca_b+Para.sigma_c*sqrt(bin_c)*randn;
lum = zeros(NT+1,1);
lum(1) = Para.alpha*calc(1)+Para.beta+Para.sigma_f*randn;

%% Calculate calcium signal
for nt = 1:NT
    tstart = (nt-1)*bin_c;
    tend = nt*bin_c;
    try
        nspikes = sum(spiketimes(find(spiketimes>tstart))<=tend);
    catch
        keyboard
    end
    delta_c = -(bin_c/Para.tau_c)*(calc(nt)-Para.Ca_b)+Para.A_c*nspikes + Para.sigma_c*sqrt(bin_c)*randn;
    calc(nt+1) = calc(nt)+delta_c; 
    lum(nt+1) = Para.alpha*calc(nt+1)+Para.beta+Para.sigma_f*randn;
end

%% Plot for checking
% if max(spiketimes)>0
%     keyboard
% end
if plotyn
    time = (0:NT)*bin_c;
    figure
    subplot(2,1,1)
    plot(spiketimes, max(calc)+ones(size(spiketimes)),'*')
    hold all
    plot(time, calc)
    legend('spikes','calcium')
    subplot(2,1,2)
    plot(time, lum)
    xlabel('time')
    title('\Delta F / F')
    ax = gca;
    ax.XTick = (1:NT)*bin_c;
    grid on
end