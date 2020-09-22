function SpikeTrainStruct = kernel_recording_to_spiketrain(WhiskerStruct, KernelStruct, SpikeGenStruct, savename, plotyn)
% generate PSTHs and spike trains from whisker recordings and kernels 
% Recordings and kernels can be any unit, as long as kernels and traces
% have the same unit
% Recordings and kernels can be multidimensional (for example simultaneous
% base angle and deflection angle), but activation functions need to be
% defined accordingly!


% INPUT
% * WhiskerStruct: struct with recordings, fields:
%       .Recording: NdimxNtrace cell array (all dimensions should have the
%       same Ntrace traces, and for each trace the two dims should have the same length Ntime each)
%       .binsize: in ms (should also be the same for all dimensions)
% * KernelStruct
%       .Kernels{Nkernel, Ndim}: matrix of size (Nkernel x Nkerneltime) with kernels 
%       .ActivationFunction: function handle with activation function (of
%       all dimensions combined). NB this function should be on the path!
%           .Params{Nkernel, Ndim}: parameters of the activation function
%       .kerneltime: time array (in ms), so that it is clear what is the causal and
%       what the acausal part of the kernel 
% * SpikeGenStruct (Poisson spike train generation):
%       *.refra             = refractory period (ms)
%       *.Ntrial_pertrace   = # realizations for each recording
%       *.binsize           = binsize spike trains (ms)
%       *.delay             = (ms)
%       *.scaling     
%       *.seed              = (optional) Ntrace x Ntrial_pertrace array with seeds for Poisso spiking
% * plotyn (optional, 1 or 0): if 1: make plots for checking


% OUTPUT
% * SpikeTrainStruct
%       .PSTH: (Nkernel x Ntrace) cell, with in each entry (1 x Ntime) array with a PSTH (# spikes / sec) for each cell for this deflection trace; 
%       .SpikeTimes: (Nkernel x Ntrace) cells with each a (1xNtrial_pertrace) cell with a with spike times (ms)
%       .SpikeCount: (Nkernel x Ntrace) cell
%       .seed: 

% NB Needs a separate activation function on the path! (handle passed here)
% This activation function should look like the following: 
% PSTH = activation_function(ConvTrace, Params)
% Unit PSTH = (# spikes / sec). The probability of a spike for every time
% step Pspike = PSTH*binsize/1000;
% Input to the activation function: 
% * ConvTrace: ConvTrace = NdimxNtime array with a recording convolved with a kernel
% * Params{dimension}: Parameters of the activation function for the
% kernels for this dimension

if nargin == 3
    plotyn = 0;
    savename = [];
elseif nargin == 4
    if isnumeric(savename)
        plotyn = savename;
        savename = [];
    else
        plotyn = 0;
    end
end

[Ndimw, Ntrace] = size(WhiskerStruct.Recording);
[Nkernel, Ndimk] = size(KernelStruct.Kernels);
binsize_kernels = KernelStruct.kerneltime(2)-KernelStruct.kerneltime(1);
%% Check
if Nkernel < Ndimk
    ip = input([num2str(Nkernel) ' kernels with ' num2str(Ndimk) ' dimensions each. Is this ok, or transpose? (ok/tp)'], 's');
    if strcmp(ip, 'tp')
        KernelStruct.Kernels = KernelStruct.Kernels';
        [Nkernel, Ndimk] = size(KernelStruct.Kernels);
    end
end
        

if ~(Ndimk == Ndimw)
    error('Kernels should have the same dimensions as recordings')
end

for nt = 1:Ntrace
    for nd = 1:Ndimw
        if ~(length(WhiskerStruct.Recording{1,nt}) == length(WhiskerStruct.Recording{nd,nt}))
            error(['Error in trace ' num2str(nt) ', dimension ' num2str(nd) ': Make sure the length of the recording is the same for all dimensions'])
        end
    end
end

for nd = 1:Ndimk
    for nk = 1:Nkernel
        if ~(length(KernelStruct.Kernels{nk,nd}) == length(KernelStruct.kerneltime))
            error(['Error in kernel ' num2str(nk) ': Make sure all kernels have the same length as the kerneltime array'])
        end
    end
end


if WhiskerStruct.binsize == binsize_kernels
    disp('Sampling rate kernels the same as whisker recordings')
else
    if WhiskerStruct.binsize > binsize_kernels
        disp('Sampling rate kernels higher than whisker recordings: upsample whisker recordings')
    elseif WhiskerStruct.binsize < binsize_kernels
        disp('Sampling rate kernels lower than whisker recordings: downsample whisker recordings')
    end
    for nd = 1:Ndimw
        for nt = 1:Ntrace
            f1 = 1/binsize_kernels;
            f2 = 1/WhiskerStruct.binsize; 
            while ~((rem(f1,1)==0) && (rem(f2,1)==0))
                f1 = f1*10;
                f2 = f2*10;
            end
            WhiskerStruct.Recording{nd,nt} = resample(WhiskerStruct.Recording{nd,nt},f1 ,f2);
        end
    end
    WhiskerStruct.binsize = binsize_kernels;
end

%% Convolve recordings with relevant kernel
la = sum(KernelStruct.kerneltime<0);
lc = sum(KernelStruct.kerneltime>0);
ConvTrace = cell(Nkernel, Ntrace);
for nk = 1:Nkernel
    for nt = 1:Ntrace
        ConvTrace{nk,nt} = nan*ones(Ndimw, length(WhiskerStruct.Recording{nd,nt}));
        for nd = 1:Ndimw  
            ConvTrace{nk,nt}(nd,:) = convolve_kernel_acausal( WhiskerStruct.Recording{nd,nt}, KernelStruct.Kernels{nk, nd}, la, lc);
        end
    end
end
%% Apply activation function 
SpikeTrainStruct.PSTH = cell(Nkernel, Ntrace);
for nk = 1:Nkernel
    for nt = 1:Ntrace
        SpikeTrainStruct.PSTH{nk, nt} = feval(KernelStruct.ActivationFunction.function, ConvTrace{nk,nt}, KernelStruct.ActivationFunction.Params{nk});
    end
end

%% Make PSTHs into spike trains
% assuming Possion spike trains
if isfield(SpikeGenStruct, 'seed')
    [Ntrace_seed, Ntrial_pertrace_seed]=size(SpikeGenStruct.seed);
    if (Ntrace_seed == Ntrace) && (Ntrial_pertrace_seed == SpikeGenStruct.Ntrial_pertrace)
        % everything ok, do nothing
        seed = SpikeGenStruct.seed;
    elseif ~(Ntrace_seed == Ntrace) 
        disp('Number of seeds given does not match the number of trials. Please give a (Ntrace x Ntrial_pertrace) array (or nothing) as seed for spike train construction')
        disp('Ignoring seed')
        seed = [];
    elseif (Ntrace_seed == Ntrace) && ~(Ntrial_pertrace_seed == SpikeGenStruct.Ntrial_pertrace)
        % This makes no sense
        disp('Number of seeds per recording is not equal to Ntrial_pertrace. Please give a (Ntrace x Ntrial_pertrace) array (or nothing) as seed for spike train construction')
        disp('Ignoring seed')
        seed = [];
    end
else
    seed = [];
end

SpikeTrainStruct.SpikeTimes = cell(Nkernel,Ntrace);
SpikeTrainStruct.SpikeCount = cell(Nkernel,Ntrace);
for nt = 1:Ntrace     
    if ~isempty(seed)
        seednow = seed(nt,:);
    else
        seednow = [];
    end
    for nk = 1:Nkernel
        % creating spike trains 1 by 1
        disp(['Generating ',num2str(SpikeGenStruct.Ntrial_pertrace),' spike train(s) for neuron ',num2str(nk),'/',num2str(Nkernel),' for trace ',num2str(nt),'/',num2str(Ntrace) ])
        [spt, sc] = PSTH_to_SpikeTrain(SpikeTrainStruct.PSTH{nk,nt}, WhiskerStruct.binsize, SpikeGenStruct.refra, SpikeGenStruct.Ntrial_pertrace, SpikeGenStruct.binsize, SpikeGenStruct.delay, SpikeGenStruct.scaling, seednow);
        SpikeTrainStruct.SpikeTimes{nk,nt} = spt;
        SpikeTrainStruct.SpikeCount{nk,nt} = sc;
    end 
end

%% Save
if ~isempty(savename)
    save(savename, 'SpikeTrainStruct');
end

%% plot to check
maxct = 0;
minct = 0;
colorvec = {'b','r'};
for nk = 1:Nkernel
    for nt = 1:Ntrace
        maxt = max(max(ConvTrace{nk,nt}));
        mint = min(min(ConvTrace{nk,nt}));
        if maxt>maxct
            maxct = maxt;
        end
        if mint<minct
            minct = mint;
        end
    end
end
if maxct<1
    maxct = (10^floor(log10(maxct)))*round(maxct/(10^floor(log10(maxct))));
    minct = -maxct;
else
    maxct = ceil(maxct);
    minct = floor(minct);
end
kk = minct:(maxct-minct)/100:maxct;
nkk = length(kk);

if plotyn
    % Kernels
    figure
    plot_thalamic_kernels(KernelStruct)

    % Activation functions
    figure
    nsubplotx = floor(sqrt(Nkernel));
    nsubploty = ceil(sqrt(Nkernel));
    for nk = 1:Nkernel
        subplot(nsubplotx, nsubploty, nk)
        hold all
        if Ndimk == 1
            plot(kk,  feval(KernelStruct.ActivationFunction.function, kk, KernelStruct.ActivationFunction.Params{nk}))
            xlabel('Recording*kernel')
            ylabel(['Spikes / sec'])
            xlim([kk(1), kk(end)])
        elseif Ndimk == 2
            [X,Y] = meshgrid(kk,kk);
            Z = nan*ones(nkk,nkk);
            for nnkk = 1:nkk
                Z(nnkk,:) = feval(KernelStruct.ActivationFunction.function,[X(nnkk,:);Y(nnkk,:)] , KernelStruct.ActivationFunction.Params{nk});
            end
            surf(X,Y, Z, 'LineStyle','none')
            xlabel('Recording*kernel, dimension 1')
            ylabel('Recording*kernel, dimension 2')
            xlim([minct, maxct])
            ylim([minct, maxct])
            zlabel('Spikes / sec')
            view(gca, [-37.5 30]);

        else
            disp('Cannot display activation functions: too many dimensions')
        end
        title(['Kernel ' num2str(nk)])
        box on
        grid on
    end   
    
    % Recordings, PSTHs, Spike Trains
    maxr = 0;
    maxt = 0;
    for nt = 1:Ntrace
        figure
        for nk = 1:Nkernel
            for ndim = 1:Ndimk
                subplot(4,Ndimw,ndim)
                hold all
                if length(WhiskerStruct.Recording{ndim,nt})*WhiskerStruct.binsize>maxt
                    maxt = length(WhiskerStruct.Recording{ndim,nt});
                end
                plot((1:length(WhiskerStruct.Recording{ndim,nt}))*WhiskerStruct.binsize, WhiskerStruct.Recording{ndim,nt})
                title(['Recording dimension ' num2str(ndim)])
                xlim([0 maxt])
                box on
                grid on
                
                subplot(4,Ndimw,Ndimw+ndim)
                hold all
                plot(binsize_kernels*(1:length(ConvTrace{nk,nt}(ndim,:))),ConvTrace{nk,nt}(ndim,:))
                title('Recording * kernel')
                xlim([0 maxt])
                box on 
                grid on                
            end
            
            subplot(4,1,3)
            hold all            
            plot((1:length(SpikeTrainStruct.PSTH{nk, nt}))*binsize_kernels, SpikeTrainStruct.PSTH{nk, nt})
            if max(SpikeTrainStruct.PSTH{nk, nt})>maxr
                maxr = max(SpikeTrainStruct.PSTH{nk, nt});
            end
            title('PSTH')
            ylabel(['Spikes / sec'])
            try
                ylim([0 ceil(maxr)])
                xlim([0 maxt])
            end
            box on
            grid on
            
            subplot(4,1,4)
            hold all
            plot(SpikeTrainStruct.SpikeTimes{nk, nt}, (nk)*ones(size(SpikeTrainStruct.SpikeTimes{nk, nt})), '.')
            xlim([0 maxt])
            ylim([0 nk+1])
            title('Spike times')
            xlabel('time (ms)')
            ylabel('Neuron number')
            box on 
            grid on
        end
    end
end

%% Helper functions
function s = convolve_kernel_acausal( signal, filter, la, lc)
% Convolve signal with kernel with causal and acausal part, give back on same time
% scale.

% Function input:
%   signal(vector) = spike train or input signal
%   filter (vector) = filter
%   la = length acausal part filter (tbegin(negative):0-dt)
%   lc = length causal part filter (dt:tmax/dt)
%   NB lc+la+1 = length(filter) for filters with both causal and acausal
%   parts

% Function output: convolved signal s

cf = 0;
ca = 0;
c0 = 0;

[~, Ny] = size(signal);
if Ny == 1
    signal = signal';
end

if ~(la+lc+1==length(filter))
    if la == 0 && length(filter)==lc
        disp('purely causal filter')
        cf = 1;
    elseif lc == 0 && length(filter)==la
        disp('purely acausal filter')
        ca = 1;
    elseif (la+lc == length(filter))
        disp('filter does not contain t=0 value, making small delay')
        c0 = 1;
    else
        error('give correct lengths of causal and acausal part filter')
    end
end

if cf
    s = conv(signal, [0 filter]);
    s = s(1:length(signal));
elseif ca
    s = conv([signal zeros(1,la)], [filter 0], 'valid');
elseif c0
    s = conv([zeros(1,lc) signal zeros(1,la-1)], filter, 'valid');
else
    s = conv([zeros(1,lc) signal zeros(1,la)], filter, 'valid');
end



end

function [spiketimes, spikecount] = PSTH_to_SpikeTrain(psth, binsize_psth, refra, Ntrial, binsize_spikes, delay, scaling, seed)
% generate Ntrial spike trains for a given psth 
% each of the spike trains is one realization of the psth
% refra: absosolute refractory period (in ms) for the neuron, to prevent bursting
% Ntrial: number of different spike trains to generate
% binsize: binsize (in ms) for generated spike timing data

if nargin == 7
    seed = [];
end
% NB psth should have as units # spikes / sec (so independent of binsize)
nrefra = round(refra/binsize_spikes);
ndelay = round(delay/binsize_spikes);
Ntime = length(psth);

%% Resample and scale
pspike = scaling.*psth*binsize_psth/1000;
if binsize_psth == binsize_spikes
    disp('Sampling rate spikes the same as PSTH')    
else
    if binsize_psth > binsize_spikes
        disp('Sampling rate spikes higher than PSTH: upsample psth')
    elseif binsize_psth < binsize_spikes
        disp('Sampling rate spikes lower than PSTH: downsample psth')
    end
    ff1 = 1/binsize_spikes;
    ff2 = 1/binsize_psth;
    while ~((rem(ff1,1)==0) && (rem(ff2,1)==0))
        ff1 = ff1*10;
        ff2 = ff2*10;
    end
    pspike = resample(pspike,ff1 ,ff2);
end
if max(pspike)>1
    disp('Warning: Pspike has values larger than 1; should not be possible! Reduce scaling.')
    disp(['max(Pspike) = ', num2str(max(pspike))])
end
if min(pspike)<0
    disp('Warning: Pspike has values smaller than 0; should not be possible!')
    disp(['min(Pspike) = ', num2str(min(pspike))])
    keyboard
end
    
%% Make spikes
for ntr = 1:Ntrial
    if ~isempty(seed)
        rng(seed(ntr));
    else
        rng('shuffle');
    end
    p = rand(size(pspike));
    spikeids = find(p<pspike);
    Nspikes = length(spikeids);
    spike_true = ones(size(spikeids));
    % check refractory period
    for ns = 2:Nspikes
        % find last spike
        previous_spikes = find(spike_true(1:ns-1) == 1);
        lastspikeid = spikeids(previous_spikes(end));
        if ((spikeids(ns)-lastspikeid)<nrefra) 
            % spike withing refractory period previous spike
            spike_true(ns) = 0;
        end            
    end

    spikeids = spikeids(find(spike_true==1));
    
    % apply delay
    spikeids = spikeids+ndelay;
    spikeids(find(spikeids>Ntime))=[];
end

spiketimes = spikeids*binsize_spikes;
spikecount = length(spikeids);

    
    

end

function plot_thalamic_kernels(KernelStruct)
    [Nbx, Nby] = size(KernelStruct);
    Nbarrel = Nbx*Nby;
    nb = 0;
    for nbx = 1:Nbx
        for nby=1:Nby
            nb = nb+1;
            [Nkernel, Ndimk]  = size(KernelStruct.Kernels);
            for nd = 1:Ndimk
                subplot(Nbarrel,Ndimk, (nb-1)*Ndimk+nd)
                title(['Dimension ', num2str(nd) ', barrel ' num2str(nb)])
                hold all
                for nk = 1:Nkernel            
                    plot(-1*fliplr(KernelStruct.kerneltime), fliplr(KernelStruct.Kernels{nk, nd}))
                    xlabel('time before spike (ms)')
                    ylabel('stimulus amplitude')
                    grid on
                    box on
                end
            end
        end
    end
end

end

