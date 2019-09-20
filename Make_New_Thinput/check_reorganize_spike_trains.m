function Input_spike_trains = check_reorganize_spike_trains(SpikeTrainStruct,barrelstruct)

[Nbx, Nby] = size(barrelstruct);
Nbarrel = Nbx*Nby;

% Check number of barrels
[Nbx_st, Nby_st] = size(SpikeTrainStruct);
if (~(Nbx == Nbx_st) || ~(Nby_st==Nby))
    error('Number of barrels in Thalamic spike trains does not fit number of barrels in connectivity data')
end
Nthneuron = 0;
for nbx = 1:Nbx
    for nby = 1:Nby
        [Nneuron, Ntrial] = size(SpikeTrainStruct{nbx,nby}.SpikeTimes);
        if ~(Nneuron == barrelstruct{nbx,nby}.Nthalamic)
            error(['Number of thalamic neurons in barrel ' num2str(nbx) ' , ' num2str(nby) ' in the spike trains does not fit the connectivity data'])
        else
            Nthneuron = Nthneuron + Nneuron;
        end
    end
end


% Reorganize
Input_spike_trains = cell(1,Ntrial);
for nt = 1:Ntrial
    Input_spike_trains{nt} = zeros(Nthneuron, 1000);
    nspikesmax = 0;
    NN = 0;
    for nbarrel = 1:Nbarrel
        [NNeuron, ~] = size(SpikeTrainStruct{nbarrel}.SpikeTimes);
        for nneuron = 1:NNeuron
            nspikes = length(SpikeTrainStruct{nbarrel}.SpikeTimes{nneuron, nt});
            if nspikes>nspikesmax
                nspikesmax = nspikes;
            end
            Input_spike_trains{nt}(nneuron+NN,1:nspikes) = SpikeTrainStruct{nbarrel}.SpikeTimes{nneuron, nt};
        end
        NN = NN+NNeuron;
    end
    Input_spike_trains{nt} = Input_spike_trains{nt}(:,1:nspikesmax);
end
