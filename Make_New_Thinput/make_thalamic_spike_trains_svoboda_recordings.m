function SpikeTrainStruct = make_thalamic_spike_trains_svoboda_recordings(savefolder, savename_input, SvobodaStruct, Barrelstruct, make_new_thalamic_kernels)
% Make Thalamic spike trains from Svoboda recording.

% Input:
% * savefolder & savename: folder where to store file with spike trains
% * SvobodaStruct with fields
%       * loadfolder (string): folder with Svoboda data files
%       * dataname (string): file name suffix
%       * volume (integer > 2): trials from which volume to use
%       * animal (string): animal name
%       * sessionvec (string array): array of sessions to load
%       * window (struct) with fields
%           *.start: 'pole in reach' or 'first touch' or 'first'
%           *.window (2x1 array with window (ms) around start, first touch or pole in reach trial
%       * trialvec (optional): vector of which trials to use
% * BarrelStruct{barrelidX, barrelIDY} with fields
%       * .mainbarrel (1,2,3): main, secondary or tertiary barrel
% * make_new_thalamic_kernels: 0 (use existing file) or 1 (make new kernels)

disp('Making new Thalamic spike trains')
    
f = filesep;
    %% Load Svoboda recordings
    addpath(genpath(SvobodaStruct.loadfolder));

    plotcheck = 0;
    [whiskermat, ~, ~, binsize_whisker] = load_data_across_sessions(SvobodaStruct.loadfolder, SvobodaStruct.animal, SvobodaStruct.sessionvec, SvobodaStruct.dataname, SvobodaStruct.volume, SvobodaStruct.window, plotcheck);
    [~, ~, Ntrialtot] = size(whiskermat);
    if isfield(SvobodaStruct, 'trialvec')
        if ~isempty(SvobodaStruct.trialvec)
            trialvec = SvobodaStruct.trialvec;
        else
            trialvec = 1:Ntrialtot;
        end
    else
      trialvec = 1:Ntrialtot;
    end
    Ntrial = length(trialvec);
    for nt = 1:Ntrial
        WhiskerTrace.Recording{1,nt} = pi*squeeze(whiskermat(1,:,trialvec(nt)))/180;    % from degree to radian (base angle)
        WhiskerTrace.Recording{2,nt} = squeeze(whiskermat(2,:,trialvec(nt)));           % curvature
    end
    WhiskerTrace.binsize = binsize_whisker; % ms
    WhiskerTrace.quantity = {'base_angle','curvature'};                                 % See Whisker_Recording class for recognized quantities and units
    WhiskerTrace.unit = {'radian','mm-1'};
    
    %% Make / load kernels
    if make_new_thalamic_kernels
        disp('Making new Thalamic kernels')
        % Make new kernels
        thalamickernelfolder = ['..' f 'Make_New_Thinput' f];
        addpath(genpath(thalamickernelfolder));
        
        binsize = 1; % ms
        [Nbx, Nby] = size(Barrelstruct);
        Nbarrel = Nbx*Nby;
        KernelStruct = cell(Nbx, Nby);
        for nbx = 1:Nbx
            for nby = 1:Nby
                KernelStruct{nbx,nby} = make_kernels_angle_curvature_Svobodaexp(SvobodaStruct.Nkernel_ba, SvobodaStruct.Nkernel_c, SvobodaStruct.Nkernel_m, binsize, [savefolder SvobodaStruct.savename '_Thalamic_Kernels_barrel_' num2str(nbx) '_' num2str(nby)]);
            end
        end
        save([savefolder SvobodaStruct.savename '_Thalamic_Kernels'], 'KernelStruct')
        % delete temp files
        for nbx = 1:Nbx
            for nby = 1:Nby 
                delete([savefolder SvobodaStruct.savename '_Thalamic_Kernels_barrel_' num2str(nbx) '_' num2str(nby) '.mat'])
            end
        end   
    else
        % Load kernels from file
         load([savefolder SvobodaStruct.savename '_Thalamic_Kernels']);
         [Nbx, Nby] = size(KernelStruct);
         Nbarrel = Nbx*Nby;
    end
    
    %% Make spike trains for each barreloid
    SpikeTrainStruct = cell(Nbx,Nby);
    SpikeGenStruct = cell(Nbx,Nby);
    nb = 0;
    for nbx = 1:Nbx
        for nby = 1:Nby    
            nb = nb+1;
            disp(['Making Thalamic spike trains for barrel ' num2str(nb) '/' num2str(Nbarrel)])
            SpikeGenStruct{nbx,nby}.refra             = 3; % refractory period (ms)
            SpikeGenStruct{nbx,nby}.Ntrial_pertrace   = 1; % # trials for each Deflection trace
            SpikeGenStruct{nbx,nby}.binsize           = 1; % binsize spike trains (ms)
            if Barrelstruct{nbx,nby}.mainbarrel == 1
                % principal barreloid
                SpikeGenStruct{nbx,nby}.delay             = 0; % (ms)
                SpikeGenStruct{nbx,nby}.scaling           = 1; % Scaling of PSTH
            elseif Barrelstruct{nbx,nby}.mainbarrel == 2
                % secondary barreloid
                SpikeGenStruct{nbx,nby}.delay             = 2.5; % (ms)
                SpikeGenStruct{nbx,nby}.scaling           = .3;
            elseif Barrelstruct{nbx,nby}.mainbarrel == 3
                % tertiary barreloid: no spikes
                SpikeGenStruct{nbx,nby}.delay             = 0; % (ms)
                SpikeGenStruct{nbx,nby}.scaling           = 0;
            end
            % Generate spike trains
            plotyn = 1;
            SpikeTrainStruct{nbx,nby} = kernel_recording_to_spiketrain(WhiskerTrace, KernelStruct{nbx,nby}, SpikeGenStruct{nbx,nby}, [savefolder savename_input '_Thalamic_Spike_Trains_barrel_' num2str(nbx) '_' num2str(nby)], plotyn);
        end
    end
    save([savefolder savename_input '_Thalamic_Spike_Trains'], 'SpikeTrainStruct', 'WhiskerTrace', 'SvobodaStruct')
    
    %% delete temp files
    for nbx = 1:Nbx
        for nby = 1:Nby 
            delete([savefolder savename_input '_Thalamic_Spike_Trains_barrel_' num2str(nbx) '_' num2str(nby) '.mat'])
        end
    end    
end