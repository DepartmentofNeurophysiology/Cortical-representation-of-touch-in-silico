function [whiskermat, luminesmat, upsample_rate, binsize_whisker, dtstart, validtrials] = load_data_across_sessions(loadfolder, animal, sessionvec, dataname, volume, window, plotcheck)
% Load data from the same animal, across multiple sessions
% INPUT:
% * loadfolder (string): folder name of where the session data files are
% * animal (string): animal name
% * sessionvec (string array): array of sessions to load
% * dataname (string): name of session data structs
% * volume (int): which volume to load 
% * window (struct, optional) with fields
% *     .start: 'pole in reach' or 'first touch' or 'first' 
% *     .window (2x1 array with window (ms) around start, first touch or pole in reach trial
% * plotcheck (logical, optional): whether to pause&plot every trial for visual inspection

% OUTPUT
% * whiskermat: 2xNtimexNtrial double array with whisker recordings
% * luminesmat: NneuronxNtimexNtrial array with calcium recordings
% * upsample_rate: sampling rate whisker recordings / sampling rate calcium recordings
% * binsize_whisker (ms)
% * dtstart: starttime_whisker-starttime_luminescence;

% NB 
% * Ntime is not the same for whisker and calcium recordings, due to different sampling rates
% * Only neurons that are recorded across all sessions are included
% * Only trials with whisker recordings are included, and where not all calcium recordings are 'NaN' 


if nargin == 5
    window.window = [];
    window.start = [];
    plotcheck = 0;
elseif nargin == 6
    if isnumeric(window)
        window.window = [];
        window.start = [];
        plotcheck = window;
    else
        plotcheck = 0;
    end
end

if (strcmp(window.start, 'first') && window.window(1)<0)
    disp('Cannot include data before the start of the trial. Starting at 0.')
    window.window(1) = 0;
end

Nsession = length(sessionvec);
validsessionvec = ones(1,Nsession);
Nneuron_max = 0;
Ntrial_tot = 0;     % total number of trials
Ntrial_totkeep = 0; % total number of kept trials
nllmax = 0;         % maximum number of data points luminescence recordings
nwlmax = 0;         % maximum number of data points whisker recordings
nantrials_tot = 0;
notouchtrials_tot = 0;
nowhiskertrials_tot = 0;
for ns = 1:Nsession
    %% Load 
    filename = [loadfolder animal '_' sessionvec{ns} '_' dataname '.mat'];
    disp(['Loading file ' filename])
    load(filename) % any file from https://crcns.org/data-sets/ssc/ssc-2
    nVolume = length(s.timeSeriesArrayHash.value);
    if volume>nVolume
        % volume not present in this session
        disp(['Volume ' num2str(volume) ' does not exist for session ' sessionvec{ns} '; skipping session.'])
        validsessionvec(ns) = 0;
        continue
    else    
        %% Checks
        % check if multiple whiskers were recorded
        [nwhiskertrace, ~] = size(s.timeSeriesArrayHash.value{1}.valueMatrix);
        if nwhiskertrace>2
            disp('More than 2 whisker traces detected:')
            for nw = 1:nwhiskertrace
                disp(s.timeSeriesArrayHash.value{1}.idStrs{nw})
            end
            whiskertracen = input('Which traces should be used? (give 2d array)');
        else
            whiskertracen = [1,2];
        end

        % number of trials for this volume
        trialvec = unique(s.timeSeriesArrayHash.value{volume}.trial);
        Ntrial_temp = length(trialvec);

        binsize_lum = s.timeSeriesArrayHash.value{volume}.time(2)-s.timeSeriesArrayHash.value{volume}.time(1);
        binsize_whisker = s.timeSeriesArrayHash.value{1}.time(2)-s.timeSeriesArrayHash.value{1}.time(1);
        upsample_rate_new = ((binsize_lum)/(binsize_whisker));

        if exist('upsample_rate_old','var')
            if ~(upsample_rate_new == upsample_rate_old)
                error('Files do not have the same sampling rates')
            end
        end

        % Check if the number of recorded neurons corresponds to other sessions
        [Nneuron_new, ~] = size(s.timeSeriesArrayHash.value{volume}.valueMatrix);
        if Nneuron_new>Nneuron_max
            Nneuron_max = Nneuron_new;
        end
        ids_new = s.timeSeriesArrayHash.value{volume}.ids;
        if exist('Nneuron_old', 'var')
            if ~(Nneuron_new == Nneuron_old)
                disp('Not the same number of neurons in new file. Checking ids')
                neuronkeepvec_new = [];
                neuronkeepvec_old = ones(1,Nneuron_old);
                for nn = 1:Nneuron_old
                    corresponding_neuron = find(ids_new == ids_old(nn));
                    if isempty(corresponding_neuron)
                        neuronkeepvec_old(nn) = 0;
                    else
                        neuronkeepvec_new = [neuronkeepvec_new; corresponding_neuron];
                    end
                    corresponding_neuron = [];
                end
                if isempty(neuronkeepvec_new)
                    disp('No neurons left, skipping this session')
                    validsessionvec(ns) = 0;
                    continue
                end
                neuronkeepvec_old = find(neuronkeepvec_old==1); % for later use at concatenation

                Nneuron_new = length(neuronkeepvec_new);
                ids_new = ids_new(neuronkeepvec_new);
                if Nneuron_new<Nneuron_old
                    keep = input([num2str(Nneuron_old-Nneuron_new) ' neurons discarded; ' num2str(Nneuron_new) ' neurons kept. Keep session? (y/n)'],'s');
                    if strcmp(keep, 'n')
                        disp(['Discarding session ' sessionvec{nsession}])
                        validsessionvec(nsession) = 0;
                        continue
                    end
                end

            else
                disp('The same number of neurons in new file. Not checking ids')
                neuronkeepvec_new = 1:Nneuron_new;
                neuronkeepvec_old = neuronkeepvec_new;
            end
        else
            disp('First file. Not checking ids')
            neuronkeepvec_new = 1:Nneuron_new;
            neuronkeepvec_old = neuronkeepvec_new;
        end


        %% Put (valid) trials in matrix

        if strcmp(window.start, 'pole in reach')
            Ntime_w = 1000;
            Ntime_l = 100;
        elseif (~isempty(window.window) && ~strcmp(window.start, 'pole in reach'))
            whiskertime = window.window(1):binsize_whisker:window.window(end);
            Ntime_w = length(whiskertime);

            lumtime = window.window(1):binsize_lum:window.window(end);
            Ntime_l = length(lumtime);
        else
            error('Please give an appropriate time window')
        end

        whiskermat_temp = NaN*ones(2,Ntime_w,Ntrial_temp);
        luminesmat_temp = NaN*ones(Nneuron_new,Ntime_l,Ntrial_temp);

        validvec = ones(size(trialvec));
        dtstart_temp = nan*ones(1,Ntrial_temp);

        disp(['Number of trials for volume ' num2str(volume) ': ' num2str(Ntrial_temp)])

        nantrials_temp = 0;
        notouchtrials_temp = 0;
        nowhiskertrials_temp = 0;

        for nt = 1:Ntrial_temp
            disp(['Trial ' num2str(nt) '/' num2str(Ntrial_temp) ': trial id ' num2str(trialvec(nt))])
            % nt_tot = nt+Ntrial_tot;
            whiskertrialvec = find(s.timeSeriesArrayHash.value{1}.trial == trialvec(nt));

            if isempty(whiskertrialvec)
                disp(['No whisker data, skipping trial with id ' num2str(trialvec(nt))])
                validvec(nt)=0;
                nowhiskertrials_temp = nowhiskertrials_temp+1;

            else
                luminestrialvec = find(s.timeSeriesArrayHash.value{volume}.trial == trialvec(nt));

                %% Make appropriate time window
                whiskertime_thistrial = s.timeSeriesArrayHash.value{1}.time(whiskertrialvec);
                luminestime_thistrial = s.timeSeriesArrayHash.value{volume}.time(luminestrialvec);
                
                whiskertrace = s.timeSeriesArrayHash.value{1}.valueMatrix(whiskertracen,whiskertrialvec);

                luminestrace = s.timeSeriesArrayHash.value{volume}.valueMatrix(neuronkeepvec_new,luminestrialvec);
                tstart = max(whiskertime_thistrial(1), luminestime_thistrial(1));
                tend = min(whiskertime_thistrial(end), luminestime_thistrial(end));
                if isempty(window.start)
                   % nothing needed, use calculated tstart and tend;
                elseif strcmp(window.start, 'pole in reach')
                    if isfield(window, 'window')
                        disp('Using times pole in reach; ignoring given window')
                    end
                    tpole  = s.eventSeriesArrayHash.value{1}.eventTimes(s.eventSeriesArrayHash.value{1}.eventTrials==trialvec(nt));
                    tstart = tpole(1);
                    tend   = tpole(2);
                elseif strcmp(window.start, 'first touch')
                    ttouchpro = s.eventSeriesArrayHash.value{2}.eventTimes{1}(s.eventSeriesArrayHash.value{2}.eventTrials{1}==trialvec(nt));
                    ttouchre  = s.eventSeriesArrayHash.value{2}.eventTimes{2}(s.eventSeriesArrayHash.value{2}.eventTrials{2}==trialvec(nt));
                    if ~isempty(ttouchpro) && ~isempty(ttouchre)
                        ttouchpro = ttouchpro(1);
                        ttouchre  = ttouchre(1);
                        ttouch = min(ttouchpro, ttouchre);
                    elseif ~isempty(ttouchpro) && isempty(ttouchre)
                        ttouch = ttouchpro(1);
                    elseif isempty(ttouchpro) && ~isempty(ttouchre)
                        ttouch = ttouchre(1);
                    else
                        disp(['No touch data, skipping trial with id ' num2str(trialvec(nt))])
                        validvec(nt)=0;
                        ttouch = tstart;
                        notouchtrials_temp = notouchtrials_temp+1;                    
                    end                    
                    tstart = ttouch+window.window(1);
                    tend   = ttouch+window.window(2);
                elseif strcmp(window.start, 'first')
                    tstart = tstart+window.window(1);
                    tend   = tstart+window.window(2);
                end


                %% Find appropriate whisker recordings
                if strcmp(window.start, 'pole in reach')
                    % variable length
                    [~, nstart] = min(abs(whiskertime_thistrial-tstart));
                    [~, nend] = min(abs(whiskertime_thistrial-tend));
                    whiskertrace = whiskertrace(:,nstart:nend);
                else
                    % fixed length
                    if ((tstart>whiskertime_thistrial(1)) && (tend<whiskertime_thistrial(end)))
                        % window fits in trial
                        [~, nstart] = min(abs(whiskertime_thistrial-tstart));
                        whiskertrace = whiskertrace(:,nstart:nstart+Ntime_w-1);
                        whiskertime_thistrial = whiskertime_thistrial(nstart:nstart+Ntime_w-1);
                    elseif (tstart<whiskertime_thistrial(1) && ~strcmp(window.start, 'pole in reach'))
                        % no recording in beginning of window, add NaN
                        [~, nend] = min(abs(whiskertime_thistrial-tend));
                        whiskertrace_temp = whiskertrace(:,1:nend);
                        whiskertime_thistrial = whiskertime_thistrial(1:nend);
                        lw = length(whiskertrace_temp(1,:));
                        whiskertrace = [nan*ones(2,Ntime_w-lw), whiskertrace_temp];
                    elseif (tend>whiskertime_thistrial(end) && ~strcmp(window.start, 'pole in reach'))
                        % no recording in end of window, add NaN at the end
                        [~, nstart] = min(abs(whiskertime_thistrial-tstart));
                        whiskertrace_temp = whiskertrace(:,nstart:end);
                        whiskertime_thistrial = whiskertime_thistrial(nstart:end);
                        lw = length(whiskertrace_temp(1,:));
                        whiskertrace = [whiskertrace_temp, nan*ones(2,Ntime_w-lw)];
                    else
                        error('Chosen window too large for trials, please chose a smaller window')
                    end
                end
                if strcmp(window.start, 'pole in reach')
                    nwl = length(whiskertrace(1,:));
                    if nwl>nwlmax
                        nwlmax = nwl;
                    end
                    whiskermat_temp(:,1:nwl,nt) = whiskertrace;
                else
                    try
                        whiskermat_temp(:,:,nt) = whiskertrace;
                    catch
                        keyboard
                    end
                end



                %% Find appropriate neural recordings
                if strcmp(window.start, 'pole in reach')
                    % variable length
                    [~, nstart] = min(abs(luminestime_thistrial-tstart));
                    [~, nend] = min(abs(luminestime_thistrial-tend));
                    luminestrace = luminestrace(:,nstart:nend);
                else
                    % fixed length
                    if ((tstart>=luminestime_thistrial(1)) && (tend<=luminestime_thistrial(end)))
                        % window fits in trial   
                        [~, nstart] = min(abs(luminestime_thistrial-tstart));
                        if ((luminestime_thistrial(nstart)<tstart) && ((nstart+Ntime_l)<=length(luminestime_thistrial)))
                            % causality: always align to next neural recording
                            nstart = nstart+1;
                        end
                        luminestrace = luminestrace(:,nstart:nstart+Ntime_l-1);
                        luminestime_thistrial = luminestime_thistrial(nstart:nstart+Ntime_l-1);
                    elseif (tstart<luminestime_thistrial(1))
                        % no recording in beginning of window, add NaN
                        [~, nend] = min(abs(luminestime_thistrial-tend));
                        if ((luminestime_thistrial(nend)<tend) && (nend+1<= (length(luminestime_thistrial))))
                            % causality: always align to next neural recording
                            nend = nend+1;
                        end
                        if nend<Ntime_l
                            luminestrace_temp = luminestrace(:,1:nend); 
                            luminestime_thistrial = luminestime_thistrial(1:nend);
                            lw = length(luminestrace_temp(1,:));
                            luminestrace = [nan*ones(Nneuron_new,Ntime_l-lw), luminestrace_temp];
                        else
                            luminestrace = luminestrace(:,nend-Ntime_l+1:nend);
                        end
                    elseif (tend>luminestime_thistrial(end))
                        % no recording in end of window, add NaN at the end
                        [~, nstart] = min(abs(luminestime_thistrial-tstart));
                        if ((luminestime_thistrial(nstart)<tstart) && ((nstart+Ntime_l)<=length(luminestime_thistrial)))
                            % causality: always align to next neural recording
                            nstart = nstart+1;
                        end
                        luminestrace_temp = luminestrace(:,nstart:end);
                        luminestime_thistrial = luminestime_thistrial(nstart:end);
                        lw = length(luminestrace_temp(1,:));
                        luminestrace = [luminestrace_temp, nan*ones(Nneuron_new,Ntime_l-lw)];
                    else
                        error('Chosen window too large for trials, please chose a smaller window')
                    end
                end
                if strcmp(window.start, 'pole in reach')
                    nll = length(luminestrace(1,:));
                    if nll>nllmax
                        nllmax = nll;
                    end
                    luminesmat_temp(:,1:nll,nt) = luminestrace;
                else
                    try
                        luminesmat_temp(:,:,nt) = luminestrace;
                    catch
                        keyboard
                    end
                end

                dtstart_temp(nt) = whiskertime_thistrial(1)-luminestime_thistrial(1);

            end

            if Nneuron_new*Ntime_l == sum(sum(isnan(luminesmat_temp(:,:,nt))))
                disp(['Only NaNs in neural recording, skipping trial with id ' num2str(trialvec(nt))])
                validvec(nt)=0;
                nantrials_temp = nantrials_temp + 1;

            end
            if plotcheck
                f = figure;
                subplot(3,1,1)
                plot(whiskertime, squeeze(whiskermat_temp(1,:,nt)))
                box on
                grid on
                title('Whisker angle')
                subplot(3,1,2)
                plot(whiskertime, squeeze(whiskermat_temp(2,:,nt)))
                box on
                grid on
                title('Whisker curvature')
                subplot(3,1,3)
                plot(lumtime, squeeze(luminesmat_temp(:,:,nt)))
                box on
                grid on
                title('\Delta F / F')
                xlabel('time (ms)')
                pause
                close(f)
            end
        end
        Ntrial_tempkeep = sum(validvec);
        whiskermat_temp = whiskermat_temp(:,:,validvec == 1);
        luminesmat_temp = luminesmat_temp(:,:,validvec == 1);
        dtstart_temp    = dtstart_temp(validvec == 1);
        validtrials_temp = [trialvec; validvec];

        %% Concatenate with previous sessions
        if exist('whiskermat','var')
            if (~isempty(luminesmat_temp)) && (~isempty(whiskermat_temp))
                whiskermat = cat(3,whiskermat, whiskermat_temp);
            end
        else
            whiskermat = whiskermat_temp;
        end
        if exist('luminesmat','var')
            if (~isempty(luminesmat_temp)) && (~isempty(whiskermat_temp))
                luminesmat = cat(3,luminesmat(neuronkeepvec_old,:,:), luminesmat_temp);
            end
        else
            luminesmat = luminesmat_temp;
        end

        if exist('dtstart','var')
            if (~isempty(luminesmat_temp)) && (~isempty(whiskermat_temp))
                dtstart = [dtstart dtstart_temp];
            end
        else
            dtstart = dtstart_temp;
        end
        
        if exist('validtrials','var')
            if (~isempty(luminesmat_temp)) && (~isempty(whiskermat_temp))
                validtrials = [validtrials validtrials_temp];
            end
        else
            validtrials = validtrials_temp;
        end

        %% For the next round
        Nneuron_old = Nneuron_new;
        ids_old = ids_new;
        upsample_rate_old = upsample_rate_new;
        if ~(length(ids_old)==Nneuron_old)
            error('Number of neuron ids does not correspond to number of recordings')
        end
        Ntrial_tot = Ntrial_tot+Ntrial_temp;
        Ntrial_totkeep = Ntrial_totkeep+Ntrial_tempkeep;
        nantrials_tot = nantrials_tot+nantrials_temp;
        notouchtrials_tot = notouchtrials_tot+notouchtrials_temp;
        nowhiskertrials_tot = nowhiskertrials_tot+nowhiskertrials_temp;
    end
end

%% Discard padded zeros
if strcmp(window.start, 'pole in reach')
    whiskermat = whiskermat(:,1:nwlmax,:);
    luminesmat = luminesmat(:,1:nllmax,:);
end

upsample_rate = upsample_rate_new;

%% Print summary
disp(['In a total of ' num2str(Nsession) ' sessions, ' num2str(sum(validsessionvec)) ' were valid (overlapping neuron ids with first session).'])
disp(['There were ' num2str(Ntrial_tot) ' trials, of which ' num2str(Ntrial_totkeep) ' were kept.'])
disp(['There were ' num2str(nantrials_tot) ' trials with only NaNs in the recordings'])
disp(['There were ' num2str(notouchtrials_tot) ' trials with no first touch data'])
disp(['There were ' num2str(nowhiskertrials_tot) ' trials with no whisker data'])
disp('NB Overlap possible')
disp(['Out of a maximum of ' num2str(Nneuron_max) ' neurons, ' num2str(Nneuron_old) ' were kept (shared across all sessions)'])