function [f_example_1vol_l23, f_average_1vol_l23,f_example_1vol_l4, f_average_1vol_l4, zrange_vol_l23, zrange_vol_l4] = load_plot_Svoboda_data_to_compare(InputDataStruct, vol_l23, vol_l4, depth_l23, depth_l4, edges, edgesav)

load([InputDataStruct.loadfolder InputDataStruct.animal '_' InputDataStruct.session{1} '_' InputDataStruct.dataname])
nvol = length(s.timeSeriesArrayHash.value); % # volumes in data set

disp('Loading Svoboda data')
M = get_first_touch(s, 0); % matrix with 1) trial numbers, 2) times (whisker base), 3) protraction (1) or retraction (2), 4) volume nr of first touches and 5) times (lum base)
LMat = cell(1,nvol-1);
luminesmat = cell(1,nvol);
whiskermat = cell(1,nvol);
alltrials = cell(1,nvol);
validtrials = cell(1,nvol); % 0/1 trials that have recordings
touchtrials = cell(1,nvol); % trials that have touch data
Ntrial = nan(1,nvol);
Ncell = nan(1,nvol);
for nv = 2:nvol
    % NB first 'volume' is whisker data

    % Cell locations
    LMat{nv} = find_xyz_roi( s, nv );

    % trials
    alltrials{nv} = unique(s.timeSeriesArrayHash.value{nv}.trial);
    [whiskermat{nv}, luminesmat{nv}, ~, ~, ~, validtrials{nv}] = load_data_across_sessions(InputDataStruct.loadfolder, ...
        InputDataStruct.animal, InputDataStruct.session, InputDataStruct.dataname, nv, InputDataStruct.window, 0);
    [Ncell(nv),~,Ntrial(nv) ] = size(luminesmat{nv-1});
    relmn = find(M(4,:) == nv); % relevant columns in M
    touchtrials{nv} = [relmn; M(1,relmn)]; 
end

zrange_vol_l23 = [min(LMat{vol_l23}(4,:)) max(LMat{vol_l23}(4,:))];
zrange_vol_l4 = [min(LMat{vol_l4}(4,:)) max(LMat{vol_l4}(4,:))];

% Check if trials have both touch data and recordings
% NB Note that load_data_across_sessions (whiskermat, luminesmat, validtrials) gives all trials with recordings,
% % whereas get_first_touch (M) gives all trials with touch data. 
trials_recording_touch = cell(1,nvol);
for nv = 2:nvol
    trials_now = alltrials{nv};
    if length(trials_now) == length(validtrials{nv})
        Nkepttrials = sum(validtrials{nv}(2,:));
        keepvec = nan(1,Nkepttrials);
        keepvec_counter = 0;
        for nt = 1:length(trials_now)
            trial_id = trials_now(nt);
            if (validtrials{nv}(1,nt)==trial_id && validtrials{nv}(2,nt)==1) 
                % trial has recording
                keepvec_counter = keepvec_counter+1;
                if sum(touchtrials{nv}(2,:) == trial_id)==1
                    % trial has touch data 
                    trials_recording_touch{nv} = [trials_recording_touch{nv} trial_id];
                    keepvec(keepvec_counter) = 1;
                else
                    % trial has no touch data
                    keepvec(keepvec_counter) = 0;
                    whiskermat{nv} = whiskermat{nv}(:,:, keepvec);
                    luminesmat{nv} = luminesmat{nv}(:,:,keepvec);
                end
            end
        end
    else
        disp('Check valid trials')
        keyboard
    end
    if sum(keepvec) == length(keepvec)
        disp(['All recorded trials have touch data for volume ' num2str(nv)])
    else
        keyboard        
    end
end


%% Get&plot relevant neuron trials and neurons for each layer
disp('Plotting Svoboda data')
% L23:  
% depth
[f_example_l23, f_average_l23, volumes_l23, trials_l23] = plot_activity_depth(InputDataStruct, depth_l23, LMat, s,edges, 0);
savefig(f_example_l23,[InputDataStruct.savefolder f_example_l23.Name '.fig'],'compact')
savefig(f_average_l23,[InputDataStruct.savefolder 'average_l23.fig'],'compact')
% single volume
InputDataStruct.volume = vol_l23;
[f_example_1vol_l23, f_average_1vol_l23, trials_1vol_l23] = plot_activity_single_volume(InputDataStruct, LMat, s, edges, edgesav, 0);
savefig(f_example_1vol_l23,[InputDataStruct.savefolder f_example_1vol_l23.Name '.fig'],'compact')
savefig(f_average_1vol_l23,[InputDataStruct.savefolder 'average_v' num2str(vol_l23) '.fig'],'compact')

% L4: 
% depth
[f_example_l4, f_average_l4, volumes_l4, trials_l4] = plot_activity_depth(InputDataStruct, depth_l4, LMat, s,edges, 0);
savefig(f_example_l4,[InputDataStruct.savefolder f_example_l4.Name '.fig'],'compact')
savefig(f_average_l4,[InputDataStruct.savefolder 'average_l4.fig'],'compact')
% single volume
InputDataStruct.volume = vol_l4;
[f_example_1vol_l4, f_average_1vol_l4, trials_1vol_l4] = plot_activity_single_volume(InputDataStruct, LMat, s, edges, edgesav, 0);
savefig(f_example_1vol_l4,[InputDataStruct.savefolder f_example_1vol_l4.Name '.fig'],'compact')
savefig(f_average_1vol_l4,[InputDataStruct.savefolder 'average_v' num2str(vol_l4) '.fig'],'compact')