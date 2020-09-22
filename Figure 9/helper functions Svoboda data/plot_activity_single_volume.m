function [f_trial1, f_average, reltrials] = plot_activity_single_volume(InputDataStruct, LMat, s, edges, edgesav, showyn)


%% find relevant trials and get data
% make normalized (cells) histogram over cells for each trial, and take
% average



totcounts = zeros(1,length(edges)-1);
totcountsdiff= zeros(1,length(edges)-1);
savetrial = 1;
counter = 0;
[~, Ncell] = size(LMat{InputDataStruct.volume});

[SampleMatlum, reltrials, time] = get_plot_window_lum_first_touch(s, InputDataStruct.window.window, InputDataStruct.volume, 0);
timebin0 = find(time == 0);
[~, totaltrials] = size(reltrials);
differencemat = nan(totaltrials, Ncell);
activitymat = nan(totaltrials, Ncell);
meanlummat = nan(totaltrials, Ncell);
stdlummat = nan(totaltrials, Ncell);
for nt = 1:totaltrials
    counter = counter+1;
    if counter == savetrial
        f_trial1 = figure('Name',['Volume ' num2str(InputDataStruct.volume) ' trial id ' num2str(reltrials(2,nt))]);
        showyn_now = 1;
    else
        showyn_now = showyn;
        f = figure('Name',['Volume ' num2str(InputDataStruct.volume) ' trial id ' num2str(reltrials(2,nt))]);
    end
    [activitymat(nt,:), differencemat(nt,:), counts, countsdiff] = plot_single_trial(LMat{InputDataStruct.volume}(2,:), LMat{InputDataStruct.volume}(3,:), squeeze(SampleMatlum(nt,:,:)), timebin0, edges);
    [meanlummat(nt,:), stdlummat(nt,:)] =  get_mean_std_lum_per_neuron(squeeze(SampleMatlum(nt,:,:)), 0);
    totcounts = totcounts+counts;
    totcountsdiff = totcountsdiff+countsdiff;

    
    if showyn_now
        pause
    end
    if counter == savetrial

    else
        close(f)
    end

end

%% Plot averages
% NB Note that averaging over histograms is not the same as averaging over
% neurons and then making a histogram!
f_average = figure('Name',['Volume ' num2str(InputDataStruct.volume) ' average' ]);
plot_average(LMat{InputDataStruct.volume}(2,:), LMat{InputDataStruct.volume}(3,:),activitymat, differencemat, edgesav)

