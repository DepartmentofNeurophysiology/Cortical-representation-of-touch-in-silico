function [f_trial1, f_average, rel_volumes, reltrials] = plot_activity_depth(InputDataStruct, depth, LMat, s, edges, showyn)

bins = edges+(edges(2)-edges(1))/2;
bins = bins(1:end-1);

%% find cell ids
nvol = length(s.timeSeriesArrayHash.value); % # volumes in data set
rel_volumes = [];
rel_cells = [];
for nv = 2:nvol
    relcel_temp = [];
    relcel_temp = LMat{nv}(1,(LMat{nv}(4,:) >=depth(1) & LMat{nv}(4,:) <=depth(2)));
    if ~isempty(relcel_temp)
        rel_volumes = [rel_volumes, nv];
        rel_cells = [rel_cells, relcel_temp];
    end
end

%% find relevant trials
% make normalized (cells) histogram over cells for each trial, and take
% average

SampleMatlum = cell(1,length(rel_volumes));
reltrials = cell(1,length(rel_volumes));
totaltrials = 0;
totcounts = zeros(1,length(edges)-1);
totcountsdiff= zeros(1,length(edges)-1);
savetrial = 1;
counter = 0;
for nrv = 1:length(rel_volumes)
    volume_now = rel_volumes(nrv);
    [SampleMatlum{nrv}, reltrials{nrv}, time] = get_plot_window_lum_first_touch(s, InputDataStruct.window.window, volume_now, 0);
    timebin0 = find(time == 0);
    [~,ntrial_now] = size(reltrials{nrv});
    totaltrials = totaltrials+ntrial_now;
    for nt = 1:ntrial_now
        counter = counter+1;
        if counter == savetrial && nrv == 1
            f_trial1 = figure('Name',['Volume ' num2str(volume_now) ' trial id ' num2str(reltrials{nrv}(2,nt))]);
            showyn_now = 1;
        else
            showyn_now = showyn;
            f = figure('Name',['Volume ' num2str(volume_now) ' trial id ' num2str(reltrials{nrv}(2,nt))]);
        end
        [~, ~, counts, countsdiff] = plot_single_trial(LMat{volume_now}(2,:), LMat{volume_now}(3,:), squeeze(SampleMatlum{nrv}(nt,:,:)), timebin0, edges);
    
        totcounts = totcounts+counts;
        totcountsdiff = totcountsdiff+countsdiff;
        
        if showyn_now
%             pause
        end
        if counter == savetrial
            
        else
            close(f)
        end
            
    end
    ntrial_now = nan;
end
totcounts = totcounts/totaltrials;
totcountsdiff = totcountsdiff/totaltrials;

f_average = figure;
subplot(1,2,1)
plot(bins, totcounts)
title('Activity')
xlabel('\Delta F / F')
ylabel('fraction of cells')

subplot(1,2,2)
plot(bins, totcountsdiff)
title('Difference in Activity')
xlabel('\Delta (\Delta F / F)')
ylabel('fraction of cells')