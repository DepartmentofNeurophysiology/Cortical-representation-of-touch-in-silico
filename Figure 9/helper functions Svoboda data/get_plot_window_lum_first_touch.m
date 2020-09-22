function [SampleMatlum, reltrials, time] = get_plot_window_lum_first_touch(s, lumwindow, volume, plotyn)
% Get from a Svoboda data file (s) the luminescence data in window (lumwindow) of given volume (2-Nvol) around the first touch times 
% OUTPUT:
% * Samplelum = matrix of size NTrial x Ncell x size window with all samples of luminescence traces of the given volume around the first touch times
% * reltrials = vector with trials in this volume with both luminescence and whiskertraces
%               * relevant trial nrs (1) in valid trials (see M)
%               * trial ids (2) 
% * time = relevant time base
% if plotyn=1 give 'pepper and salt'-plot of evolution over time

% Example use:
% [SampleMatlum, reltrials] = get_plot_window_lum(s, [-200,500], 3, 1);

% NB uses get_first_touch

if nargin == 3
    plotyn = 0;
end

%% Get general data
M = get_first_touch(s, 0); % matrix with 1) trial numbers, 2) times (whisker base), 3) protraction (1) or retraction (2), 4) volume nr of first touches and 5) times (lum base)
[~, Ntrial] = size(M);

dtlum = s.timeSeriesArrayHash.value{volume}.time(2) - s.timeSeriesArrayHash.value{volume}.time(1);
lumwindown = round(lumwindow/dtlum);
nlumwindow = length(lumwindown(1):lumwindown(2));
n0lum = -round(lumwindow(1)/dtlum)+1;
[Ncell, ~] = size(s.timeSeriesArrayHash.value{volume}.valueMatrix);

%% For later classification
dtwh = s.timeSeriesArrayHash.value{1}.time(2) - s.timeSeriesArrayHash.value{1}.time(1);
a_trace = s.timeSeriesArrayHash.value{1}.valueMatrix(1, :); % whisker angle
v_trace = (a_trace(2:end)-a_trace(1:end-1))/dtwh;
v_trace = [nan, v_trace]; % whisker velocity


%% get luminescence data
relmn = find(M(4,:) == volume); % relevant columns in M
reltrials = [relmn; M(1,relmn)]; % trials for this 
[~, NT] = size(relmn);

SampleMatlum = nan*ones(NT, Ncell, nlumwindow);
for nt = 1:NT
    % find relevant time point in luminescence time
    ntouchl = find(s.timeSeriesArrayHash.value{volume}.time == M(5,relmn(nt)));
    if (ntouchl+lumwindown(1)>0) && (ntouchl+lumwindown(2)<=size(s.timeSeriesArrayHash.value{volume}.valueMatrix,2))
        SampleMatlum(nt,:,:) = s.timeSeriesArrayHash.value{volume}.valueMatrix(:,ntouchl+lumwindown(1):ntouchl+lumwindown(2));  
    elseif (ntouchl+lumwindown(1)<=0)
        % recording first trial does not start early enough for start window
        trace = s.timeSeriesArrayHash.value{volume}.valueMatrix(:,1:ntouchl+lumwindown(2));
        SampleMatlum(nt,:,end-size(trace,2)+1:end) = trace;
    elseif (ntouchl+lumwindown(2)>size(s.timeSeriesArrayHash.value{volume}.valueMatrix,2))
        % recording last trial stops before end window
        trace = s.timeSeriesArrayHash.value{volume}.valueMatrix(:,ntouchl+lumwindown(1):end);
        SampleMatlum(nt,:,1:size(trace,2)) = trace; 
    end
end

time = (lumwindown(1)+(0:nlumwindow-1))*dtlum;
% (lumwindown(1)+(nlt-1))*dtlum
if plotyn
    %% plot luminescence data

    nxy = ceil(sqrt(Ncell));
    naxis = ceil(sqrt(nlumwindow));

    for nt = 1:NT
        figure('Name',['Trial ',num2str(reltrials(nt))])
        for nlt = 1:nlumwindow
            Zt = nan*ones(nxy,nxy);
            Zt(1:Ncell) = SampleMatlum(nt,:,nlt);

            subplot(naxis, naxis,nlt)
            imagesc(squeeze(Zt)); 
            colormap jet; 
%             colorbar
            if lumwindown(1)+(nlt-1) == 0
                % Give extra info in title
                ntouchw = find(s.timeSeriesArrayHash.value{1}.time == M(2,relmn(nt)));
                if a_trace(ntouchw)>0 && v_trace(ntouchw)>0
                    title({'a+ v+';['T = ',num2str(time(nlt))]})
                elseif a_trace(ntouchw)>0 && v_trace(ntouchw)<0
                    title({'a+ v-';['T = ',num2str(time(nlt))]})
                elseif a_trace(ntouchw)<0 && v_trace(ntouchw)>0
                    title({'a- v+';['T = ',num2str(time(nlt))]})
                elseif a_trace(ntouchw)<0 && v_trace(ntouchw)<0
                    title({'a- v-';['T = ',num2str(time(nlt))]})
                else
                    disp('no angle/velocity data')
                    keyboard
                end
            else
                title(['T = ',num2str(time(nlt))])
            end
            caxis([0 4])
        end
    end
end