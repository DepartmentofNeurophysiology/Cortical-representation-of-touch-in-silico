function FirstTouchMat = get_first_touch(s, plotyn)
% Get from a Svoboda data file (s) the first touch times and plot(plotyn=1)
% or not (plotyn = 0)
% OUTPUT: FirstTouchMat = 3xNtrial file, where
% * Ntrial = number of valid trials (i.e. there was an object touch)
% * first row: trial numbers (ids)
% * second row: times
% * third row: protraction (1) or retraction (2)
% * fourth row: volume nr
% * fifth row: times in luminescance time

if nargin == 1
    plotyn = 0;
end

%[~, Ntrial] = size(s.trialIds);

[~, NV] = size(s.timeSeriesArrayHash.value);
NV = NV-1; % first is whisker data, rest is recorded 'volumes'
%% Load first touch times
toucht_pro = s.eventSeriesArrayHash.value{2}.eventTimes{1}; % touch times protraction
toucht_re  = s.eventSeriesArrayHash.value{2}.eventTimes{2}; % touch times retraction

trialn_pro = s.eventSeriesArrayHash.value{2}.eventTrials{1}; % corresponding trial nrs
trialn_re  = s.eventSeriesArrayHash.value{2}.eventTrials{2};

firstn_pro = unique(trialn_pro); % unique trial nrs
firstn_re  = unique(trialn_re);

ntrial_pro = length(firstn_pro); % number of unique trials
ntrial_re  = length(firstn_re);

% first touch times per trial
firstt_pro = nan*ones(1,ntrial_pro);  % protraction
firstt_re  = nan*ones(1,ntrial_re);   % retraction

for nt = 1:ntrial_pro
    ttemp = toucht_pro(trialn_pro==firstn_pro(nt));
    firstt_pro(nt) = ttemp(1);
end

for nt = 1:ntrial_re
    ttemp = toucht_re(trialn_re==firstn_re(nt));
    firstt_re(nt) = ttemp(1);
end

trials_uni = unique([firstn_pro, firstn_re]); % NB unique sorts too
N = length(trials_uni);
FirstTouchMat = nan*ones(5,N);
for nn=1:N
    current_trial = trials_uni(nn);
    cp = find(firstn_pro== current_trial);
    cr = find(firstn_re == current_trial);
    if isempty(cp) && isempty(cr)
        error('something went wrong: no first touch times')
    elseif ~isempty(cp) && isempty(cr)
        % only 1 first touch time on protraction
        FirstTouchMat(1,nn) = current_trial;
        FirstTouchMat(2,nn) = firstt_pro(cp);
        FirstTouchMat(3,nn) = 1;
    elseif isempty(cp) && ~isempty(cr)
        % only 1 first touch time on retraction
        FirstTouchMat(1,nn) = current_trial;
        FirstTouchMat(2,nn) = firstt_re(cr);
        FirstTouchMat(3,nn) = 2;
    elseif ~isempty(cp) && ~isempty(cr)
        % both a first protraction and retraction time exist: pick first
        [mintime, minpr] = min([firstt_pro(cp), firstt_re(cr)]);
        FirstTouchMat(1,nn) = current_trial;
        FirstTouchMat(2,nn) = mintime;
        FirstTouchMat(3,nn) = minpr;
    end
    
    % find volume and relevant time point in luminescence time
    v = [];
    for ii=1:NV
       t = find(s.timeSeriesArrayHash.value{ii+1}.trial == current_trial);
       if ~isempty(t)
           v = [v,ii+1];
       end
    end
    if isempty(v)
        disp(['no luminescence recording for trial ', num2str(current_trial)])
        FirstTouchMat(4,nn) = nan;
        FirstTouchMat(5,nn) = nan;
    elseif length(v) == 1
        FirstTouchMat(4,nn) = v;
        [~, nlum] = min(abs(s.timeSeriesArrayHash.value{v}.time - FirstTouchMat(2,nn)));
        FirstTouchMat(5,nn) = s.timeSeriesArrayHash.value{v}.time(nlum);
    else
        disp(['multiple recordings for trial ', num2str(current_trial)])
        keyboard
    end
    
    
    
    clear cp cr mintime minpr
        
end

if plotyn
    figure
    plot(toucht_pro, trialn_pro, '.')
    hold all
    plot(toucht_re, trialn_re, '.')
    plot(firstt_pro, firstn_pro, 'o')
    plot(firstt_re, firstn_re, 'o')
    plot(FirstTouchMat(2,:), FirstTouchMat(1,:), 'x')
    xlabel('touch time')
    ylabel('trial')
    legend('protraction touch times','retraction touch times','first protraction touch times','first retraction touch times','first touch times')
end

end
