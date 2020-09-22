function WhPara = WhiskerPara_direct_modulation(WhiskerTrace, barrelstruct, timestep_simulations, savefolder, savename)

% get base angles and calculate phase and amplitude for direct modulation



[Nbx, Nby] = size(barrelstruct);
Nbarrel = Nbx*Nby;
[Ndim, Ntrial] = size(WhiskerTrace.Recording);

%% check units, convert to degrees if in radians
Whiskerangles = cell(Ntrial,1);
baseanglepresent = 0;
for nd = 1:Ndim
    if strcmp(WhiskerTrace.quantity{nd}, 'base_angle')
        baseanglepresent  = 1;
        if strcmp(WhiskerTrace.unit{nd}, 'degree')
            for nt=1:Ntrial
                Whiskerangles{nt} = WhiskerTrace.Recording{nd,nt};
            end                
        elseif strcmp(WhiskerTrace.unit{nd}, 'radian')
            % convert to degrees
            for nt = 1:Ntrial
                Whiskerangles{nt} = WhiskerTrace.Recording{nd,nt}*180/pi;
            end
        else
            unit_angle = input('Were the base angles recorded in degrees (1) or radians (2)? ');
            if unit_angle == 1
            elseif unit_angle == 2
            else
                error('Direct modulation cannot be implemented: unit of base angle trace recordings not known');
            end
        end
    end
end
if ~baseanglepresent
    error('Direct modulation cannot be implemented: no base angle trace present in recordings');
end
    
WhPara = cell(Ntrial,1);
for nt = 1:Ntrial
    % resample the base angle to the step size used in simulation)
    f1 = 1/timestep_simulations;
    f2 = 1/WhiskerTrace.binsize; 
    while ~((rem(f1,1)==0) && (rem(f2,1)==0))
        f1 = f1*10;
        f2 = f2*10;
    end
    BaseAngleRS = resample(Whiskerangles{nt},f1 ,f2);
    % find location of positive peaks
    [~,Plocs] = findpeaks(BaseAngleRS, 'MinPeakProminence', 15, 'MinPeakDistance', 250, 'MinPeakWidth', 80);
    % location of negtive peaks
    [~,Nlocs] = findpeaks(max(BaseAngleRS) - BaseAngleRS, ...
        'MinPeakProminence', 15, 'MinPeakDistance', 250, 'MinPeakWidth', 80);
    
    % find each cycle; cycle starts from most retracted point
    % artifacts from resampling
    tmp = BaseAngleRS;
    tmp(end-20:end) = NaN;
    tmp(1:20) = NaN;
    % parameters to calculate
    Ams = zeros(size(tmp));
    AmsD = zeros(size(tmp));
    Pha = zeros(size(tmp));
    RetraSet = zeros(size(tmp));
    if Plocs(1) < Nlocs(1) % positive peak first; need look back to get starting point
        idx = find(tmp(1:Plocs(1)) == min(tmp(1:Plocs(1))), 1, 'first');
        % amplitude in full circle
        Ams(idx:Nlocs(1)-1) = tmp(Plocs(1)) - tmp(idx);
        % amplitude calculated seperately for protraction and retraction
        AmsD(idx:Plocs(1)-1) = tmp(Plocs(1)) - tmp(idx);
        AmsD(Plocs(1):Nlocs(1)-1) = tmp(Plocs(1)) - tmp(Nlocs(1));
        RetraSet(idx:Nlocs(1)-1) = tmp(idx);
        % phase is evenly distributed in time in a circle; by convention
        % protraction end points are pi
        Pha(idx:Plocs(1)) = linspace(-pi, 0, Plocs(1)-idx+1);
        Pha(Plocs(1):Nlocs(1)) = linspace(0, pi, Nlocs(1)-Plocs(1)+1);
        % remove first Plocs, to simplify the coding
        Plocs(1) = [];
    else % negtive peak first; can ignore the data before first negtive peak
        
        
    end
    
    for j = 1:length(Nlocs) - 1
        Ams(Nlocs(j):Nlocs(j+1)-1) = tmp(Plocs(j)) - tmp(Nlocs(j));
        AmsD(Nlocs(j):Plocs(j)-1) = tmp(Plocs(j)) - tmp(Nlocs(j));
        AmsD(Plocs(j):Nlocs(j+1)-1) = tmp(Plocs(j)) - tmp(Nlocs(j+1));
        RetraSet(Nlocs(j):Nlocs(j+1)-1) = tmp(Nlocs(j));
        Pha(Nlocs(j):Plocs(j)) = linspace(-pi, 0, Plocs(j)-Nlocs(j)+1);
        Pha(Plocs(j):Nlocs(j+1)) = linspace(0, pi, Nlocs(j+1)-Plocs(j)+1);
    end
    % last imcomplete circle
    idx = Nlocs(end) + find(tmp(Nlocs(end):end) == max(tmp(Nlocs(end):end)), 1, 'first');
    Ams(Nlocs(end):end) = tmp(idx) - tmp(Nlocs(end));
    AmsD(Nlocs(end):idx-1) = tmp(idx) - tmp(Nlocs(end));
    AmsD(idx:end) = tmp(idx) - min(tmp(idx:end));
    RetraSet(Nlocs(end):end) = tmp(Nlocs(end));
    Pha(Nlocs(end):idx) = linspace(-pi, 0, idx-Nlocs(end)+1);
    Pha(idx:end) = linspace(0, pi, length(tmp)-idx+1);
    
    WhPara{nt}.Phase = Pha;
    WhPara{nt}.RetraSetPoints = RetraSet;
    WhPara{nt}.Amplitude = Ams;
    WhPara{nt}.Amplitude_PR = AmsD;
    WhPara{nt}.Baseangles_degrees = BaseAngleRS;
    
end
save([savefolder savename '_WhiskerModulation'], 'WhPara');
