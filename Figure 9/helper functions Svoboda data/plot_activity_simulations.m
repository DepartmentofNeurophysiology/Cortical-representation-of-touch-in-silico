function [f_trial1_l23, f_average_l23, f_trial1_l4, f_average_l4] = plot_activity_simulations(CalciumPara, InputDataStruct, savefolder, savename, savename_sims, animaltype, edges, edgesav, showyn)

f = filesep;

%% find cell relevant simulations
list_temp = dir(savefolder);
list = {list_temp.name};
list = list([list_temp.isdir]==0);
loadlist = {};
counter = 0;
for nl=1:length(list)
    if length(list{nl})>=length(savename_sims)
        if strcmp(list{nl}(1:length(savename_sims)), savename_sims)
            if strcmp(list{nl}(end-14:end), 'calciumdata.mat')
                counter = counter+1;
                loadlist{counter} = list{nl};
            end
        end
    end
end

% sort loadlist
totaltrials = length(loadlist);
loadlist_temp = cell(size(loadlist));
for nl=1:totaltrials
    simplace = strfind(loadlist{nl}, 'simulation');
    underscoreplace = strfind(loadlist{nl}, '_');
    underscoreplace = underscoreplace(find(underscoreplace>simplace));
    trialnr = str2num(loadlist{nl}(underscoreplace(1)+1:underscoreplace(2)-1));
    loadlist_temp{trialnr} = loadlist{nl};
end
loadlist = loadlist_temp;
clear loadlist_temp

%% relevant cell info
load([savefolder 'cellinfo_' savename]);
[Nl4tot,~] =size(l4info);
[Nl23tot,~] =size(l23info);
% NB L4 cells are stored first
if strcmp(animaltype, 'mouse')
    l23info = rat_to_mouse_locations(l23info);
    l4info = rat_to_mouse_locations(l4info);
elseif strcmp(animaltype, 'rat')
    % do nothing
elseif isempty(animaltype)
    disp('no animal type defined, assume rat')
else
    disp('animal type not recognized, assume rat')
end

if isfield(InputDataStruct, 'depth')
    if ~isempty(InputDataStruct.depth)
        cellidtokeepl4 = (InputDataStruct.depth(1)<=l4info(:,3) & l4info(:,3)<=InputDataStruct.depth(2));       
        l4info = l4info(cellidtokeepl4,:);  
        cellidtokeepl23 = (InputDataStruct.depth(1)<=l23info(:,3) & l23info(:,3)<=InputDataStruct.depth(2));        
        l23info = l23info(cellidtokeepl23,:);           
    end
else
    cellidtokeepl23 = ones(size(l23info(:,3)));
    cellidtokeepl4 = ones(size(l4info(:,3)));
end
[Nl4,~] =size(l4info);
[Nl23,~] =size(l23info);

%% load relevant trials
% make normalized (cells) histogram over cells for each trial, and take
% average



savetrial = 1;
counter = 0;
deltat = 1./CalciumPara.frame_rate_c;
time = InputDataStruct.window.window(1):deltat:InputDataStruct.window.window(2);
totcountsl23 = zeros(1,length(edges)-1);
totcountsdiffl23= zeros(1,length(edges)-1);
totcountsl4 = zeros(1,length(edges)-1);
totcountsdiffl4= zeros(1,length(edges)-1);

[~, timebin0] = min(abs(time));
differencematl23 = nan(totaltrials, Nl23);
activitymatl23 = nan(totaltrials, Nl23);
differencematl4 = nan(totaltrials, Nl4);
activitymatl4 = nan(totaltrials, Nl4);
for nl = 1:length(loadlist)
    load([savefolder loadlist{nl}])
    
    lum_l4 = lummat(1:Nl4tot,:);
    lum_l4 = lum_l4(cellidtokeepl4,:);
    lum_l23 = lummat(Nl4tot+1:end,:);
    lum_l23 = lum_l23(cellidtokeepl23,:);
    
    counter = counter+1;
          
    if counter == savetrial
        f_trial1_l23 = figure('Name',['L23, Volume ' num2str(InputDataStruct.volume) ' trial ' num2str(nl)]);
        showyn_now = 1;
    else
        showyn_now = showyn;
        f_l23 = figure('Name',['L23, Volume ' num2str(InputDataStruct.volume) ' trial ' num2str(nl)]);
    end
    if sum(cellidtokeepl23)>0  
        [activitymatl23(nl,:), differencematl23(nl,:), countsl23, countsdiffl23] = plot_single_trial(l23info(:,1), l23info(:,2), lum_l23, timebin0, edges);
        totcountsl23 = totcountsl23+countsl23;
        totcountsdiffl23 = totcountsdiffl23+countsdiffl23;
    end
    
    
    if counter == savetrial
        f_trial1_l4 = figure('Name',['L4, Volume ' num2str(InputDataStruct.volume) ' trial ' num2str(nl)]);
        showyn_now = 1;
    else
        showyn_now = showyn;
        f_l4 = figure('Name',['L4, Volume ' num2str(InputDataStruct.volume) ' trial ' num2str(nl)]);
    end
    if sum(cellidtokeepl4)>0
        [activitymatl4(nl,:), differencematl4(nl,:), countsl4, countsdiffl4] = plot_single_trial(l4info(:,1), l4info(:,2), lum_l4, timebin0, edges);
        totcountsl4 = totcountsl4+countsl4;
        totcountsdiffl4= totcountsdiffl4+countsdiffl4;
    end
    
    
    if showyn_now
%         pause
    end
    if sum(cellidtokeepl23)>0
        if counter == savetrial
        else
            close(f_l23)
        end
    end
    if sum(cellidtokeepl4)>0
        if counter == savetrial
        else
            close(f_l4)
        end
    end

end
% NB Note that averaging over histograms is not the same as averaging over
% neurons and then making a histogram!
f_average_l23 = figure('Name',['L23 sims, Average Volume ' num2str(InputDataStruct.volume) ]);
if sum(cellidtokeepl23)>0  
    plot_average(l23info(:,1), l23info(:,2),activitymatl23, differencematl23, edgesav)
end
f_average_l4 = figure('Name',['L4 sims, Average Volume ' num2str(InputDataStruct.volume) ]);
if sum(cellidtokeepl4)>0
    plot_average(l4info(:,1), l4info(:,2),activitymatl4, differencematl4, edgesav)
end


    
