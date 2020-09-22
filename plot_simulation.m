function plot_simulation(Nsim, initdata, ConData, savename_input, varargin)
% Loads the results of simulations, and plots
% INPUT:
% * initdata: structure with initial settings, in '-filename--_initialsettings.mat';
% * ConData: structure with connectivity data, in '-filename--_ConData.mat'
% * Nsim: # simulations to plot
% * 'simvec' (optional): vector with which simulations to use

% exaple use: 
% * plot_simulation(10, initdata, ConData) plots the first 10 simulations
% * plot_simulation(10, initdata, ConData, )

f = filesep;

%% Plot example
if ~exist('Nsim','var')
    Nsim = input('Which simulation(s) do you want to plot? (give array)');
end

simvec = 1:Nsim;
% Look for 'varargin' inputs
len = length(varargin);
% check "len" for even number
if mod(len,2) > 0
    error('Wrong arguments: must be name-value pairs.');
end
for i = 1:2:len
    switch lower(varargin{i})
        case 'simvec'
            simvec=varargin{i+1};
            Nsim = length(simvec);
        otherwise
            % neglect invalid option
            disp(['Ignoring invalid input ' varargin{i}])
    end
end


whiskers = ones(1,Nsim);
for nt = 1:Nsim
    disp(['Loading simulation ' num2str(simvec(nt)) ' for plotting'])
    try
        if initdata.Vthresdyn
            thresholdname = 'dynthreshold_set';
        else
            thresholdname = 'fixthreshold_set';
        end
        load([ConData.savefolder ConData.FnametoSave '_Simcolumn_' thresholdname '_' initdata.setVthres.type '_simulation_' num2str(simvec(nt))  ]);
        load([ConData.savefolder savename_input '_Thalamic_Spike_Trains.mat'  ]);
    catch
        if nt>1
            keyboard
        end
        savefolder = input('What is the folder where the simulations were saved?', 's');
        if ~strcmp(savefolder(end), f)
            savefolder(end+1) = f ;
        end
        animal = input('What was the animal name?','s');
        savename = input('What were the file names?', 's');
        thresholdtype = input('What kind of threshold was used? (d/f)', 's');
        if strcmp(thresholdtype, 'd')
            thresholdname = 'dynthreshold_set';
        else
            thresholdname = 'fixthreshold_set';
        end
        thresholdset = input('How were the thresholds set? (d/t/i)','s');
        if strcmp(thresholdset, 'd')
            thresholdsetname = 'distribution';
        elseif strcmp(thresholdset, 't')
            thresholdsetname = 'pertype';
        else
            thresholdsetname = 'individual';
        end
        try
            if ~isempty(animal)
                ConData.savefolder = savefolder;
                ConData.FnametoSave = [savename '_' animal];
                initdata.setVthres.type = thresholdsetname;
                load([ConData.savefolder ConData.FnametoSave '_Simcolumn_' thresholdname '_' initdata.setVthres.type '_simulation_' num2str(simvec(nt))  ]);
                load([ConData.savefolder ConData.FnametoSave '_Thalamic_Spike_Trains.mat'  ], 'WhiskerTrace');
            else
                ConData.savefolder = savefolder;
                ConData.FnametoSave = savename;
                initdata.setVthres.type = thresholdsetname;
                load([ConData.savefolder ConData.FnametoSave '_Simcolumn_' thresholdname '_' initdata.setVthres.type '_simulation_' num2str(simvec(nt))  ]);
                load([ConData.savefolder ConData.FnametoSave '_Thalamic_Spike_Trains.mat'  ], 'WhiskerTrace');
            end
        catch
            load([savefolder savename num2str(simvec(nt)) '.mat'  ]);
            savenameth = input('What were the file names for the thalamic spike trains?', 's');
            load([savefolder savenameth], 'WhiskerTrace');
        end
    end
    
    [barrelnr, ind_barrel] = sort(cellinfo_all(:,5)); % sort by barrel
    celltypevec = cellinfo_all(ind_barrel,4);
    Nbarrel = max(barrelnr);
    Nall = length(modelsc);
    
    neuronlist = randi(Nall, [1,4]);
%     neuronlist = [5138,672,676,6720];
    disp(['Plotting simulation ' num2str(simvec(nt)) ])
    figure
    for nn = 1:4
        time = (1:length(V(ind_barrel(neuronlist(nn)),:)))*simdata.timestep;
        subplot(2,2,nn)

        plot(time, V(ind_barrel(neuronlist(nn)),:))
        hold all
        try
            plot(time,U(ind_barrel(neuronlist(nn)),:))
        catch
            disp('Parameter U not saved')
        end
        if initdata.Vthresdyn
            try
                plot(time, VT(ind_barrel(neuronlist(nn)),:))
            catch
                disp('Parameter V_T not saved')
            end
            legend('V', 'u','V_T')
        else
            legend('V', 'u')
        end
        ylim([-80, 20])
        grid on
        box on
        title(['Neuron ',num2str(ind_barrel(neuronlist(nn)))])
    end
    
    
    figure
    if whiskers(nt)==1
        for na=1:2
            subplot(6,1,na)
            plot((1:length(WhiskerTrace.Recording{na,simvec(nt)}))*WhiskerTrace.binsize, WhiskerTrace.Recording{na,simvec(nt)}, 'k', 'LineWidth',2)
            xlim([1, length(WhiskerTrace.Recording{na,simvec(nt)})*WhiskerTrace.binsize])
            if na==1
                title('Whisker Angle')
            else
                title('Whisker Curvature')
            end
            set(gca, 'XGrid','on')
        end
    end

    [Nneuronth, ~] = size(inputspikes);
    subplot(3,1,2)
    mint = 1;
    maxt = 0;
    for nn = 1:Nneuronth
        hold all
        plot(inputspikes(nn,:), nn*ones(size(inputspikes(nn,:))), '.k')
        u = unique(inputspikes(nn,:)); % remove 0
        if length(u)>1
            if u(2)<mint
                mint = u(2);
            end
        end
        if max(inputspikes(nn,:))>maxt
            maxt = max(inputspikes(nn,:));
        end
    end
    NCB = 0;
    [barrelnrth, ~] = sort(cellinfo_input(:,5)); % sort by barrel
    for nb = 1:Nbarrel
        % solid lines between barrels
        NCB = NCB+sum(barrelnrth==nb);
        try
            plot([mint maxt],[NCB,NCB], '-b')
        catch
            plot([0 maxt],[NCB,NCB], '-b')
        end
    end
    ylim([1 Nneuronth])
    title('Thalamic input spikes')
    box on
    
    subplot(3,1,3)
    [Nneuron, ~] = size(modelspt);
    for nn = 1:Nneuron
        hold all
        plot(modelspt(ind_barrel(nn),:), nn*ones(size(modelspt(ind_barrel(nn),:))), '.k')
        u = unique(modelspt(ind_barrel(nn),:)); % remove 0
        if length(u)>1
            if u(2)<mint
                mint = u(2);
            end
        end
        if max(modelspt(ind_barrel(nn),:))>maxt
            maxt = max(modelspt(ind_barrel(nn),:));
        end
    end
    
    NCB = 0;
    NCB_old = 0;
    for nb = 1:Nbarrel
        % solid lines between barrels
        NCB = NCB+sum(barrelnr==nb);
        try
            plot([mint maxt],[NCB,NCB], '-b')
        catch
            plot([0 maxt],[NCB,NCB], '-b')
        end
        % dotted lines between cell-types
        NCT = NCB_old;
        for celltype = 1:15
            NCT = NCT+sum(celltypevec(NCB_old+1:NCB)==celltype);
            try
                plot([mint maxt],[NCT,NCT], '--r')
            catch
                plot([0 maxt],[NCT,NCT], '--r')
            end
        end
        NCB_old = NCB;
    end
    ylim([1 Nneuron])
    title('Model spikes')
    
    for na=2:3
        subplot(3,1,na)
        xlim([mint, maxt])
        xlabel('time (ms)')
        ylabel('Neuron number')
        set(gca, 'XGrid','on')
        box on
    end
end

end