function run_trials(networkData, initdata, simdata, inputspikes, seeds, WhPara, varargin)
% run_save_trials:
% * iterates over the conditions (Vrest, Vthreshold, V0, u0, # iterations per condition) and runs the simulation

% INPUTS:
% * networkData: connectivity data structure from reorganize_conmat
% * initdata: structure with initial conditions: Vrest(mean resting
% membrane potential; can be list); stdVrest; Vthres (mean threshold
% potential; stdVthres; V0 (mean starting membrane potential);v stdV0; u0
% (mean starting u); stdu0; Vthresdyn (dynamic threshold (1) or not (0));
% Vthresvar (whether to use a distribution around the mean (1) or not (0) for excitatory neurons with dynamic threshold)
% * simdata: structure with simulation data: Tsim (total length
% simulation), timestep, both in ms and Trials: iterations per condition (so 1:4 runs 4 iterations, 4
% (re)runs iteration 4)
% * inputspikes: cell structure with input spikes inputspikes{# iterations per condition}(Ncell x Nspikes_max) 
% * seeds (optional): structure with seeds to re-do simulations: initseed = array
% (same size as iterations, seed for initial conditions); runseed = array
% (same size as iterations, seed for running: synaptic failures and
% amplitude)
% * WhPara: structure with the input for direct whisker modulation for each
% trial 

if isfield(simdata, 'inputvec')
    if ~isempty(simdata.inputvec)
        inputvec = simdata.inputvec;
    else
        inputvec = 1:length(inputspikes);
        simdata.inputvec = inputvec;
    end
else
    % use all input traces
    inputvec = 1:length(inputspikes);
    simdata.inputvec = inputvec;
end

whattosave = {'initdata','simdata'};
    
nwts = length(whattosave);
if exist('seeds', 'var')
    nwts = nwts+1;
    whattosave{nwts} = 'seeds';
end

save([networkData.savefolder networkData.FnametoSave '_initialsettings'], whattosave{:}, '-v7.3');

% total # simulations
Nsim = length(inputvec)*length(initdata.Vrest)*initdata.setVthres.nsim*length(initdata.V0)*length(initdata.u0)*length(simdata.Trials); 
if exist('seeds', 'var')
    if isfield(seeds, 'initseed')
        if ~(length(seeds.initseed) == Nsim)
            disp('Warning: number of seeds for initial conditions does not match number of simulations; random seeds will be generated')
        end
    end
    if isfield(seeds, 'runseed')
        if ~(length(seeds.runseed) == Nsim)
            disp('Warning: number of seeds for running the simulation does not match number of simulations; random seeds will be generated')
        end
    end
end

%% loop over different conditions to run the simulation
nsim = 0;
for lp_ipspikes = inputvec
    disp(['Input trace ' num2str(lp_ipspikes)])
    for lp_vr=1:length(initdata.Vrest)
        for lp_vt=1:initdata.setVthres.nsim
            for lp_v0 = 1:length(initdata.V0)
                for lp_u0 = 1:length(initdata.u0)
                    for lp_nit=1:length(simdata.Trials)
                        nsim = nsim+1;
                        
                        inputspikes_now = inputspikes{lp_ipspikes};
                        if isempty(WhPara)
                            WhPara_thissim = [];
                        else
                            WhPara_thissim = WhPara{lp_nit};
                        end

                        [Ninst, ~] = size(inputspikes_now);
                        if ~(Ninst == networkData.NIn)
                            error('Number of input spike trains does not correspond to number of input cells')
                        end
                        
                        initdata_thissim.nsim           = nsim;
                        initdata_thissim.Vrest          = initdata.Vrest(lp_vr);
                        initdata_thissim.stdVrest       = initdata.stdVrest(lp_vr);  
                        initdata_thissim.V0             = initdata.V0(lp_v0);
                        initdata_thissim.stdV0          = initdata.stdV0(lp_v0);
                        initdata_thissim.u0             = initdata.u0(lp_u0);
                        initdata_thissim.stdu0          = initdata.stdu0(lp_u0);
                        initdata_thissim.Vthresdyn      = initdata.Vthresdyn;
                        initdata_thissim.niterate       = simdata.Trials(lp_nit);
                        initdata_thissim.inputtracenr   = lp_ipspikes;
                        
                        disp(['Vr = ' num2str(initdata_thissim.Vrest) '; std Vr = ',num2str(initdata_thissim.stdVrest)])
                        disp(['V0 = ' num2str(initdata_thissim.V0) '; std V0 = ',num2str(initdata_thissim.stdV0)])
                        disp(['u0 = ' num2str(initdata_thissim.u0) '; std u0 = ',num2str(initdata_thissim.stdu0)])
                        disp(['trial ' num2str(initdata_thissim.niterate)]);
                        
                        if isfield(initdata, 'ExternalInput')
                            initdata_thissim.ExternalInput  = squeeze(initdata.ExternalInput(nsim,:,:));
                        end

                        if isfield(initdata, 'setindcell')
                            namevec = fieldnames(initdata.setindcell);
                            for ii = 1:length(namevec)
                                eval(['initdata_thissim.setindcell.' namevec{ii} '= initdata.setindcell.' namevec{ii} '(:,nsim);'])
                                if strcmp(namevec{ii}, 'Vthres')
                                    disp (['Vt= set per cell individually']);
                                end
                            end
                        end

                        initdata_thissim.setVthres.type = initdata.setVthres.type;
                        if strcmp(initdata.setVthres.type , 'distribution')
                            initdata_thissim.Vthres     = initdata.Vthres(lp_vt); 
                            initdata_thissim.stdVthres  = initdata.stdVthres(lp_vt);
                            disp([' Vt= ' num2str(initdata_thissim.Vthres)])
                        elseif strcmp(initdata.setVthres.type , 'pertype')
                            initdata_thissim.Vthresvar  = initdata.Vthresvar;
                            disp ('Vt= set per type')
                        end
                        
                        simdata_thissim = simdata;
                        if length(simdata.TSim)>1
                            try simdata_thissim.Tsim = simdata.TSim(lp_ipspikes);
                            catch
                                error('Please give TSim as an 1xNtrace array, or as a single value');
                            end
                        end
                            


                        if exist('seeds','var')
                            try
                                seeds_thissim.initseed = seeds.initseed(nsim);
                                seeds_thissim.runseed  = seeds.runseed(nsim);
                            catch
                                disp('Number of seeds does not fit number of simulations; adding random seed')
                                seeds_thissim.initseed = randi(2^32);
                                seeds_thissim.runseed = randi(2^32);
                                seeds.initseed(nsim) = seeds_thissim.initseed;
                                seeds.runseed(nsim)  = seeds_thissim.runseed;

                            end
                            % NB STDP changes the connectivity. The last
                            % version of the connectivity is stored in
                            % networkData and given back
                            networkData = simcolumn_indtrial(networkData, initdata_thissim, simdata_thissim, inputspikes_now, seeds_thissim, WhPara_thissim);
                        else
                            networkData = simcolumn_indtrial(networkData, initdata_thissim, simdata_thissim, inputspikes_now, [], WhPara_thissim);
                        end

                        close;
                    end % lp_nit
                end % lp_u0
            end % lp_v0
        end % lp_vt
    end % lp_vr
end % lp_ipspikes
