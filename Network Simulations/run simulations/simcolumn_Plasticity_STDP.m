function [spikepairs, Am] = simcolumn_Plasticity_STDP(spikepairs, spikeidx, celltypes, Am, alpha)
% implementing STDP plasticity rule for synapse
% currently dont model gain/loss of synaptic contact
% the parameter, spikepairs, is N-by-2 matrix, in which N is number of
% synapses to be simulated
% preSpikes is N-by-1 binary vector indicating current presynaptic activity
% postSpike is N-by-1 indication postsynaptic activity
% Am is synaptic weight matrix (N-by-1)
% celltypes is a N-by-4 matrix, with postType-postEI-preType-preEI as each
% column
% alpha is used to modulate learning rate; actual learning rate is A*alpha,
% with A as original learning rate

% an easy way of implementing is set two exponential decay function to
% model increase/decrease of synaptic weights; however need to keep track
% of pre/post synaptic firing

% check for spikes in both pre- and post- synaptic site
idx = spikeidx((spikepairs(spikeidx,3) + spikepairs(spikeidx,4)) == 2);

% modify Am based on plasticity rules
%% spike-timing plasticity is cell-type specific; set specific rules for different types of connections
% celltype is organized as postType, postEI, preType, preEI
if ~isempty(idx)
    Ampost = Am(idx);
    Ampre = Am(idx);
    % timing difference between post and pre synaptic spikes
    % timing difference is post - pre
    diffT = spikepairs(idx, 2) - spikepairs(idx, 1);
    
    % excitatory synapses onto excitatory neurons
    % L4-L4 connections, Egger et al 1999
    id_EE_L4L4 = find(celltypes(idx, 1) < 3 & celltypes(idx, 3) < 3);
    Ampost(id_EE_L4L4) = STDP_EE_L4spiny(diffT(id_EE_L4L4), Ampre(id_EE_L4L4), alpha);
    
    % L4-L2/3 connections, celikel et al 2004
    id_EE_L4L23 = find(celltypes(idx, 1) == 5 & celltypes(idx, 3) < 3);
%     if isempty(id_EE_L4L23) == 0
%         keyboard
%     end
    Ampost(id_EE_L4L23) = STDP_EE_celikel2004(diffT(id_EE_L4L23), Ampre(id_EE_L4L23), alpha);
    
    % L2/3-L2/3 connections, Banerjee et al 2014
    id_EE_L23L23 = find(celltypes(idx, 1) == 5 & celltypes(idx, 3) == 5);
    Ampost(id_EE_L23L23) = STDP_EE_Banerjee2014(diffT(id_EE_L23L23), Ampre(id_EE_L23L23), alpha);
    
%     % excitatory synapses onto inhibitory neurons
%     % excitatory to FS interneuron
%     id_EI_FS = find((celltypes(idx, 1) == 3 | celltypes(idx, 1) == 6) & celltypes(idx, 4) == 1);
%     Ampost(id_EI_FS) = STDP_EI_Fs_Lu2007(diffT(id_EI_FS), Ampre(id_EI_FS), alpha);
%     
%     % excitatory to LTS (SOM+) interneuron
%     id_EI_LTS = find((celltypes(idx, 1) == 4 | celltypes(idx, 1) == 9) & celltypes(idx, 4) == 1);
%     Ampost(id_EI_LTS) = STDP_EI_LTS_Lu2007(diffT(id_EI_LTS), Ampre(id_EI_LTS), alpha);
%     
%     % no data available on other types of interneurons; given their firing
%     % properties, probabily similar to FS interneurons?
%     
%     % inhibitory synapses onto all other cells
%     id_I = find(celltypes(idx, 4) == -1);
%     % Haas et al 2006, entorhinal cortex
%     Ampost(id_I) = STDP_I_haas2006(diffT(id_I), Ampre(id_I), alpha);
    
    
    % modify Am matrix
    Am(idx) = Ampost;
    % remove the early spike from the spike pairs register
    % this is the implementation of Izhikevich 2003, where only neigboring
    % spike pairs are considered
    timecheck = spikepairs(idx, 2) - spikepairs(idx, 1);
    spikepairs(idx(timecheck <=0), 2) = 0;
    spikepairs(idx(timecheck <=0), 4) = 0;
    spikepairs(idx(timecheck > 0), 1) = 0;
    spikepairs(idx(timecheck > 0), 3) = 0;
end


%% different plasticity rules
    % about additive and multiplicative updating rules: additive
    % potentiation and depression cause competetion (song et al 2000);
    % additive potentiation and multiplicative depression do not cause
    % competetion (Rossum et al 2000); 
    
    % about using all spikes or nearest neigbor spikes: all spikes cause
    % monotonic changes depending on whether potentiation dominant or
    % depression dominat; nearest neigbor spikes leads to convex evolution
    % more similar to BCM rule
    
    % the implementation here is additive potentiation, multiplicative
    % depression and nearest neigbor spikes

    % use alpha to modulate learning rate; actual learning rate is A*alpha,
    % with A as original learning rate
    % excitatory to excitatory 
    % celikel et al 2004; L4-L2/3 connections
    function Am_post = STDP_EE_celikel2004(dTs, Ampre, alpha)
        % gmax is used in additive case; maybe better to set it to
        % different synapse outside of STDP function
        gmax = 3;
        Am_post = Ampre;
        % using rule from celikel et al 2004
        % pre-lead-post: 50-100 pairings; A_post = 0.4, tau = 18
        idd = find(dTs > 0);
%         % multiplicative
%         Am_post(idd) = Am_post(idd).*(1 + alpha*0.4/75*exp(-dTs(idd)/18));
        % additive
        Am_post(idd) = Am_post(idd) + gmax*alpha*0.4/75*exp(-dTs(idd)/18);
        % post-lead-pre: 50-100 pairings; A_pre = -0.24, tau = 49
        idd = find(dTs < 0);
        % multiplicative
        Am_post(idd) = Am_post(idd).*(1 - alpha*0.24/75*exp(dTs(idd)/49));
%         % additive
%         Am_post(idd) = Am_post(idd) - gmax*alpha*0.24/75*exp(dTs(idd)/49);
        
        % set min-max of allowed Am values
        Am_post = setMinMax(Am_post, 0, 3.5);
    end

    % L2/3 - L2/3 connections
    function Am_post = STDP_EE_Banerjee2014(dTs, Ampre, alpha)
        gmax = 3;
        Am_post = Ampre;
        % using rule from Banerjee et al 2014
        % pre-lead-post: 100 pairings; A_post = 0.53, tau = 18
        idd = find(dTs > 0);
%         % multiplicative
%         Am_post(idd) = Am_post(idd).*(1 + alpha*0.5/100*exp(-dTs(idd)/18));
        % additive
        Am_post(idd) = Am_post(idd) + gmax*alpha*0.5/100*exp(-dTs(idd)/18);
        % post-lead-ppre: 100 pairings; A_pre = -0.32, tau = 18
        idd = find(dTs < 0);
        % multiplicative
        Am_post(idd) = Am_post(idd).*(1 - alpha*0.32/100*exp(dTs(idd)/18));
%         % additive
%         Am_post(idd) = Am_post(idd) - gmax*alpha*0.32/100*exp(dTs(idd)/18);
        
        % set min-max of allowed Am values
        Am_post = setMinMax(Am_post, 0, 4.5);
    end


    function Am_post = STDP_EE_Froekeme2002(dTs, Ampre, alpha)
        gmax = 3;
        Am_post = Ampre;
        % using rule from Froekeme 2002
        % pre-lead-post: potentiation; A_post = 1.03, tau = 14
        idd = find(dTs > 0);
        % multiplicative
        Am_post(idd) = Am_post(idd).*(1 + alpha*1.03*exp(-dTs(idd)/14));
        % additive
        Am_post(idd) = Am_post(idd) + gmax*alpha*0.5/100*exp(-dTs(idd)/18);
        % post-lead-pre: depression; A_pre = -0.51, tau = 31
        idd = find(dTs < 0);
        Am_post(idd) = Am_post(idd).*(1 - alpha*0.51*exp(dTs(idd)/31));
        % additive
        Am_post(idd) = Am_post(idd)  - gmax*alpha*0.51*exp(dTs(idd)/31);
        
        % set min-max of allowed Am values
        Am_post = setMinMax(Am_post, 0, 5);
    end

    % L4 spiny-spiny connections
    function Am_post = STDP_EE_L4spiny(dTs, Ampre, alpha)
        % specific to L4 spiny stellate neuron connections; only show 
        % symmetric depressing
        Am_post = Ampre;
        gmax = 3;
        % symmetric; with 5*10 pairings A = -0.2, tau = 10
        % multiplicative
        Am_post = Am_post.*(1 - alpha*0.2/50*exp(-abs(dTs)/12));
%         % additive
%         Am_post = Am_post - gmax*alpha*0.2/50*exp(-abs(dTs)/12);
        
        % set min-max of allowed Am values
        Am_post = setMinMax(Am_post, 0, 4.5);
    end

    % excitatory to inhibitory
    % excitatory to fast spiking interneuron
    function Am_post = STDP_EI_Fs_Lu2007(dTs, Ampre, alpha)
        % asymmetric depression
        Am_post = Ampre;
        gmax = 5;
        % pre-lead-post: depression; with 12*5 parings A_post = -0.56, tau = 39
        idd = find(dTs > 0);
        % multiplicative
        Am_post(idd) = Am_post(idd).*(1 - alpha*0.56/60*exp(-dTs(idd)/39));
        % additive
%         Am_post(idd) = Am_post(idd) - gmax*alpha*0.56/60*exp(-dTs(idd)/39);
%         % post-lead-pre: depression; with 12*5 parings A_pre = -0.46, tau = 39
        idd = find(dTs < 0);
        % multiplicative
        Am_post(idd) = Am_post(idd).*(1 - alpha*0.46/60*exp(dTs(idd)/39));
%         % additive
%         Am_post(idd) = Am_post(idd) - alpha*0.46/60*exp(dTs(idd)/39);
        
        % set min-max of allowed Am values
        Am_post = setMinMax(Am_post, 0, 5);
    end

    % excitatory to LTS interneuron (most likely SOM+, given the
    % facilitating synapse)
    function Am_post = STDP_EI_LTS_Lu2007(dTs, Ampre, alpha)
        % seems very similar to EE synapse; not much data in the paper
        Am_post = STDP_EE_celikel2004(dTs, Ampre, alpha);
    end

    
    % inhibitory to all
    % woodin et al 2003 in hippicampus; mexican-hat function
    function Am_post = STDP_I_woodin2003(dTs, Ampre, alpha)
        % 150 parings
        
        
    end

    % Haas et al 2006 in entorhinal cortex; similar to EE but
    % near-synchronized spikes not affect the synapse
    function Am_post = STDP_I_haas2006(dTs, Ampre, alpha)
        % 600 pairings
        gmax = 5;
        Am_post = Ampre;
        % pre-lead-post: potentiation; using LSRE fit in the paper
        idd = find(dTs > 0);
%         % multiplicative
%         Am_post(idd) =  Am_post(idd).*(1 + alpha/600*2.29e-6*dTs(idd).^10.*exp(-1.10*(dTs(idd))));
        % additive
        Am_post(idd) =  Am_post(idd) + gmax*alpha/600*2.29e-6*dTs(idd).^10.*exp(-1.10*(dTs(idd)));
        % post-lead-pre: depression; A_pre
        idd = find(dTs < 0);
        % mulitiplicative
        Am_post(idd) =  Am_post(idd).*(1 - alpha/600*2.6e-7 * dTs(idd).^10.*exp(0.94*dTs(idd)));
%         % additive
%         Am_post(idd) =  Am_post(idd) - gmax*alpha/600*2.6e-7 * dTs(idd).^10.*exp(0.94*dTs(idd));
        
        % set min-max of allowed Am values
        Am_post = setMinMax(Am_post, -6, 0);
    end

    % perform max-min set
    function Am_post = setMinMax(Am_post, minv, maxv)
        Am_post(Am_post > maxv) = maxv;
        Am_post(Am_post < minv) = minv;
    end
    
end