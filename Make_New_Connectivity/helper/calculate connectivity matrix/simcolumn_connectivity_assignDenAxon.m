function Neuron = simcolumn_connectivity_assignDenAxon(Neuron)
% assign the axon and dendritic parameters to each cell
% the axon-dendritic parameters are also in this function
% need to specify subfields as constrains:
% 1: 'home': the distribution is confined in home barrel column
% 2: 'soma': the axon mainly target postsynaptic neuron somata
% for axon and dendritic distributions,
% for axonal and dendritic distribution, use N-by-6 maxtrix, with each row
% specifys a semi-gaussian distribution with center (x-y-z), lateral span
% (x, y) and vertical span (z)
% maximum number of connections need to be set specificially for each type
% of pre- and post- synaptic neurons

% global Model_Space;
Model_Space = make_Model_Space();
L2 = Model_Space.L2;
L4 = Model_Space.L4;
Loca = Neuron.Loca; % location of cell soma

% type of cells in layer 2/3
if Neuron.Layer == 2
%%    
    if Neuron.Type == 1 % pyramidal neurons
        % currently only model basal dendrite
        if Neuron.subType == 11 % L2 pyramidal
            rdx = 40*rand(1) - 20;
            rdy = 40*rand(1) - 20;
            rdz = 40*getPositive((L2(2) - L2(1))/6 - (Loca(3) - L2(1)))/(L2(2) - L2(1))*6 + 20*rand(1); % move the center toward L2 center
            Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) + rdz, 210 + 15*randn(1), 210 + 15*randn(1), 70 + 5 * randn(1)], ...
                [Loca(1) + rdx, Loca(2) + rdy, Loca(3) + rdz, 150 + 15*randn(1), 150 + 15*randn(1), 70 + 5 * randn(1)], ...
                [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), L4(2) + 80 + 50*rand(1), 100 + 10*randn(1), 100 + 10*randn(1), 85 + 8 * randn(1)]};
            Neuron.AxonTarget = [2, 4, 5]; % specific dirtibution has specific target
            Neuron.AxonRestriction = {'none', 'home', 'none'};
            Neuron.AxonTargetRestriction = {'none', 'none', 'none'};
            Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1), 50 + 5*randn(1)]);
            rdz = 30*getPositive((L2(2) - L2(1))/6 - (Loca(3) - L2(1)))/(L2(2) - L2(1))*6 + 20*rand(1); % move the center toward L2 center
            Neuron.Dendrite = {[Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz, 30 + 3*randn(1), 30 + 3*randn(1), 35 + 3 * randn(1)], ...
                [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 30 + 3*randn(1), 30 + 3*randn(1), 35 + 3 * randn(1)], ...
                [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 30 + 3*randn(1), 30 + 3*randn(1), 35 + 3 * randn(1)], ...
                };
            Neuron.DendriteTarget = [2, 4, 0];
            Neuron.DendriteRestriction = {'none', 'none', 'none'};
            Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
            
        elseif Neuron.subType == 12 % L3 pyramidal
            rdx = 40*rand(1) - 20;
            rdy = 40*rand(1) - 20;
            rdz = 60*getPositive((Loca(3) - L2(1)) - 4*(L2(2) - L2(1))/5)/(L2(2) - L2(1))/5*4 - 10*rand(1); % move the center toward L2 center
            Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 110 + 10*randn(1), 130 + 10*randn(1), 100 + 10 * randn(1)], ...
                [Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 90 + 10*randn(1), 90 + 10*randn(1), 100 + 10 * randn(1)], ...
                [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), L4(2) + 80 + 50*rand(1), 100 + 10*randn(1), 100 + 10*randn(1), 85 + 8 * randn(1)]};
            Neuron.AxonTarget = [2, 4, 5]; % specific dirtibution has specific target
            Neuron.AxonRestriction = {'none', 'home', 'none'};
            Neuron.AxonTargetRestriction = {'none', 'none', 'none'};
            Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1), 50 + 5*randn(1)]);
            rdz = 60*getPositive((Loca(3) - L2(1)) - 4*(L2(2) - L2(1))/5)/(L2(2) - L2(1))/5*4 + 10*rand(1); % move the center toward L2 center
            Neuron.Dendrite = {[Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz, 30 + 3*randn(1), 30 + 3*randn(1), 35 + 3 * randn(1)], ...
                [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 30 + 3*randn(1), 30 + 3*randn(1), 35 + 3 * randn(1)], ...
                 [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 30 + 3*randn(1), 30 + 3*randn(1), 35 + 3 * randn(1)], ...
                };
            Neuron.DendriteTarget = [2, 4, 0];
            Neuron.DendriteRestriction = {'none', 'none', 'none'};
            Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
        else
            error(['Neuron subtype type ' num2str(Neuron.subType) 'in layer ' num2str(Neuron.Layer) ' is not defined'])
        end
%%
    elseif Neuron.Type == 2 % fast spikting PV+ interneuron
        rdx = 40*rand(1) - 20;
        rdy = 40*rand(1) - 20;
        rdz = 30*getPositive((Loca(3) - L2(1)) - 4*(L2(2) - L2(1))/5)/(L2(2) - L2(1))/5*4 + 10*randn(1); % move the center toward L2 center
        Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 95 + 15*randn(1), 95 + 15*randn(1), 110 + 12 * randn(1)], ...
            [Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz + 5*randn(1), 100 + 15*randn(1), 100 + 15*randn(1), 90 + 12 * randn(1)], ...
        };
        Neuron.AxonTarget = [2, 4]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'none', 'none'};
        Neuron.AxonTargetRestriction = {'soma', 'soma'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        rdz = 30*getPositive((Loca(3) - L2(1)) - 4*(L2(2) - L2(1))/5)/(L2(2) - L2(1))/5*4 + 10*randn(1); % move the center toward L2 center
        Neuron.Dendrite = {[Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz, 40 + 3*randn(1), 40 + 3*randn(1), 50 + 3 * randn(1)], ...
            [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 50 + 3 * randn(1)], ...
            [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 50 + 3 * randn(1)], ...
            };
        Neuron.DendriteTarget = [2, 4, 0];
        Neuron.DendriteRestriction = {'none', 'none', 'none'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
%%        
    elseif Neuron.Type == 3 % chandlier cell
        rdx = 40*rand(1) - 20;
        rdy = 40*rand(1) - 20;
        rdz = -50*getPositive((L2(2) - L2(1))/6 - (Loca(3) - L2(1)))/(L2(2) - L2(1))*6 + 10*randn(1); % move the center toward L4 center
        Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) + rdz, 70 + 5*randn(1), 70 + 5*randn(1), 65 + 5 * randn(1)], ...
            [Loca(1) + rdx, Loca(2) + rdy, Loca(3) + rdz + 5*randn(1), 60 + 5*randn(1), 60 + 5*randn(1), 55 + 5 * randn(1)], ...
        };
        Neuron.AxonTarget = [2, 4]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'none', 'none'};
        Neuron.AxonTargetRestriction = {'soma', 'soma'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        rdz = 50*getPositive((L2(2) - L2(1))/6 - (Loca(3) - L2(1)))/(L2(2) - L2(1))*6 + 10*randn(1); % move the center toward L2 center
        Neuron.Dendrite = {[Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz, 40 + 3*randn(1), 40 + 3*randn(1), 50 + 3 * randn(1)], ...
            [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 50 + 3 * randn(1)], ...
            [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 50 + 3 * randn(1)], ...
            };
        Neuron.DendriteTarget = [2, 4, 0];
        Neuron.DendriteRestriction = {'none', 'none', 'none'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
%%
    elseif Neuron.Type == 4 % bursting PV+ neuron, as in Blatow et al
        % however, there is debating on if this group of cell exist in
        % mouse
        rdx = 40*rand(1) - 20;
        rdy = 40*rand(1) - 20;
        rdz = 20*getPositive((L2(2) - L2(1))/6 - (Loca(3) - L2(1)))/(L2(2) - L2(1))*6 + 10*randn(1); % move the center toward L2 center
        Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 120 + 15*randn(1), 120 + 15*randn(1), 110 + 12 * randn(1)], ...
            [Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz + 5*randn(1), 120 + 15*randn(1), 120 + 15*randn(1), 110 + 12 * randn(1)], ...
        };
        Neuron.AxonTarget = [2, 4]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'none', 'none'};
        Neuron.AxonTargetRestriction = {'none', 'none'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        rdz = 10*getPositive((L2(2) - L2(1))/6 - (Loca(3) - L2(1)))/(L2(2) - L2(1))*6 + 10*randn(1); % move the center toward L2 center
        Neuron.Dendrite = {[Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz, 50 + 5*randn(1), 50 + 5*randn(1), 60 + 5 * randn(1)], ...
            [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 50 + 3 * randn(1)], ...
            [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 50 + 3 * randn(1)], ...
            };
        Neuron.DendriteTarget = [2, 4, 0];
        Neuron.DendriteRestriction = {'none', 'none', 'none'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
%%
    elseif Neuron.Type == 5 % SOM+ martinotii cell
        rdx = 40*rand(1) - 20;
        rdy = 40*rand(1) - 20;
        rdz = 60*getPositive((Loca(3) - L2(1)) - 4*(L2(2) - L2(1))/5)/(L2(2) - L2(1))/5*4 + 10*randn(1); % move the center toward L2 border
        Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 110 + 9*randn(1), 110 + 9*randn(1), 90 + 12 * randn(1)], ...
            [Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz + 5*randn(1), 70 + 7*randn(1), 70 + 7*randn(1), 90 + 9 * randn(1)], ...
            };
        Neuron.AxonTarget = [2, 4]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'none', 'none'};
        Neuron.AxonTargetRestriction = {'none', 'none'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        rdz = 10*getPositive((L2(2) - L2(1))/6 - (Loca(3) - L2(1)))/(L2(2) - L2(1))*6 + 10*randn(1); % move the center toward L2 center
        Neuron.Dendrite = {[Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz, 40 + 3*randn(1), 40 + 3*randn(1), 40 + 3 * randn(1)], ...
            };
        Neuron.DendriteTarget = [2];
        Neuron.DendriteRestriction = {'none'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
        
%%
    elseif Neuron.Type == 6 % SOM+ bitufted cell; note that it seems they are the same group as SOM+ Mar cells electrophysiologically, but may differ anatomically
        % currently not different from group 5
        rdx = 40*rand(1) - 20;
        rdy = 40*rand(1) - 20;
        rdz = 60*getPositive((Loca(3) - L2(1)) - 4*(L2(2) - L2(1))/5)/(L2(2) - L2(1))/5*4 + 10*randn(1); % move the center toward L2 border
        Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 110 + 9*randn(1), 110 + 9*randn(1), 90 + 12 * randn(1)], ...
            [Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz + 5*randn(1), 70 + 7*randn(1), 70 + 7*randn(1), 90 + 9 * randn(1)], ...
            };
        Neuron.AxonTarget = [2, 4]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'none', 'none'};
        Neuron.AxonTargetRestriction = {'none', 'none'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        rdz = 10*getPositive((L2(2) - L2(1))/6 - (Loca(3) - L2(1)))/(L2(2) - L2(1))*6 + 10*randn(1); % move the center toward L2 center
        Neuron.Dendrite = {[Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz, 40 + 3*randn(1), 40 + 3*randn(1), 50 + 3 * randn(1)], ...
            };
        Neuron.DendriteTarget = [2];
        Neuron.DendriteRestriction = {'none'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
        
%%
    elseif Neuron.Type == 7 % VIP+ bipolar/double bounquit cell (CR-)
        rdx = 40*rand(1) - 20;
        rdy = 40*rand(1) - 20;
        rdz = 40*getPositive((Loca(3) - L2(1)) - 4*(L2(2) - L2(1))/5)/(L2(2) - L2(1))/5*4 + 10*randn(1); % move the center toward L2 border
        Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 80 + 9*randn(1), 80 + 9*randn(1), 240 + 20 * randn(1)], ...
            [Loca(1) + rdx, Loca(2) + rdy, (L4(2) + L4(1))/2 + 20*randn(1), 60 + 7*randn(1), 60 + 7*randn(1), 120 + 10 * randn(1)], ...
            };
        Neuron.AxonTarget = [2, 4]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'none', 'none'};
        Neuron.AxonTargetRestriction = {'none', 'none'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        rdz = 10*getPositive((L2(2) - L2(1))/6 - (Loca(3) - L2(1)))/(L2(2) - L2(1))*6 + 10*randn(1); % move the center toward L2 center
        Neuron.Dendrite = {[Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz, 30 + 3*randn(1), 30 + 3*randn(1), 50 + 3 * randn(1)], ...
            [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz, 30 + 3*randn(1), 30 + 3*randn(1), 50 + 3 * randn(1)], ...
             [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 50 + 3 * randn(1)], ...
            };
        Neuron.DendriteTarget = [2, 4, 0];
        Neuron.DendriteRestriction = {'none', 'none', 'none'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
        
%%
    elseif Neuron.Type == 8 % VIP+ bipolar/double bounquit cell; merged with group 7
        rdx = 40*rand(1) - 20;
        rdy = 40*rand(1) - 20;
        rdz = 40*getPositive((Loca(3) - L2(1)) - 4*(L2(2) - L2(1))/5)/(L2(2) - L2(1))/5*4 + 10*randn(1); % move the center toward L2 border
        Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 80 + 9*randn(1), 80 + 9*randn(1), 200 + 45 * randn(1)], ...
            [Loca(1) + rdx, Loca(2) + rdy, (L4(2) + L4(1))/2 + 20*randn(1), 60 + 7*randn(1), 60 + 7*randn(1), 120 + 10 * randn(1)], ...
            };
        Neuron.AxonTarget = [2, 4]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'none', 'none'};
        Neuron.AxonTargetRestriction = {'none', 'none'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        rdz = 10*getPositive((L2(2) - L2(1))/6 - (Loca(3) - L2(1)))/(L2(2) - L2(1))*6 + 10*randn(1); % move the center toward L2 center
        Neuron.Dendrite = {[Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz, 30 + 3*randn(1), 30 + 3*randn(1), 50 + 3 * randn(1)], ...
            [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz, 30 + 3*randn(1), 30 + 3*randn(1), 50 + 3 * randn(1)], ...
             [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 50 + 3 * randn(1)], ...
            };
        Neuron.DendriteTarget = [2, 4, 0];
        Neuron.DendriteRestriction = {'none', 'none', 'none'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
        
%%
    elseif Neuron.Type == 9 % VIP+/CR+ bipolar cell
        rdx = 40*rand(1) - 20;
        rdy = 40*rand(1) - 20;
        rdz = 40*getPositive((Loca(3) - L2(1)) - 4*(L2(2) - L2(1))/5)/(L2(2) - L2(1))/5*4 + 10*randn(1); % move the center toward L2 border
        Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 70 + 9*randn(1), 70 + 9*randn(1), 190 + 40 * randn(1)], ...
            [Loca(1) + rdx, Loca(2) + rdy, (L4(2) + L4(1))/2 + 20*randn(1), 60 + 7*randn(1), 60 + 7*randn(1), 120 + 10 * randn(1)], ...
            };
        Neuron.AxonTarget = [2, 4]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'none', 'none'};
        Neuron.AxonTargetRestriction = {'none', 'none'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        rdz = 10*getPositive((L2(2) - L2(1))/6 - (Loca(3) - L2(1)))/(L2(2) - L2(1))*6 + 10*randn(1); % move the center toward L2 center
        Neuron.Dendrite = {[Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz, 30 + 3*randn(1), 30 + 3*randn(1), 50 + 3 * randn(1)], ...
            [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) - rdz, 30 + 3*randn(1), 30 + 3*randn(1), 50 + 3 * randn(1)], ...
             [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) - rdz + 5*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 50 + 3 * randn(1)], ...
            };
        Neuron.DendriteTarget = [2, 4, 0];
        Neuron.DendriteRestriction = {'none', 'none', 'none'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
        
%%
    elseif Neuron.Type == 10 % CR+/VIP- multipolar cell
        rdx = 40*rand(1) - 20;
        rdy = 40*rand(1) - 20;
        rdz = 20*getPositive((Loca(3) - L2(1)) - 4*(L2(2) - L2(1))/5)/(L2(2) - L2(1))/5*4 + 10*randn(1); % move the center toward L2 center
        Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 120 + 15*randn(1), 120 + 15*randn(1), 120 + 12 * randn(1)], ...
            [Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz + 5*randn(1), 110 + 15*randn(1), 110 + 15*randn(1), 100 + 12 * randn(1)], ...
        };
        Neuron.AxonTarget = [2, 4]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'none', 'none'};
        Neuron.AxonTargetRestriction = {'none', 'none'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        rdz = 10*getPositive((L2(2) - L2(1))/6 - (Loca(3) - L2(1)))/(L2(2) - L2(1))*6 + 10*randn(1); % move the center toward L2 center
        Neuron.Dendrite = {[Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz, 40 + 5*randn(1), 40 + 5*randn(1), 50 + 5 * randn(1)], ...
            [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 50 + 3 * randn(1)], ...
             [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 50 + 3 * randn(1)], ...
            };
        Neuron.DendriteTarget = [2, 4, 0];
        Neuron.DendriteRestriction = {'none', 'none', 'none'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
        
%%
    elseif Neuron.Type == 11 % neuroglia form cell
        rdx = 40*rand(1) - 20;
        rdy = 40*rand(1) - 20;
        rdz = 60*getPositive((Loca(3) - L2(1)) - 2*(L2(2) - L2(1))/3)/(L2(2) - L2(1))/5*4 + 10*randn(1); % move the center toward L2 center
        Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 90 + 10*randn(1), 90 + 10*randn(1), 80 + 7 * randn(1)], ...
            [Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz + 5*randn(1), 80 + 10*randn(1), 80 + 10*randn(1), 70 + 7 * randn(1)], ...
        };
        Neuron.AxonTarget = [2, 4]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'none', 'none'};
        Neuron.AxonTargetRestriction = {'none', 'none'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        rdz = 10*getPositive((L2(2) - L2(1))/6 - (Loca(3) - L2(1)))/(L2(2) - L2(1))*6 + 10*randn(1); % move the center toward L2 center
        Neuron.Dendrite = {[Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz, 50 + 5*randn(1), 50 + 5*randn(1), 60 + 5 * randn(1)], ...
            [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 40 + 3 * randn(1)], ...
            [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + rdz + 5*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 40 + 3 * randn(1)], ...
            };
        Neuron.DendriteTarget = [2, 4, 0];
        Neuron.DendriteRestriction = {'none', 'none', 'none'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
        
%%
    else
        error(['Neuron type ' num2str(Neuron.Type) 'in layer ' num2str(Neuron.Layer) ' is not defined'])
    end
    
    % type of neurons in layer 4
elseif Neuron.Layer == 4
%%     
    Bloca = Model_Space.barrels;
    Brow = Bloca.row{Bloca.lable == Neuron.BarrelLoca};
    Barc = Bloca.arc{Bloca.lable == Neuron.BarrelLoca};
    if Neuron.Type == 1 % spiny pyramidal cell
        % move x,y location toward barrel center
%         keyboard
        rdx = 10*randn(1);
        rdy = 10*randn(1);
        rdz = 205 + 50*(Loca(3) - L4(1))/(L4(2) - L4(1)) + 10*randn(1); % move the center toward L2 center
        Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 115 + 10*randn(1), 115 + 10*randn(1), 165 + 12 * randn(1)], ...
            [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + 10*randn(1), 135 + 10*randn(1), 135 + 10*randn(1), 95 + 10 * randn(1)], ...
            [Loca(1) + 10*randn(1), Loca(2) + 10*randn(1), Loca(3) + 250 + 10*randn(1), 90 + 10*randn(1), 90 + 10*randn(1), 80 + 10 * randn(1)], ...
            };
        Neuron.AxonTarget = [2, 4, 5]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'none', 'none', 'none'};
        Neuron.AxonTargetRestriction = {'none', 'none', 'none'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        
        rdx =  - (0.5 + 0.5*rand(1))*(Brow(2) - Brow(1))/4*...
            sign(Loca(1) - mean(Brow))*getPositive(abs(Loca(1) - mean(Brow)) - (Brow(2) - Brow(1))/4)/(Brow(2) - Brow(1))*4;
        rdy =  - (0.5 + 0.5*rand(1))*(Barc(2) - Barc(1))/4*...
            sign(Loca(2) - mean(Barc))*getPositive(abs(Loca(2) - mean(Barc)) - (Barc(2) - Barc(1))/4)/(Barc(2) - Barc(1))*4;
        rdz =  - 0.15*(Loca(3) - mean(L4));
        Neuron.Dendrite = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) + rdz + 10*randn(1), 40 + 5*randn(1), 40 + 5*randn(1), 40 + 5 * randn(1)], ...
            [Loca(1) + rdx + 10*randn(1), Loca(2) + rdy + 10*randn(1), Loca(3) + rdz, 40 + 3*randn(1), 40 + 3*randn(1), 40 + 3 * randn(1)], ...
            };
        Neuron.DendriteTarget = [0, 4];
        Neuron.DendriteRestriction = {'none', 'none'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
%%
    elseif Neuron.Type == 2 % spiny stallete cell
        rdx = 10*randn(1) - 0.15*(Loca(1) - mean(Brow));
        rdy = 10*randn(1) - 0.15*(Loca(2) - mean(Barc));
        rdz = 205 + 50*(Loca(3) - L4(1))/(L4(2) - L4(1)) + 10*randn(1);
        Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 110 + 10*randn(1), 110 + 10*randn(1), 155 + 12 * randn(1)], ...
            [Loca(1) + rdx + 10*randn(1), Loca(2) + rdy + 10*randn(1), Loca(3) + 10*randn(1), 120 + 10*randn(1), 120 + 10*randn(1), 85 + 10 * randn(1)], ...
            };
        Neuron.AxonTarget = [2, 4]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'seminone', 'home'};
        Neuron.AxonTargetRestriction = {'none', 'none'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        
        rdx =  - 20 - (0.5 + 0.5*rand(1))*(Brow(2) - Brow(1))/4*...
            sign(Loca(1) - mean(Brow))*getPositive(abs(Loca(1) - mean(Brow)) - (Brow(2) - Brow(1))/4)/(Brow(2) - Brow(1))*4;
        rdy =  - (0.5 + 0.5*rand(1))*(Barc(2) - Barc(1))/4*...
            sign(Loca(2) - mean(Barc))*getPositive(abs(Loca(2) - mean(Barc)) - (Barc(2) - Barc(1))/4)/(Barc(2) - Barc(1))*4;
        rdz = - 0.15*(Loca(3) - mean(L4));
        Neuron.Dendrite = {[Loca(1) + rdx + 10*randn(1), Loca(2) + rdy + 10*randn(1), Loca(3) + rdz + 10*randn(1), 40 + 5*randn(1), 40 + 5*randn(1), 40 + 5 * randn(1)], ...
            [Loca(1) + rdx +  10*randn(1), Loca(2) + rdx + 10*randn(1), Loca(3) + rdz, 40 + 3*randn(1), 40 + 3*randn(1), 40 + 3 * randn(1)], ...
            };
        Neuron.DendriteTarget = [0, 4];
        Neuron.DendriteRestriction = {'home', 'home'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
%%
    elseif Neuron.Type == 3 % PV+ fast spiking cell
        Neuron.TypeName = 'FsPV';
        % - subclassification required
        if Neuron.subType == 31 % this type has projection into layer 2/3
            rdx = - (0.5 + 0.5*rand(1))*(Brow(2) - Brow(1))/4*...
                sign(Loca(1) - mean(Brow))*getPositive(abs(Loca(1) - mean(Brow)) - (Brow(2) - Brow(1))/4)/(Brow(2) - Brow(1))*4;
            rdy = - (0.5 + 0.5*rand(1))*(Barc(2) - Barc(1))/4*...
                sign(Loca(2) - mean(Barc))*getPositive(abs(Loca(2) - mean(Barc)) - (Barc(2) - Barc(1))/4)/(Barc(2) - Barc(1))*4;
            rdz = 140 + 80*(Loca(3) - L4(1))/(L4(2) - L4(1)) + 10*randn(1); % move the center toward L2 center
            Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 90 + 10*randn(1), 90 + 10*randn(1), 150 + 12 * randn(1)], ...
                [Loca(1) + rdx + 10*randn(1), Loca(2) + rdy + 10*randn(1), Loca(3) + 10*randn(1), 115 + 10*randn(1), 115 + 10*randn(1), 85 + 10 * randn(1)], ...
                };
            Neuron.AxonTarget = [2, 4]; % specific dirtibution has specific target
            Neuron.AxonRestriction = {'semihome', 'home'};
            Neuron.AxonTargetRestriction = {'none', 'none'};
            Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
            rdz = - 0.15*(Loca(3) - mean(L4));
            Neuron.Dendrite = {[Loca(1) + 10*randn(1) + rdx, Loca(2) + 10*randn(1) + rdy, Loca(3) + rdz + 10*randn(1), 40 + 5*randn(1), 40 + 5*randn(1), 40 + 5 * randn(1)], ...
                [Loca(1) + rdx + 10*randn(1), Loca(2) + rdy + 10*randn(1), Loca(3) + rdz + 10*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 40 + 3 * randn(1)], ...
                };
            Neuron.DendriteTarget = [0, 4];
            Neuron.DendriteRestriction = {'home', 'home'};
            Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
            
        elseif Neuron.subType == 32 %projection is confined in the barrel
            rdx = - (0.5 + 0.5*rand(1))*(Brow(2) - Brow(1))/4*...
                sign(Loca(1) - mean(Brow))*getPositive(abs(Loca(1) - mean(Brow)) - (Brow(2) - Brow(1))/4)/(Brow(2) - Brow(1))*4;
            rdy = - (0.5 + 0.5*rand(1))*(Barc(2) - Barc(1))/4*...
                sign(Loca(2) - mean(Barc))*getPositive(abs(Loca(2) - mean(Barc)) - (Barc(2) - Barc(1))/4)/(Barc(2) - Barc(1))*4;
            rdz = - 0.15*(Loca(3) - mean(L4)); % move the center toward L2 center
            Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 110 + 10*randn(1), 110 + 10*randn(1), 85 + 12 * randn(1)], ...
                };
            Neuron.AxonTarget = [4]; % specific dirtibution has specific target
            Neuron.AxonRestriction = {'home'};
            Neuron.AxonTargetRestriction = {'none'};
            Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
            rdz = - 0.15*(Loca(3) - mean(L4));
            Neuron.Dendrite = {[Loca(1) + rdx + 10*randn(1), Loca(2) + rdy + 10*randn(1), Loca(3) + 10*randn(1), 40 + 5*randn(1), 40 + 5*randn(1), 40 + 5 * randn(1)], ...
                [Loca(1) + rdx + 10*randn(1), Loca(2) + rdy + 10*randn(1), Loca(3) + rdz + 10*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 40 + 3 * randn(1)], ...
                };
            Neuron.DendriteTarget = [0, 4];
            Neuron.DendriteRestriction = {'home', 'home'};
            Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
            
        elseif Neuron.subType == 33 %projection is confined in the barrel, but connects to L2/3
            rdx = - (0.5 + 0.5*rand(1))*(Brow(2) - Brow(1))/4*...
                sign(Loca(1) - mean(Brow))*getPositive(abs(Loca(1) - mean(Brow)) - (Brow(2) - Brow(1))/4)/(Brow(2) - Brow(1))*4;
            rdy = - (0.5 + 0.5*rand(1))*(Barc(2) - Barc(1))/4*...
                sign(Loca(2) - mean(Barc))*getPositive(abs(Loca(2) - mean(Barc)) - (Barc(2) - Barc(1))/4)/(Barc(2) - Barc(1))*4;
            rdz = - 0.15*(Loca(3) - mean(L4)); % move the center toward L2 center
            Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 110 + 10*randn(1), 110 + 10*randn(1), 85 + 12 * randn(1)], ...
                [Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 110 + 10*randn(1), 110 + 10*randn(1), 85 + 12 * randn(1)]};
            Neuron.AxonTarget = [4,2]; % specific dirtibution has specific target
            Neuron.AxonRestriction = {'home','semihome'};
            Neuron.AxonTargetRestriction = {'none','none'};
            Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
            rdz = - 0.15*(Loca(3) - mean(L4));
            Neuron.Dendrite = {[Loca(1) + rdx + 10*randn(1), Loca(2) + rdy + 10*randn(1), Loca(3) + 10*randn(1), 40 + 5*randn(1), 40 + 5*randn(1), 40 + 5 * randn(1)], ...
                [Loca(1) + rdx + 10*randn(1), Loca(2) + rdy + 10*randn(1), Loca(3) + rdz + 10*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 40 + 3 * randn(1)], ...
                };
            Neuron.DendriteTarget = [0, 4];
            Neuron.DendriteRestriction = {'home', 'home'};
            Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
            
        else
             error(['Neuron subtype ' num2str(Neuron.subType) 'in layer ' num2str(Neuron.Layer) ' is not defined'])
        end
            
%%
    elseif Neuron.Type == 4 % regular/bursting, PV- cell
        rdx =  - (0.5 + 0.5*rand(1))*(Brow(2) - Brow(1))/4*...
            sign(Loca(1) - mean(Brow))*getPositive(abs(Loca(1) - mean(Brow)) - (Brow(2) - Brow(1))/4)/(Brow(2) - Brow(1))*4;
        rdy =  - (0.5 + 0.5*rand(1))*(Barc(2) - Barc(1))/4*...
            sign(Loca(2) - mean(Barc))*getPositive(abs(Loca(2) - mean(Barc)) - (Barc(2) - Barc(1))/4)/(Barc(2) - Barc(1))*4;
        rdz =  - 0.15*(Loca(3) - mean(L4)); % move the center toward L2 center
        Neuron.Axon = {[Loca(1) + rdx, Loca(2) + rdy, Loca(3) - rdz, 110 + 10*randn(1), 110 + 10*randn(1), 85 + 12 * randn(1)], ...
            };
        Neuron.AxonTarget = [4]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'home'};
        Neuron.AxonTargetRestriction = {'none'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        rdz =  - 0.15*(Loca(3) - mean(L4));
        Neuron.Dendrite = {[Loca(1) + rdx + 10*randn(1), Loca(2) + rdy + 10*randn(1), Loca(3) + 10*randn(1), 40 + 5*randn(1), 40 + 5*randn(1), 40 + 5 * randn(1)], ...
            [Loca(1) + rdx + 10*randn(1), Loca(2) + rdy + 10*randn(1), Loca(3) + rdz + 10*randn(1), 40 + 3*randn(1), 40 + 3*randn(1), 40 + 3 * randn(1)], ...
            };
        Neuron.DendriteTarget = [0, 4];
        Neuron.DendriteRestriction = {'home', 'home'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
        
%%
    else
        error(['Neuron type ' num2str(Neuron.Type) 'in layer ' num2str(Neuron.Layer) ' is not defined'])
    end
    
elseif Neuron.Layer == 0 % use Layer 0 as label of thalamus neuron
    % thalamic axons projection to L3, L4 and L5/6
    % they are only used as input cells, so do not need to assign cell
    % location
    % there is sub-classification of thalamic cells, as discussed in
    % Oberlaender et al 2012; see Furuta et al 2011, 
    % currently not implementing different subtypes
    
    % no need to specify dendrite
    % need barrel information
    Bloca = Model_Space.barrels;
    Brow = Bloca.row{Bloca.lable == Neuron.BarrelLoca};
    Barc = Bloca.arc{Bloca.lable == Neuron.BarrelLoca};
    
    % randomly generate center of projection based on barrel location and
    % L4 depth
    x = (Brow(2) - Brow(1))*rand(1) + Brow(1);
    y = (Barc(2) - Barc(1))*rand(1) + Barc(1);
    z = (L4(2) - L4(1))*rand(1) + L4(1); 
%   keyboard  
    % move the distribution located at the border towards barrel center
    x = x - (0.2 + 0.8*rand(1)).*(Brow(2) - Brow(1))/4.*...
        sign(x - mean(Brow)).*getPositive(abs(x - mean(Brow)) - (Brow(2) - Brow(1))/4)/(Brow(2) - Brow(1))*4;
    y = y - (0.2 + 0.8*rand(1)).*(Barc(2) - Barc(1))/4.*...
        sign(y - mean(Barc)).*getPositive(abs(y - mean(Barc)) - (Barc(2) - Barc(1))/4)/(Barc(2) - Barc(1))*4;
    z = z - 0.15*(z - mean(L4));
    
    % some thalamic axons are confined into one particular barrel
   if rand(1) < 0.85
%      keyboard
        Neuron.Axon = {[x, y, z, 115 + 8*randn(1), 115 + 8*randn(1), 70 + 5 * randn(1)], ...
            [x, y, z, 115 + 8*randn(1), 115 + 8*randn(1), 75 + 5 * randn(1)], ...
            };
        Neuron.AxonTarget = [2, 4]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'home', 'home'};
        Neuron.AxonTargetRestriction = {'none', 'none'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        Neuron.Dendrite = {[0, 0, 0, 50, 50, 50]
            };
        Neuron.DendriteTarget = [0];
        Neuron.DendriteRestriction = {'home'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
    
    else % some thalamic axon projects into neigboring barrel
        
        Neuron.Axon = {[x, y, z, 115 + 8*randn(1), 115 + 8*randn(1), 70 + 5 * randn(1)], ...
            [x, y, z, 135 + 8*randn(1), 135 + 8*randn(1), 75 + 5 * randn(1)], ...
            };
        Neuron.AxonTarget = [2, 4]; % specific dirtibution has specific target
        Neuron.AxonRestriction = {'none', 'none'};
        Neuron.AxonTargetRestriction = {'none', 'none'};
        Neuron.AxonMaxConn = round([100 + 10*randn(1), 10 + 2*randn(1)]);
        Neuron.Dendrite = {[0, 0, 0, 50, 50, 50]
            };
        Neuron.DendriteTarget = [0];
        Neuron.DendriteRestriction = {'home'};
        Neuron.DendriteMaxConn = round([100 + 10*randn(1), 250 + 10 * randn(1)]);
    end
end


%% internal functions
    function xpos = getPositive(x)
        xpos = x;
        xpos(xpos < 0) = 0;
    end
end