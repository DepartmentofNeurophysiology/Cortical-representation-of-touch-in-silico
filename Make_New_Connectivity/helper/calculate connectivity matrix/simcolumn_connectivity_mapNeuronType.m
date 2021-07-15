function Neuron = simcolumn_connectivity_mapNeuronType(Neuron)
% map numeric name and text-based name of each type of cell 

% global Model_Space;
Model_Space = make_Model_Space();

% type of cells in layer 2/3
if Neuron.Layer == 2
    if Neuron.Type == 1 % pyramidal neurons
        Neuron.TypeName = 'Pyr'; 
        % - further classification is required: L2 and L3 pyramidal
        L2 = Model_Space.L2;
        if Neuron.Loca(3) < L2(1) + (L2(2)-L2(1))/3
            Neuron.subType = 11;
            Neuron.subTypeName = 'L2Pyr';
        else
            Neuron.subType = 12;
            Neuron.subTypeName = 'L3Pyr';
        end
        
    elseif Neuron.Type == 2 % fast spikting PV+ interneuron
        Neuron.TypeName = 'FsPV';
    elseif Neuron.Type == 3 % chandlier cell
        Neuron.TypeName = 'Chan';
    elseif Neuron.Type == 4 % bursting PV+ neuron, as in Blatow et al
        % however, there is debating on if this group of cell exist in
        % mouse
        Neuron.TypeName = 'BsPV';
    elseif Neuron.Type == 5 % SOM+ martinotii cell
        Neuron.TypeName = 'Mar';
    elseif Neuron.Type == 6 % SOM+ bitufted cell; note that it seems they are the same group as SOM+ Mar cells electrophysiologically, but may differ anatomically
        % currently not different of group 5
        Neuron.TypeName = 'BitSOM';
    elseif Neuron.Type == 7 % VIP+ bipolar/double bounquit cell (CR-)
        Neuron.TypeName = 'BipVIP';
    elseif Neuron.Type == 8 % VIP+ bipolar/double bounquit cell; merged with group 7
        Neuron.TypeName = 'BipVIP';
    elseif Neuron.Type == 9 % VIP+/CR+ bipolar cell
        Neuron.TypeName = 'BipVC';
    elseif Neuron.Type == 10 % CR+/VIP- multipolar cell
        Neuron.TypeName = 'MulCR';
    elseif Neuron.Type == 11 % neuroglia form cell
        Neuron.TypeName = 'NG';
    else
        error(['Neuron type ' num2str(Neuron.Type) 'in layer ' num2str(Neuron.Layer) ' is not defined'])
    end

% type of neurons in layer 4
% distribution of ss and sp neurons seems depth-dependent; need to include
% this factor in future work
elseif Neuron.Layer == 4
    if Neuron.Type == 1 % spiny pyramidal cell
        Neuron.TypeName = 'sp';
    elseif Neuron.Type == 2 % spiny stallete cell
        Neuron.TypeName = 'ss';
    elseif Neuron.Type == 3 % PV+ fast spiking cell
        Neuron.TypeName = 'FsPV';
        % - subclassification required: barrel-L4 confined, and barrel
        % confined but projecting to L2/3 and L5
        if rand(1) < 1/4
            Neuron.subType = 31;
            Neuron.subTypeName = 'FsPV_Cin';
        elseif rand(1) > 2/3
            Neuron.subType = 33;% no seperate projection, but connnects to L2/3
            Neuron.subTypeName = 'FsPV_Cin1';
        else
            Neuron.subType = 32;
            Neuron.subTypeName = 'FsPV_Bin';
        end
        
    elseif Neuron.Type == 4 % regular/bursting, PV- cell
        Neuron.TypeName = 'RSNP';
    else 
        error(['Neuron type ' num2str(Neuron.Type) 'in layer ' num2str(Neuron.Layer) ' is not defined'])
    end
    
    
elseif Neuron.Layer == 0 % thalamus cell
    % thalamic axons projection to L3, L4 and L5/6
    % they are only used as input cells, so do not need to assign cell
    % location
    % there is sub-classification of thalamic cells, as discussed in
    % Oberlaender et al 2012; see Furuta et al 2011, 
    % currently not implementing different subtypes
    Neuron.TypeName = 'Tha';
    
end