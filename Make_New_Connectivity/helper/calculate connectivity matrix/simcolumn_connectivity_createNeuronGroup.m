function Neuron = simcolumn_connectivity_createNeuronGroup(cellinfo, layer)
% create data struction, Neuron, in which each element is a neuron, with
% location and type information taken from anatomical reconstruction

% Model_Space = make_Model_Space();
Neuron = {};
% keyboard
for i = 1:size(cellinfo, 1)
    Neuron{i}.Loca = cellinfo(i, 1:3);
    Neuron{i}.Type = cellinfo(i, 4);
    % and also the layer information
    Neuron{i}.Layer = layer;
    % also assign neuron type name
    Neuron{i} = simcolumn_connectivity_mapNeuronType(Neuron{i});
    % maybe barrel identity as well
    Neuron{i}.BarrelLoca = cellinfo(i, 5);
    % assign axonal and dendritic parameters
    Neuron{i} = simcolumn_connectivity_assignDenAxon(Neuron{i});
    % calculate percent of axon/dendrite in model region
    Neuron{i} = simcolumn_connectivity_maxConn(Neuron{i});
end