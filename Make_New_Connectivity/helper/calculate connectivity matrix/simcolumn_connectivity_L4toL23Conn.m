function [CM, DM] = simcolumn_connectivity_L4toL23Conn(preCells, postCells)
% calculate connectivity matrix, based on information provided in structure
% array Neurons

% first calcualte Axon/Dendrite overlapping for each potential presynaptic
% neurons, for a given post synaptic cell
% global Model_Space;
Model_Space = make_Model_Space();
step_size = Model_Space.stepsize;

AxonDendOverlap = zeros(length(postCells), length(preCells));
DistMat = zeros(length(postCells), length(preCells));
% PostLayer = zeros(length(Neurons), 1); % layer identifier of postsynaptic cells
% PreLayer = zeros(length(Neurons), 1); % layer identifier of presynaptic cells
PostType = zeros(length(postCells), 1); % type identifier of postsynaptic cells
PreType = zeros(1, length(preCells)); % type identifier of presynaptic cells

% keyboard

% use parellel computing
if verLessThan('matlab','8.2')
    if matlabpool('size') == 0
        matlabpool open 7
    end
else
    poolobj = gcp('nocreate'); % If no pool, do not create new one.
    if isempty(poolobj)
        pp = parcluster;
        nw = pp.NumWorkers;
        if nw>=7
            p = parpool(7);
        else
            p = parpool(nw);
        end        
    end
end

tic;
for i = 1:length(postCells) % post synaptic cell
%     disp(i)
    postcell = postCells{i};
    % calculate dendritic field, centered around the dendrite distribution
    % center
    Denfield = {};
    DendTarget = postcell.DendriteTarget;
    [X, Y, Z] = meshgrid(-200:step_size:200, -200:step_size:200, -200:step_size:200);
    for j = 1:length(postcell.Dendrite)
        Ref_point = postcell.Dendrite{j}(1:3);
        Denfield{j} = simcolumn_connectivity_DendFunc(Ref_point, postcell.Dendrite{j}, X, Y, Z);
        Denfield{j}(Denfield{j} < 0) = 0;
    end
%     PostLayer(i) = postcell.Layer;
    PostType(i) = postcell.Type;
    
    temp_Overlap = zeros(1, length(preCells));
    temp_Dist = zeros(1, length(preCells));
    
    parfor j = 1:length(preCells) % all presynaptic cells
        
        precell = preCells{j};
        if i == 1
%             PreLayer(j) = precell.Layer;
            PreType(j) = precell.Type;
        end
        
        if ~ismember(precell.Layer, DendTarget) || ~ismember(postcell.Layer, precell.AxonTarget)%presynaptic cell does not make projection to this layer
            continue
        end
        AxonTarget = precell.AxonTarget;
        % find the right Axon field to use
        axon = find(AxonTarget == postcell.Layer);
        
        % special restrictions
        if strcmpi(precell.AxonTargetRestriction{axon}, 'soma') ~= 1
            % center the field to the dendrite center
            Ref_point = postcell.Dendrite{find(DendTarget == precell.Layer)}(1:3);
            % axon field, only for the specific target region
            AxonField = simcolumn_connectivity_AxonFunc(Ref_point, precell.Axon{axon}, X, Y, Z);
            AxonField(AxonField < 0) = 0;
            
            Overlap = Denfield{find(DendTarget == precell.Layer)}.*AxonField;
%             temp_Overlap(i, j) = sum(Overlap(:));
        else
            % when the axon only target the soma, use soma location as
            % reference point
            Ref_point = postcell.Loca;
            % only need to calcuate density at reference point
            Overlap = simcolumn_connectivity_AxonFunc(Ref_point, precell.Axon{axon}, 0, 0, 0);
        end
        
        % in case the axon only target the home barrel
        if strcmpi(precell.AxonRestriction{axon}, 'home') == 1 && precell.BarrelLoca ~= postcell.BarrelLoca
            Overlap = 0.25 * Overlap;
        end
        
        % in case the axon has a bias target the home barrel
        if strcmpi(precell.AxonRestriction{axon}, 'semihome') == 1 && precell.BarrelLoca ~= postcell.BarrelLoca
            Overlap = 0.85 * Overlap;
        end
        
        
        temp_Overlap(j) = sum(Overlap(:));
        %distance between cells
        temp_Dist(j) = (sum((postcell.Loca - precell.Loca).^2)).^0.5;
    end
    AxonDendOverlap(i, :) = temp_Overlap;
    DistMat(i, :) = temp_Dist;
end
toc
% keyboard
% matlabpool close
delete(p);
%% ratio of axon/dendrite field located within model region
AxonRatio = [];
DendRatio = [];
for i = 1:length(postCells) % get DendRatio for postsynaptic cells
    try
       DendRatio(i) = postCells{i}.DendConnRatio(postCells{i}.DendriteTarget == 4);
    catch
        DendRatio(i) = 0;
    end
end
for i = 1:length(preCells) % get AxonRatio for presynaptic cells
    try
        AxonRatio(i) = preCells{i}.AxonConnRatio(preCells{i}.AxonTarget == 2);
    catch
        AxonRatio(i) = 0;
    end
end


%% convert overlapping into binary connectivity matrix
% keyboard
% post cells in row dimension and precells in column dimension
% note: to make the arrangement simplier the cell information should be
% sorted according to cell type
PostType_Unique = unique(PostType);
PreType_Unique = unique(PreType);
postN = [];
for i = 1:length(PostType_Unique)
    postN(i) = length(find(PostType == PostType_Unique(i)));
end
preN = [];
for i = 1:length(PreType_Unique)
    preN(i) = length(find(PreType == PreType_Unique(i)));
end
preIND = cumsum([0, preN]);
postIND = cumsum([0, postN]);

% keyboard
% connectivity is calculated between each cell type 
for i = 1:length(PreType_Unique)
    switch PreType_Unique(i)  
        case(1) % presynaptic cell is spiny pyramidal
            CM(:, preIND(i)+1:preIND(i+1)) = simConn_L4toL23_SPyrtoAll(AxonDendOverlap(:, preIND(i)+1:preIND(i+1)), ...
                postIND, PostType_Unique, AxonRatio(preIND(i)+1:preIND(i+1)), DendRatio);
        case(2) % presynaptic cell is spiny stallete
            CM(:, preIND(i)+1:preIND(i+1)) = simConn_L4toL23_SStetoAll(AxonDendOverlap(:, preIND(i)+1:preIND(i+1)), ...
                postIND, PostType_Unique, AxonRatio(preIND(i)+1:preIND(i+1)), DendRatio);
        case(3) % presynaptic cell is FsPV+ basket neurons
            CM(:, preIND(i)+1:preIND(i+1)) = -1*simConn_L4toL23_FsPVtoAll(AxonDendOverlap(:, preIND(i)+1:preIND(i+1)), ...
                postIND, PostType_Unique, AxonRatio(preIND(i)+1:preIND(i+1)), DendRatio, DistMat(:, preIND(i)+1:preIND(i+1)));
        case(4) % presynaptic cell is RS interneurons
            CM(:, preIND(i)+1:preIND(i+1)) = -1*simConn_L4toL23_RSNPtoAll(AxonDendOverlap(:, preIND(i)+1:preIND(i+1)), ...
                postIND, PostType_Unique, AxonRatio(preIND(i)+1:preIND(i+1)), DendRatio, DistMat(:, preIND(i)+1:preIND(i+1)));
    end
end

DM = DistMat.*abs(CM);

%% convert to sparse matrix
CM = sparse(CM);
DM = sparse(DM);
