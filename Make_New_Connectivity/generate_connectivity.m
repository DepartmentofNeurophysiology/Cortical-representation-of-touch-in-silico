function generate_connectivity(barrelstruct, savefolder, savename)

%% use random process to generate cellinfo for new barrels
load('cellinfo_1645_assign.mat');
% % use the original cellinfo as barrel 1
% l23info(:,5) = 1;
% l4info(:,5) = 1;
% redistribute inhibitory neurons
l23info(:,4) = [];
l4info(:,4) = [];

barrels_l23 = {};
barrels_l4 = {};
paras_l23 = {};
paras_l4 = {};
thalamus = {};

[Nbx, Nby] = size(barrelstruct);
Nbarrel = Nbx*Nby;

barrels_l23 = cell(Nbarrel,1);
paras_l23   = cell(Nbarrel,1);
barrels_l4  = cell(Nbarrel,1);
paras_l4    = cell(Nbarrel,1);
thalamus    = cell(Nbarrel,1);
nb = 0;
for nbx = 1:Nbx
    for nby = 1:Nby
        nb = nb+1;
        [barrels_l23{nb}, paras_l23{nb}, barrels_l4{nb}, paras_l4{nb}, l23type, l4type] = ...
        simcolumn_connectivity_generateBarrel(l23info, l4info, [barrelstruct{nbx,nby}.xpos, barrelstruct{nbx,nby}.ypos], nb);

        % thalamic neurons
        thalamus{nb} = zeros(barrelstruct{nbx,nby}.Nthalamic, 3);
        thalamus{nb}(:,4) = 1;  % neuron type, currently not used
        thalamus{nb}(:,5) = nb; % barrel identifer
        thalamus{nb}(:,6) = 1;  % excitatory cells
    end
end    
%% group cell information together
l23info = [];
l23para = [];
l4info = [];
l4para = [];
thainfo = [];
for nb = 1:Nbarrel
    l23info = [l23info; barrels_l23{nb}];
    l23para = [l23para; paras_l23{nb}];
    l4info = [l4info; barrels_l4{nb}];
    l4para = [l4para; paras_l4{nb}];
    thainfo = [thainfo; thalamus{nb}];
end
% sort the matrix according to celltype
[~, idx] = sort(l23info(:,4));
l23info = l23info(idx, :);
l23para = l23para(idx, :);
[~, idx] = sort(l4info(:,4));
l4info = l4info(idx, :);
l4para = l4para(idx, :);
[~, idx] = sort(thainfo(:,4));
thainfo = thainfo(idx, :);

clear idx i

%% save the resulting cell information
save([savefolder 'cellinfo_' savename])

%% calculate connectivity matrix
disp('Calculating connectivity matrix')
% load model space
Model_Space = make_Model_Space();
% generate neuron groups
disp('Create neuron groups')
disp('L4')
L4 = simcolumn_connectivity_createNeuronGroup(l4info, 4);
disp('L23')
L23 = simcolumn_connectivity_createNeuronGroup(l23info, 2);
disp('Thalamus')
Tha = simcolumn_connectivity_createNeuronGroup(thainfo, 0);

% connectivity matrix
disp('Connectivity within cortex')
disp('L23 to L23')
[CMl23tol23, DMl23tol23]    = simcolumn_connectivity_L23toL23Conn(L23, L23);
disp('L4 to L23')
[CMl4tol23, DMl4tol23]      = simcolumn_connectivity_L4toL23Conn(L4, L23);
disp('L4 to L4')
[CMl4tol4, DMl4tol4]        = simcolumn_connectivity_L4toL4Conn(L4, L4);
disp('Connectivity with thalamus')
disp('Thalamus to L23')
[CMThtol23, DMThtol23]      = simcolumn_connectivity_ThtoL23Conn(Tha, L23);
disp('Thalamus to L4')
[CMThtol4, DMThtol4]        = simcolumn_connectivity_ThtoL4Conn(Tha, L4);

disp('Save')
save([savefolder 'CMDMs_' savename])

end