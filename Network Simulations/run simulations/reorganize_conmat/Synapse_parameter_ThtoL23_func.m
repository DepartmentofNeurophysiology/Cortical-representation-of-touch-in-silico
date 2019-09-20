function PThto2 = Synapse_parameter_ThtoL23_func
%required input value
% the data is currently unavailable; here use Th-L4 data instead
% ------------------------------------------------------------------------
%% TH to E
parathEto4E.I0 = [0.9,1.1];       %Bruno 2006
parathEto4E.tau1 = [15.5,3.1];
parathEto4E.tau2 = [1.1,0.2];
parathEto4E.Plas = [0.7,0.1];
parathEto4E.Fail = [0.05,0.02];
parathEto4E.Nb = 5;
parathEto4E.va = 0.4;
parathEto4E.CV = [0.41,0.2];
parathEto4E.delay = [500,1.3];

%% TH to I
parathEto4FS.I0 = [2.1,1.1];       %Beierlein et al JNP 2003
parathEto4FS.tau1 = [8,1];
parathEto4FS.tau2 = [0.41,0.15];
parathEto4FS.Plas = [0.75,0.1];
parathEto4FS.Fail = [0.05,0.02];
parathEto4FS.Nb = 5;
parathEto4FS.va = 0.9;
parathEto4FS.CV = [0.41,0.2];
parathEto4FS.delay = [500,1.3];

parathEto4RS.I0 = [0.7,0.5];       %Beierlein et al JNP 2003
parathEto4RS.tau1 = [15.5,3];
parathEto4RS.tau2 = [1.12,0.48];
parathEto4RS.Plas = [0.75,0.1];
parathEto4RS.Fail = [0.2,0.1];
parathEto4RS.Nb = 3;
parathEto4RS.va = 0.7;
parathEto4RS.CV = [0.7,0.3];
parathEto4RS.delay = [500,1.3];

%% put the parameters into a single cell
PThto2 = {};
PThto2{1, 1} = parathEto4E;    % thalamic to layer 2/3 Pyr

PThto2{2, 1} = parathEto4FS;    % thalamic to layer 2/3 FsPV

PThto2{3, 1} = parathEto4FS;    % thalamic to layer 2/3 ChPV

PThto2{4, 1} = parathEto4RS;    % thalamic to layer 2/3 RSPV

PThto2{5, 1} = parathEto4RS;    % thalamic to layer 2/3 MarSOM

PThto2{6, 1} = parathEto4RS;    % thalamic to layer 2/3 MarSOM

PThto2{7, 1} = parathEto4RS;    % thalamic to layer 2/3 BipVIP

PThto2{8, 1} = parathEto4RS;    % thalamic to layer 2/3 BipVIP

PThto2{9, 1} = parathEto4RS;    % thalamic to layer 2/3 BipCR

PThto2{10, 1} = parathEto4RS;    % thalamic to layer 2/3 SbcCR

PThto2{11, 1} = parathEto4RS;    % thalamic to layer 2/3 NG

