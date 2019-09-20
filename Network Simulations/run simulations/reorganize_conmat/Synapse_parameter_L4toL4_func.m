function P4to4 = Synapse_parameter_L4toL4_func

% synapses parameters for L4-L4 connections
% currently 4 types of cells in L4: spiny pyramidal, spiny stallet, 
% Fast-spiking interneuron and non-fast-spiking interneuron

%% E to E
para4Eto4E.I0 = [1.2,1.1];      %Feldmeyer et al JPhysio 1999
para4Eto4E.tau1 = [17.8,3.1];
para4Eto4E.tau2 = [0.92,0.18];
para4Eto4E.Plas = [0.65,0.15];
para4Eto4E.Fail = [0.05,0.078];
para4Eto4E.Nb = 5;
para4Eto4E.va = 0.2;
para4Eto4E.CV = [0.37, 0.16];
para4Eto4E.delay = [190, 0.4];

%% E to I
%Beierlein et al JNP 2003
para4Eto4FS.I0 = [2.2,2.2];          
para4Eto4FS.tau1 = [7.9,1.9];
para4Eto4FS.tau2 = [0.37,0.11];
para4Eto4FS.Plas = [0.7,0.3];
para4Eto4FS.Fail = [0.03,0.08];
para4Eto4FS.Nb = 20;
para4Eto4FS.va = 0.4;
para4Eto4FS.CV = [0.27,0.13];
para4Eto4FS.delay = [190,0.3];

para4Eto4RS.I0 = [0.3,0.5];            
para4Eto4RS.tau1 = [11.9,2.9];
para4Eto4RS.tau2 = [0.86,0.48];
para4Eto4RS.Plas = [1.3,0.1];
para4Eto4RS.Fail = [0.57,0.35];
para4Eto4RS.Nb = 5;
para4Eto4RS.va = 0.3;
para4Eto4RS.CV = [0.8,0.4];
para4Eto4RS.delay = [190,0.4];

%% I to E
para4FSto4E.I0 = [1.1,0.8];      
para4FSto4E.tau1 = [24,4.8];
para4FSto4E.tau2 = [1.5,0.6];
para4FSto4E.Plas = [0.7,0.3];
para4FSto4E.Fail = [0.03,0.07];
para4FSto4E.Nb = 20;
para4FSto4E.va = -0.25;
para4FSto4E.CV = [0.25,0.11];
para4FSto4E.delay = [190,0.3];

para4RSto4E.I0 = [0.48,0.45];   
para4RSto4E.tau1 = [22.6,6.7];
para4RSto4E.tau2 = [2.1,0.5];
para4RSto4E.Plas = [1.0,0.1];
para4RSto4E.Fail = [0.29,0.26];
para4RSto4E.Nb = 5;
para4RSto4E.va = -0.5;
para4RSto4E.cv = 0.45;
para4RSto4E.CV = [0.41,0.21];
para4RSto4E.delay = [190,0.5];

%% I to I data is missing

%% put all data to a cell to be used in generateF2
P4Eto4E={para4Eto4E};
P4Eto4I={para4Eto4FS,para4Eto4RS};
P4Ito4E={para4FSto4E,para4RSto4E};
P4Ito4FS={para4FSto4E,para4RSto4E};
P4Ito4RS={para4FSto4E,para4RSto4E};

P4Ito4I={P4Ito4FS,P4Ito4RS};

P4to4 = {};
P4to4{1, 1} = para4Eto4E; %layer 4 SPyr to layer 4 SPyr
P4to4{1, 2} = para4Eto4E; %layer 4 SSte to layer 4 SPyr
P4to4{1, 3} = para4FSto4E; %layer 4 FsPV to layer 4 SPyr
P4to4{1, 4} = para4RSto4E; %layer 4 RSNP to layer 4 SPyr

P4to4{2, 1} = para4Eto4E; %layer 4 SPyr to layer 4 SSte
P4to4{2, 2} = para4Eto4E; %layer 4 SSte to layer 4 SSte
P4to4{2, 3} = para4FSto4E; %layer 4 FsPV to layer 4 SSte
P4to4{2, 4} = para4RSto4E; %layer 4 RSNP to layer 4 SSte

P4to4{3, 1} = para4Eto4FS; %layer 4 SPyr to layer 4 FsPV
P4to4{3, 2} = para4Eto4FS; %layer 4 SSte to layer 4 FsPV
P4to4{3, 3} = para4FSto4E; %layer 4 FsPV to layer 4 FsPV
P4to4{3, 4} = para4RSto4E; %layer 4 RSNP to layer 4 FsPV

P4to4{4, 1} = para4Eto4RS; %layer 4 SPyr to layer 4 FsPV
P4to4{4, 2} = para4Eto4RS; %layer 4 SSte to layer 4 FsPV
P4to4{4, 3} = para4FSto4E; %layer 4 FsPV to layer 4 FsPV
P4to4{4, 4} = para4RSto4E; %layer 4 RSNP to layer 4 FsPV