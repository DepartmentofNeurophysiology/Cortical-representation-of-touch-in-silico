function P4to2 = Synapse_parameter_L4toL23_func
%%store synaptic efficacy parameters for the network model
% ------------------------------------------------------------------------
%% layer 4 to layer 2/3 E to E
para4E2E.I0 = [0.70,0.6];       %Feldmeyer 2002
para4E2E.tau1 = [12.7,3.5];
para4E2E.tau2 = [0.6,0.3];
para4E2E.Plas = [0.8,0.28];
para4E2E.Fail = [0.049,0.08];
para4E2E.STD.U = [0.7,0.1];
para4E2E.STD.F = [1000,100];
para4E2E.STD.D = [800,80];
para4E2E.Nb = 5;
para4E2E.va = 0.15;
para4E2E.delay = [190, 0.4];

%% E to I
para4E2I.I0 = [1.2,1.1];       %Helmstaedter 2008
para4E2I.tau1 = [14.7,7];
para4E2I.tau2 = [0.4,0.1];
para4E2I.Plas = [0.8,0.1];
para4E2I.Fail = [0.049,0.04];
para4E2I.STD.U = [0.7,0.1];
para4E2I.STD.F = [1000,100];
para4E2I.STD.D = [800,80];
para4E2I.Nb = 3;
para4E2I.va = 0.2;
para4E2I.delay = [190, 0.4];

para4E2lateral.I0 = [1.4,1.3];       %Helmstaedter 2008
para4E2lateral.tau1 = [11,4];        %CR+ basket
para4E2lateral.tau2 = [0.4,0.2];
para4E2lateral.Plas = [0.92,0.2];
para4E2lateral.Fail = [0.13,0.2];
para4E2lateral.STD.U = [0.45,0.1];
para4E2lateral.STD.F = [1000,100];
para4E2lateral.STD.D = [800,80];
para4E2lateral.Nb = 3;
para4E2lateral.va = 0.55;
para4E2lateral.delay = [190, 0.2];

para4E2lateral1.I0 = [1.1,0.7];       %Helmstaedter 2008 %lateral1 is basket cell 
para4E2lateral1.tau1 = [10.3,2];        %assigned to PV+ FS neurons
para4E2lateral1.tau2 = [0.4,0.1];
para4E2lateral1.Plas = [0.84,0.17];
para4E2lateral1.Fail = [0.03,0.01];
para4E2lateral1.STD.U = [0.6,0.1];
para4E2lateral1.STD.F = [1000,100];
para4E2lateral1.STD.D = [800,80];
para4E2lateral1.Nb = 3;
para4E2lateral1.va = 0.3;
para4E2lateral1.delay = [190, 0.2];

para4E2trans.I0 = [1.3,0.9];       %Helmstaedter 2008
para4E2trans.tau1 = [14,4];        %assigned to VIP+ and CR+ bitufted
para4E2trans.tau2 = [0.6,0.2];
para4E2trans.Plas = [1,0.6];
para4E2trans.Fail = [0.13,0.09];
para4E2trans.STD.U = [0.32,0.08];
para4E2trans.STD.F = [62,12];
para4E2trans.STD.D = [144,28];
para4E2trans.Nb = 3;
para4E2trans.va = 0.55;
para4E2trans.cv = 0.38;
para4E2trans.delay = [190, 0.5];

para4E2local.I0 = [0.96,0.6];       %Helmstaedter 2008
para4E2local.tau1 = [15,7];         %chandlier neurons
para4E2local.tau2 = [0.4,0.2];
para4E2local.Plas = [1,0.2];
para4E2local.Fail = [0.13,0.09];
para4E2local.STD.U = [0.32,0.08];
para4E2local.STD.F = [62,12];
para4E2local.STD.D = [144,28];
para4E2local.Nb = 3;
para4E2local.va = 0.6;
para4E2local.cv = 0.35;
para4E2local.delay = [190, 0.2];

para4E2local2.I0 = [0.59,0.15];       %Helmstaedter 2008 %local2 is NGF cell
para4E2local2.tau1 = [13,4];          %NGF cell
para4E2local2.tau2 = [0.42,0.09];
para4E2local2.Plas = [0.7,0.17];
para4E2local2.Fail = [0.17,0.1];
para4E2local2.STD.U = [0.7,0.1];
para4E2local2.STD.F = [1000,100];
para4E2local2.STD.D = [800,80];
para4E2local2.Nb = 3;
para4E2local2.va = 0.48;
para4E2local2.cv = 0.4;
para4E2local2.delay = [190, 0.5];

para4E2local3.I0 = [0.3,0.1];       %Helmstaedter 2008 %local3 is bitufted cell
para4E2local3.tau1 = [22,9];        %assigned to SOM+ bitufted neurons
para4E2local3.tau2 = [0.42,0.09];
para4E2local3.Plas = [1.2,0.06];
para4E2local3.Fail = [0.51,0.075];
para4E2local3.STD.U = [0.16,0.1];
para4E2local3.STD.F = [376,110];
para4E2local3.STD.D = [42,11];
para4E2local3.Nb = 3;
para4E2local3.va = 0.8;
para4E2local3.cv = 0.5;
para4E2local3.delay = [190, 0.5];

para4E2layer1.I0 = [0,0];       %Helmstaedter 2008 %not connected
para4E2layer1.tau1 = [15,7];
para4E2layer1.tau2 = [0.4,0.2];
para4E2layer1.Plas = [1,0.2];
para4E2layer1.Fail = [1,0];
para4E2layer1.STD.U = [0.16,0.1];
para4E2layer1.STD.F = [376,110];
para4E2layer1.STD.D = [42,11];
para4E2layer1.Nb = 0;
para4E2layer1.delay = [190, 0.5];

%% I to E
% only known possible connection is from one subtype of FsPV; synaptic
% parameters not known, use I-E connection properties from L4-L4
% connections
para4FSto4E.I0 = [1.1,0.8];      
para4FSto4E.tau1 = [24,4.8];
para4FSto4E.tau2 = [1.5,0.6];
para4FSto4E.Plas = [0.7,0.3];
para4FSto4E.Fail = [0.03,0.07];
para4FSto4E.Nb = 20;
para4FSto4E.va = -1;
para4FSto4E.CV = [0.25,0.11];
para4FSto4E.delay = [190, 0.7];

para4RSto4E.I0 = [0.48,0.45];   
para4RSto4E.tau1 = [22.6,6.7];
para4RSto4E.tau2 = [2.1,0.5];
para4RSto4E.Plas = [1.0,0.1];
para4RSto4E.Fail = [0.29,0.26];
para4RSto4E.Nb = 5;
para4RSto4E.va = -1;
para4RSto4E.CV = [0.41,0.21];
para4RSto4E.delay = [190, 0.7];


%% put parameters together
P4E2E={para4E2E};
P4E2I={para4E2lateral1,para4E2local,para4E2lateral1,para4E2layer1,para4E2local3,para4E2trans,para4E2trans,para4E2trans,para4E2lateral,para4E2local2};
P4I2E={para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E};
P4I2DBC_VIP={para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E};
P4I2Bip_VIP={para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E};
P4I2SBC_CR={para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E};
P4I2Bip_CR={para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E};
P4I2BS_PV={para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E};
P4I2CH_PV={para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E};
P4I2Res_PV={para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E};
P4I2Mar_SOM={para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E};
P4I2NGC_AC={para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E};
P4I2Bit_SOM={para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E,para4E2E};
P4I2I={P4I2BS_PV,P4I2CH_PV,P4I2Res_PV,P4I2Mar_SOM,P4I2Bit_SOM,P4I2DBC_VIP,P4I2Bip_VIP,P4I2Bip_CR,P4I2SBC_CR,P4I2NGC_AC};

% different type of organization
P4to2 = {};
P4to2{1, 1} = para4E2E; % layer 4 SPyr to layer 2/3 E
P4to2{1, 2} = para4E2E; % layer 4 SSte to layer 2/3 E
P4to2{1, 3} = para4FSto4E; % layer 4 FsPV to layer 2/3 E
P4to2{1, 4} = para4RSto4E; % layer 4 RsNP to layer 2/3 E

P4to2{2, 1} = para4E2local; % layer 4 SPyr to layer 2/3 FsPV
P4to2{2, 2} = para4E2local; % layer 4 SSte to layer 2/3 FsPV
P4to2{2, 3} = para4FSto4E; % layer 4 FsPV to layer 2/3 FsPV
P4to2{2, 4} = para4RSto4E; % layer 4 RsNP to layer 2/3 FsPV

P4to2{3, 1} = para4E2local; % layer 4 SPyr to layer 2/3 CHPV
P4to2{3, 2} = para4E2local; % layer 4 SSte to layer 2/3 CHPV
P4to2{3, 3} = para4FSto4E; % layer 4 FsPV to layer 2/3 CHPV
P4to2{3, 4} = para4RSto4E; % layer 4 RsNP to layer 2/3 CHPV

P4to2{4, 1} = para4E2lateral1; % layer 4 SPyr to layer 2/3 ResPV
P4to2{4, 2} = para4E2lateral1; % layer 4 SSte to layer 2/3 ResPV
P4to2{4, 3} = para4FSto4E; % layer 4 FsPV to layer 2/3 ResPV
P4to2{4, 4} = para4RSto4E; % layer 4 RsNP to layer 2/3 ResPV

P4to2{5, 1} = para4E2layer1; % layer 4 SPyr to layer 2/3 MarSOM
P4to2{5, 2} = para4E2layer1; % layer 4 SSte to layer 2/3 MarSOM
P4to2{5, 3} = para4FSto4E; % layer 4 FsPV to layer 2/3 MarSOM
P4to2{5, 4} = para4RSto4E; % layer 4 RsNP to layer 2/3 MarSOM

P4to2{6, 1} = para4E2layer1; % layer 4 SPyr to layer 2/3 DBCSOM
P4to2{6, 2} = para4E2layer1; % layer 4 SSte to layer 2/3 DBCSOM
P4to2{6, 3} = para4FSto4E; % layer 4 FsPV to layer 2/3 DBCSOM
P4to2{6, 4} = para4RSto4E; % layer 4 RsNP to layer 2/3 DBCSOM

P4to2{7, 1} = para4E2trans; % layer 4 SPyr to layer 2/3 BipVIP
P4to2{7, 2} = para4E2trans; % layer 4 SSte to layer 2/3 BipVIP
P4to2{7, 3} = para4FSto4E; % layer 4 FsPV to layer 2/3 BipVIP
P4to2{7, 4} = para4RSto4E; % layer 4 RsNP to layer 2/3 BipVIP

P4to2{8, 1} = para4E2trans; % layer 4 SPyr to layer 2/3 BipVIP
P4to2{8, 2} = para4E2trans; % layer 4 SSte to layer 2/3 BipVIP
P4to2{8, 3} = para4FSto4E; % layer 4 FsPV to layer 2/3 BipVIP
P4to2{8, 4} = para4RSto4E; % layer 4 RsNP to layer 2/3 BipVIP

P4to2{9, 1} = para4E2trans; % layer 4 SPyr to layer 2/3 BipCR
P4to2{9, 2} = para4E2trans; % layer 4 SSte to layer 2/3 BipCR
P4to2{9, 3} = para4FSto4E; % layer 4 FsPV to layer 2/3 BipCR
P4to2{9, 4} = para4RSto4E; % layer 4 RsNP to layer 2/3 BipCR

P4to2{10, 1} = para4E2lateral; % layer 4 SPyr to layer 2/3 SbcCR
P4to2{10, 2} = para4E2lateral; % layer 4 SSte to layer 2/3 SbcCR
P4to2{10, 3} = para4FSto4E; % layer 4 FsPV to layer 2/3 SbcCR
P4to2{10, 4} = para4RSto4E; % layer 4 RsNP to layer 2/3 SbcCR

P4to2{11, 1} = para4E2local2; % layer 4 SPyr to layer 2/3 NG
P4to2{11, 2} = para4E2local2; % layer 4 SSte to layer 2/3 NG
P4to2{11, 3} = para4FSto4E; % layer 4 FsPV to layer 2/3 NG
P4to2{11, 4} = para4RSto4E; % layer 4 RsNP to layer 2/3 NG