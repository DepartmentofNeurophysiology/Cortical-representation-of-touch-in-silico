function P2to2 = Synapse_parameter_L23toL23_func
% parameters are saves as mean+-std; except for delay, which in the form of
% [velocity, release delay]; and the delay between a connection is
% determined as distance/velocity + release delay
% currently all distance is in microns, thus volecity is micron/ms

%% E to E
para23E2E.I0 = [0.9,0.6];      %Feldmeyer 2006
para23E2E.tau1 = [15.7,4.5];
para23E2E.tau2 = [0.7,0.2];
para23E2E.Plas = [0.7,0.3];
para23E2E.Fail = [0.03,0.02];
para23E2E.STD.U = [0.7,0.1];
para23E2E.STD.F = [1000,100];
para23E2E.STD.D = [800,80];
para23E2E.Nb = 3;
para23E2E.va = 0.18;
para23E2E.delay = [190, 1.1];

%% Bip VIP
% data missing (VIP+/CR-)
para23E2Bip_VIP.I0 = [0.42,0.32];           %adopted from Rozov et al 2001
para23E2Bip_VIP.tau1 = [26,8.1];            %in motor cortex: Rin = 327 MO, EPSC = 18.7 pA, cv = 0.37, failure = 0.28
para23E2Bip_VIP.tau2 = [0.64,0.27];
para23E2Bip_VIP.Plas = [0.83,0.14];
para23E2Bip_VIP.Fail = [0.28,0.11];
para23E2Bip_VIP.STD.U = [1.0,0.1];
para23E2Bip_VIP.STD.F = [1000,100];
para23E2Bip_VIP.STD.D = [800,80];
para23E2Bip_VIP.Nb = 10;
para23E2Bip_VIP.va = 2;
para23E2Bip_VIP.cv = 0.37;
para23E2Bip_VIP.delay = [190, 0.6];

para23Bip_VIP2E.I0 = [0.49,0.49];            %adopted from Kapfer et al 2007
para23Bip_VIP2E.tau1 = [56.2,6.1];
para23Bip_VIP2E.tau2 = [5.2,2.2];
para23Bip_VIP2E.Plas = [1.02,0.36];
para23Bip_VIP2E.Fail = [0.09,0.05];
para23Bip_VIP2E.STD.U = [0.32,0.1];
para23Bip_VIP2E.STD.F = [62,21];
para23Bip_VIP2E.STD.D = [144,30];
para23Bip_VIP2E.Nb = 5;
para23Bip_VIP2E.va = -0.15;
para23Bip_VIP2E.delay = [190, 0.5];

para23Bip_VIP2Bip_VIP.I0 = [0.49,0.49];    %adopted from kapfer et al 2007  
para23Bip_VIP2Bip_VIP.tau1 = [43.3,6.1];
para23Bip_VIP2Bip_VIP.tau2 = [4.9,2.2];
para23Bip_VIP2Bip_VIP.Plas = [1.02,0.36];
para23Bip_VIP2Bip_VIP.Fail = [0.09,0.05];
para23Bip_VIP2Bip_VIP.STD.U = [0.32,0.1];
para23Bip_VIP2Bip_VIP.STD.F = [62,21];
para23Bip_VIP2Bip_VIP.STD.D = [144,30];
para23Bip_VIP2Bip_VIP.Nb = 5;
para23Bip_VIP2Bip_VIP.va = -0.15;
para23Bip_VIP2Bip_VIP.delay = [190, 0.5];


para23Bip_VIP2BS_PV.I0 = [0.37,0.33];    %adopted from kapfer et al 2007  
para23Bip_VIP2BS_PV.tau1 = [20.2,6.1];
para23Bip_VIP2BS_PV.tau2 = [3.1,2.0];
para23Bip_VIP2BS_PV.Plas = [1.02,0.36];
para23Bip_VIP2BS_PV.Fail = [0.09,0.05];
para23Bip_VIP2BS_PV.STD.U = [0.32,0.1];
para23Bip_VIP2BS_PV.STD.F = [62,21];
para23Bip_VIP2BS_PV.STD.D = [144,30];
para23Bip_VIP2BS_PV.Nb = 5;
para23Bip_VIP2BS_PV.va = -0.15;
para23Bip_VIP2BS_PV.delay = [190, 0.5];


para23Bip_VIP2I.I0 = [0.49,0.49];    %adopted from kapfer et al 2007  
para23Bip_VIP2I.tau1 = [33.3,6.1];
para23Bip_VIP2I.tau2 = [4.9,2.2];
para23Bip_VIP2I.Plas = [1.02,0.36];
para23Bip_VIP2I.Fail = [0.09,0.05];
para23Bip_VIP2I.STD.U = [0.32,0.1];
para23Bip_VIP2I.STD.F = [62,21];
para23Bip_VIP2I.STD.D = [144,30];
para23Bip_VIP2I.Nb = 5;
para23Bip_VIP2I.va = -0.15;
para23Bip_VIP2I.delay = [190, 0.5];

%% DBC VIP

para23E2DBC_VIP.I0 = [0.42,0.32];      % adopted from Holmgren et al 2007
para23E2DBC_VIP.tau1 = [26.0,3.1];
para23E2DBC_VIP.tau2 = [0.6,0.27];
para23E2DBC_VIP.Plas = [0.83,0.13];
para23E2DBC_VIP.Fail = [0.28,0.11];
para23E2DBC_VIP.STD.U = [0.7,0.1];
para23E2DBC_VIP.STD.F = [1000,100];
para23E2DBC_VIP.STD.D = [800,80];
para23E2DBC_VIP.Nb = 10;
para23E2DBC_VIP.va = 2;
para23E2DBC_VIP.cv = 0.37;
para23E2DBC_VIP.delay = [190, 0.6];

para23DBC_VIP2E.I0 = [0.49,0.49];      %adopted from kapfer et al 2007
para23DBC_VIP2E.tau1 = [46.2,4.1];
para23DBC_VIP2E.tau2 = [5.4,2.2];
para23DBC_VIP2E.Plas = [1.4,0.5];
para23DBC_VIP2E.Fail = [0.28,0.14];
para23DBC_VIP2E.STD.U = [0.32,0.1];
para23DBC_VIP2E.STD.F = [62,21];
para23DBC_VIP2E.STD.D = [144,30];
para23DBC_VIP2E.Nb = 5;
para23DBC_VIP2E.va = -0.6;
para23DBC_VIP2E.delay = [190, 0.5];

para23DBC_VIP2DBC_VIP.I0 = [0.49,0.56];      %adopted from kapfer et al 2007
para23DBC_VIP2DBC_VIP.tau1 = [33.3,6];
para23DBC_VIP2DBC_VIP.tau2 = [4.9,1.4];
para23DBC_VIP2DBC_VIP.Plas = [1.4,0.2];
para23DBC_VIP2DBC_VIP.Fail = [0.28,0.14];
para23DBC_VIP2DBC_VIP.STD.U = [0.32,0.1];
para23DBC_VIP2DBC_VIP.STD.F = [62,21];
para23DBC_VIP2DBC_VIP.STD.D = [144,30];
para23DBC_VIP2DBC_VIP.Nb = 5;
para23DBC_VIP2DBC_VIP.va = -0.6;
para23DBC_VIP2DBC_VIP.delay = [190, 0.5];

para23DBC_VIP2BS_PV.I0 = [0.37,0.33];      %adopted from kapfer et al 2007
para23DBC_VIP2BS_PV.tau1 = [20.3,6];
para23DBC_VIP2BS_PV.tau2 = [3.1,1.0];
para23DBC_VIP2BS_PV.Plas = [1.4,0.2];
para23DBC_VIP2BS_PV.Fail = [0.28,0.14];
para23DBC_VIP2BS_PV.STD.U = [0.32,0.1];
para23DBC_VIP2BS_PV.STD.F = [62,21];
para23DBC_VIP2BS_PV.STD.D = [144,30];
para23DBC_VIP2BS_PV.Nb = 4;
para23DBC_VIP2BS_PV.va = -0.6;
para23DBC_VIP2BS_PV.delay = [190, 0.5];


para23DBC_VIP2I.I0 = [0.49,0.56];    %adopted from kapfer et al 2007  
para23DBC_VIP2I.tau1 = [33.3,6];
para23DBC_VIP2I.tau2 = [4.9,1.4];
para23DBC_VIP2I.Plas = [1.42,0.36];
para23DBC_VIP2I.Fail = [0.27,0.14];
para23DBC_VIP2I.STD.U = [0.32,0.1];
para23DBC_VIP2I.STD.F = [62,21];
para23DBC_VIP2I.STD.D = [144,30];
para23DBC_VIP2I.Nb = 5;
para23DBC_VIP2I.va = -0.6;
para23DBC_VIP2I.delay = [190, 0.5];

%% Sbc CR
%CR+/VIP-
para23E2SBC_CR.I0 = [0.42,0.36];           % adopted from Holmgren et al 2003
para23E2SBC_CR.tau1 = [10.25,1.58];       %cv and failure from Koester and Johnston science 2005
para23E2SBC_CR.tau2 = [1.26,0.5];
para23E2SBC_CR.Plas = [1.6,0.43];
para23E2SBC_CR.Fail = [0.40,0.1];
para23E2SBC_CR.STD.U = [0.16,0.1];
para23E2SBC_CR.STD.F = [376,130];
para23E2SBC_CR.STD.D = [45,11];
para23E2SBC_CR.Nb = 3;
para23E2SBC_CR.va = 0.32;
para23E2SBC_CR.delay = [190, 0.5];

para23SBC_CR2E.I0 = [0.49,0.49];      %adopted from kapfer et al 2007
para23SBC_CR2E.tau1 = [56.2,12.1];       %cv and failure rate is from L4 RSNP neurons; sun et al 2006
para23SBC_CR2E.tau2 = [5.4,1.2];
para23SBC_CR2E.Plas = [0.7,0.3];
para23SBC_CR2E.STD.U = [0.37,0.1];
para23SBC_CR2E.STD.F = [62,21];
para23SBC_CR2E.STD.D = [144,48];
para23SBC_CR2E.Nb = 5;
para23SBC_CR2E.va = -0.44;
para23SBC_CR2E.delay = [190, 0.5];

para23SBC_CR2SBC_CR.I0 = [0.49,0.56];      %adopted from kapfer et al 2007
para23SBC_CR2SBC_CR.tau1 = [33,6];
para23SBC_CR2SBC_CR.tau2 = [4.9,2.1];
para23SBC_CR2SBC_CR.Plas = [0.74,0.3];
para23SBC_CR2SBC_CR.STD.U = [0.32,0.1];
para23SBC_CR2SBC_CR.STD.F = [62,21];
para23SBC_CR2SBC_CR.STD.D = [144,48];
para23SBC_CR2SBC_CR.Nb = 5;
para23SBC_CR2SBC_CR.va = -0.44;
para23SBC_CR2SBC_CR.delay = [190, 0.5];

para23SBC_CR2BS_PV.I0 = [0.37,0.33];      %adopted from kapfer et al 2007
para23SBC_CR2BS_PV.tau1 = [20,6];
para23SBC_CR2BS_PV.tau2 = [3.1,1.1];
para23SBC_CR2BS_PV.Plas = [1.4,0.3];
para23SBC_CR2BS_PV.STD.U = [0.32,0.1];
para23SBC_CR2BS_PV.STD.F = [62,21];
para23SBC_CR2BS_PV.STD.D = [144,48];
para23SBC_CR2BS_PV.Nb = 4;
para23SBC_CR2BS_PV.va = -0.44;
para23SBC_CR2BS_PV.delay = [190, 0.5];

para23SBC_CR2I.I0 = [0.49,0.56];    %adopted from kapfer et al 2007  
para23SBC_CR2I.tau1 = [33.3,6];
para23SBC_CR2I.tau2 = [4.9,2.1];
para23SBC_CR2I.Plas = [1.1,0.36];
para23SBC_CR2I.STD.U = [0.32,0.1];
para23SBC_CR2I.STD.F = [62,21];
para23SBC_CR2I.STD.D = [144,48];
para23SBC_CR2I.Nb = 5;
para23SBC_CR2I.va = -0.44;
para23SBC_CR2I.delay = [190, 0.5];

%% Bip CR/VIP+
para23E2Bip_CR.I0 = [0.46,0.41];           % adopted from Holmgren et al 2003
para23E2Bip_CR.tau1 = [26.25,8.1];      
para23E2Bip_CR.tau2 = [0.6,0.3];
para23E2Bip_CR.Plas = [0.83,0.14];
para23E2Bip_CR.STD.U = [0.7,0.1];
para23E2Bip_CR.STD.F = [1000,100];
para23E2Bip_CR.STD.D = [800,80];
para23E2Bip_CR.Nb = 10;
para23E2Bip_CR.va = 2;
para23E2Bip_CR.cv = 0.37;
para23E2Bip_CR.delay = [190, 0.6];

para23Bip_CR2E.I0 = [0.49,0.49];      %adopted from kapfer et al 2007
para23Bip_CR2E.tau1 = [56.2,12.1];
para23Bip_CR2E.tau2 = [5.4,2.2];
para23Bip_CR2E.Plas = [1.4,0.5];
para23Bip_CR2E.STD.U = [0.22,0.05];
para23Bip_CR2E.STD.F = [376,110];
para23Bip_CR2E.STD.D = [45,11];
para23Bip_CR2E.Nb = 5;
para23Bip_CR2E.va = -0.6;
para23Bip_CR2E.delay = [190, 0.5];

para23Bip_CR2Bip_CR.I0 = [0.49,0.56];      %adopted from kapfer et al 2007
para23Bip_CR2Bip_CR.tau1 = [33.3,6.0];
para23Bip_CR2Bip_CR.tau2 = [4.9,2.3];
para23Bip_CR2Bip_CR.Plas = [1.33,0.3];
para23Bip_CR2Bip_CR.STD.U = [0.22,0.05];
para23Bip_CR2Bip_CR.STD.F = [376,110];
para23Bip_CR2Bip_CR.STD.D = [45,11];
para23Bip_CR2Bip_CR.Nb = 5;
para23Bip_CR2Bip_CR.va = -0.6;
para23Bip_CR2Bip_CR.delay = [190, 0.5];

para23Bip_CR2I.I0 = [0.49,0.56];    %adopted from kapfer et al 2007  
para23Bip_CR2I.tau1 = [33.3,6];
para23Bip_CR2I.tau2 = [4.6,2.3];
para23Bip_CR2I.Plas = [1.8,0.3];
para23Bip_CR2I.STD.U = [0.22,0.05];
para23Bip_CR2I.STD.F = [376,110];
para23Bip_CR2I.STD.D = [45,11];
para23Bip_CR2I.Nb = 5;
para23Bip_CR2I.va = -0.6;
para23Bip_CR2I.delay = [190, 0.5];

para23Bip_CR2BS_PV.I0 = [0.37,0.33];    %adopted from kapfer et al 2007  
para23Bip_CR2BS_PV.tau1 = [20.3,6];
para23Bip_CR2BS_PV.tau2 = [3.1,1.2];
para23Bip_CR2BS_PV.Plas = [1.42,0.2];
para23Bip_CR2BS_PV.STD.U = [0.22,0.05];
para23Bip_CR2BS_PV.STD.F = [376,110];
para23Bip_CR2BS_PV.STD.D = [45,11];
para23Bip_CR2BS_PV.Nb = 5;
para23Bip_CR2BS_PV.va = -0.6;
para23Bip_CR2BS_PV.delay = [190, 0.5];


%% Fs PV

% para23E2BS_PV.I0 = [3.48,2.50];                     %Holmgren et al 2003
% para23E2BS_PV.tau1 = [19,4.1];                      %cv and failure from Koester and Johnston science 2005
% para23E2BS_PV.tau2 = [2.1,1.1];
para23E2BS_PV.I0 = [0.82,0.6];                     %Averman et al 2011
para23E2BS_PV.tau1 = [14,1.6];                      %cv and failure from Koester and Johnston science 2005
para23E2BS_PV.tau2 = [1.6,0.6];
para23E2BS_PV.Plas = [0.7,0.13];
para23E2BS_PV.STD.U = [0.7,0.1];
para23E2BS_PV.STD.F = [1000,100];
para23E2BS_PV.STD.D = [800,80];
para23E2BS_PV.Nb = 20;
para23E2BS_PV.va = 0.7;
para23E2BS_PV.delay = [190, 0.3];

% para23BS_PV2E.I0 = [2.9,1.1];                     
% para23BS_PV2E.tau1 = [44,11];
% para23BS_PV2E.tau2 = [2.1,1.1];
para23BS_PV2E.I0 = [0.56,0.54];                     
para23BS_PV2E.tau1 = [40,9];
para23BS_PV2E.tau2 = [1.4,0.4];
para23BS_PV2E.Plas = [0.7,0.2];
para23BS_PV2E.STD.U = [0.35,0.1];
para23BS_PV2E.STD.F = [21,9];
para23BS_PV2E.STD.D = [760,200];
para23BS_PV2E.Nb = 7;
para23BS_PV2E.va = -0.5;
para23BS_PV2E.delay = [190, 0.1];


% para23BS_PV2BS_PV.I0 = [2.96,2.52];               
% para23BS_PV2BS_PV.tau1 = [44,11];
% para23BS_PV2BS_PV.tau2 = [2,1.1];
para23BS_PV2BS_PV.I0 = [0.56,0.5];               
para23BS_PV2BS_PV.tau1 = [16,3];
para23BS_PV2BS_PV.tau2 = [1.8,0.6];
para23BS_PV2BS_PV.Plas = [0.7,0.15];
para23BS_PV2BS_PV.STD.U = [0.35,0.1];
para23BS_PV2BS_PV.STD.F = [21,9];
para23BS_PV2BS_PV.STD.D = [760,200];
para23BS_PV2BS_PV.Nb = 5;
para23BS_PV2BS_PV.va = -0.4;
para23BS_PV2BS_PV.delay = [190, 0.7];

% para23BS_PV2I.I0 = [2.92,2.52];               
% para23BS_PV2I.tau1 = [44,11];
% para23BS_PV2I.tau2 = [2,1.1];
para23BS_PV2I.I0 = [0.83,0.6];               
para23BS_PV2I.tau1 = [36,11];
para23BS_PV2I.tau2 = [2.9,1.1];
para23BS_PV2I.Plas = [0.7,0.15];
para23BS_PV2I.STD.U = [0.35,0.1];
para23BS_PV2I.STD.F = [21,9];
para23BS_PV2I.STD.D = [760,200];
para23BS_PV2I.Nb = 5;
para23BS_PV2I.va = -0.4;
para23BS_PV2I.delay = [190, 0.3];

%% CH

para23E2CH_PV.I0 = [0.82,0.60];                   % adopted from Holmgren et al 2003
para23E2CH_PV.tau1 = [18.25,2.58];                %cv and failure from Koester and Johnston science 2005
para23E2CH_PV.tau2 = [0.6,0.3];
para23E2CH_PV.Plas = [0.7,0.13];
para23E2CH_PV.STD.U = [0.7,0.2];
para23E2CH_PV.STD.F = [1000,100];
para23E2CH_PV.STD.D = [800,80];
para23E2CH_PV.Nb = 20;
para23E2CH_PV.va = 0.7;
para23E2CH_PV.delay = [190, 0.3];

para23CH_PV2E.I0 = [0.75,0.4];      
para23CH_PV2E.tau1 = [15,4];
para23CH_PV2E.tau2 = [1.2,0.4];
para23CH_PV2E.Plas = [0.7,0.15];
para23CH_PV2E.STD.U = [0.35,0.1];
para23CH_PV2E.STD.F = [21,9];
para23CH_PV2E.STD.D = [760,200];
para23CH_PV2E.Nb = 7;
para23CH_PV2E.va = -0.5;
para23CH_PV2E.delay = [190, 0.1];

para23CH_PV2CH_PV.I0 = [0,0];      
para23CH_PV2CH_PV.tau1 = [16,4];
para23CH_PV2CH_PV.tau2 = [2,1.1];
para23CH_PV2CH_PV.Plas = [0.7,0.15];
para23CH_PV2CH_PV.STD.U = [0.35,0.1];
para23CH_PV2CH_PV.STD.F = [21,9];
para23CH_PV2CH_PV.STD.D = [760,200];
para23CH_PV2CH_PV.Nb = 0;
para23CH_PV2CH_PV.va = 0;
para23CH_PV2CH_PV.delay = [190, 0.7];

para23CH_PV2I.I0 = [0,0];      
para23CH_PV2I.tau1 = [16,4];
para23CH_PV2I.tau2 = [2,1.1];
para23CH_PV2I.Plas = [0.7,0.15];
para23CH_PV2I.STD.U = [0.35,0.1];
para23CH_PV2I.STD.F = [21,9];
para23CH_PV2I.STD.D = [760,200];
para23CH_PV2I.Nb = 0;
para23CH_PV2I.va = 0;
para23CH_PV2I.delay = [190, 0.3];

%% Bursting PV
para23E2Res_PV.I0 = [0.38,0.25];         %Blatow et al 2003
para23E2Res_PV.tau1 = [19.25,2.2];     %cv and failure from Koester and Johnston science 2005
para23E2Res_PV.tau2 = [2.76,1.0];
para23E2Res_PV.Plas = [0.51,0.13];
para23E2Res_PV.STD.U = [0.5,0.2];
para23E2Res_PV.STD.F = [1000,100];
para23E2Res_PV.STD.D = [800,80];
para23E2Res_PV.Nb = 5;
para23E2Res_PV.va = 0.2;
para23E2Res_PV.delay = [190, 0.5];

para23Res_PV2E.I0 = [0.67,0.5];      
para23Res_PV2E.tau1 = [22.25,3.58];
para23Res_PV2E.tau2 = [2.1,1.0];
para23Res_PV2E.Plas = [1.27,0.6];
para23Res_PV2E.Fail = [0.09,0.12];
para23Res_PV2E.STD.U = [0.22,0.05];
para23Res_PV2E.STD.F = [376,110];
para23Res_PV2E.STD.D = [45,11];
para23Res_PV2E.Nb = 10;
para23Res_PV2E.va = -0.25;
para23Res_PV2E.delay = [190, 0.5];

para23Res_PV2BS_PV.I0 = [0.67,0.42];      
para23Res_PV2BS_PV.tau1 = [22.25,3.58];
para23Res_PV2BS_PV.tau2 = [2.1,1.0];
para23Res_PV2BS_PV.Plas = [0.86,0.20];
para23Res_PV2BS_PV.STD.U = [0.15,0.05];
para23Res_PV2BS_PV.STD.F = [376,110];
para23Res_PV2BS_PV.STD.D = [45,11];
para23Res_PV2BS_PV.Nb = 5;
para23Res_PV2BS_PV.va = -0.15;
para23Res_PV2BS_PV.delay = [190, 0.5];

para23Res_PV2I.I0 = [0.81,0.72];      
para23Res_PV2I.tau1 = [22.6,3.58];
para23Res_PV2I.tau2 = [2.1,1.0];
para23Res_PV2I.Plas = [1.53,0.63];
para23Res_PV2I.STD.U = [0.15,0.05];
para23Res_PV2I.STD.F = [376,110];
para23Res_PV2I.STD.D = [45,11];
para23Res_PV2I.Nb = 10;
para23Res_PV2I.va = -0.25;
para23Res_PV2I.delay = [190, 0.5];

%% Mar SOM
para23E2Mar_SOM.I0 = [0.25,0.2];               %Kapfer et al 2007
para23E2Mar_SOM.tau1 = [19.25,2.2];           %cv and failure rate from Bartley et al J Neurophysiology 2008
para23E2Mar_SOM.tau2 = [2.8,1.0];
para23E2Mar_SOM.Plas = [1.91,0.8];
para23E2Mar_SOM.STD.U = [0.20,0.05];
para23E2Mar_SOM.STD.F = [376,110];
para23E2Mar_SOM.STD.D = [45,11];
para23E2Mar_SOM.Nb = 10;
para23E2Mar_SOM.va = 1;
para23E2Mar_SOM.cv = 0.55;
para23E2Mar_SOM.delay = [190, 0.6];

para23Mar_SOM2E.I0 = [0.56,0.4];      
para23Mar_SOM2E.tau1 = [16.1,3.9];
para23Mar_SOM2E.tau2 = [0.7,0.3];
para23Mar_SOM2E.Plas = [0.8,0.13];
para23Mar_SOM2E.STD.U = [0.25,0.1];
para23Mar_SOM2E.STD.F = [21,9];
para23Mar_SOM2E.STD.D = [760,200];
para23Mar_SOM2E.Nb = 3;
para23Mar_SOM2E.va = -1.5;
para23Mar_SOM2E.delay = [190, 0.5];

para23Mar_SOM2Mar_SOM.I0 = [0.56,0.4];      
para23Mar_SOM2Mar_SOM.tau1 = [16.1,3.9];
para23Mar_SOM2Mar_SOM.tau2 = [0.7,0.3];
para23Mar_SOM2Mar_SOM.Plas = [0.8,0.3];
para23Mar_SOM2Mar_SOM.STD.U = [0.25,0.1];
para23Mar_SOM2Mar_SOM.STD.F = [21,9];
para23Mar_SOM2Mar_SOM.STD.D = [760,200];
para23Mar_SOM2Mar_SOM.Nb = 3;
para23Mar_SOM2Mar_SOM.va = -1.5;
para23Mar_SOM2Mar_SOM.delay = [190, 0.5];

para23Mar_SOM2I.I0 = [0.56,0.4];      
para23Mar_SOM2I.tau1 = [16.1,3.9];
para23Mar_SOM2I.tau2 = [0.7,0.3];
para23Mar_SOM2I.Plas = [0.8,0.3];
para23Mar_SOM2I.STD.U = [0.25,0.1];
para23Mar_SOM2I.STD.F = [21,9];
para23Mar_SOM2I.STD.D = [760,200];
para23Mar_SOM2I.Nb = 3;
para23Mar_SOM2I.va = -1.5;
para23Mar_SOM2I.delay = [190, 0.5];

%% NG

para23E2NGC_AC.I0 = [0.59,0.41];      %Wozny and Willianms 2011
para23E2NGC_AC.tau1 = [9.25,3.3];  %cv and failure from Koester and Johnston science 2005
para23E2NGC_AC.tau2 = [1.26,0.53];
para23E2NGC_AC.Plas = [0.83,0.14];
para23E2NGC_AC.STD.U = [0.4,0.1];
para23E2NGC_AC.STD.F = [1000,100];
para23E2NGC_AC.STD.D = [800,80];
para23E2NGC_AC.Nb = 7;
para23E2NGC_AC.va = 1.1;
para23E2NGC_AC.delay = [190, 0.6];

para23NGC_AC2E.I0 = [0.58,0.4];      
para23NGC_AC2E.tau1 = [100,19];
para23NGC_AC2E.tau2 = [53,6];
para23NGC_AC2E.Plas = [0.51,0.13];
para23NGC_AC2E.STD.U = [0.25,0.1];
para23NGC_AC2E.STD.F = [21,9];
para23NGC_AC2E.STD.D = [760,200];
para23NGC_AC2E.Nb = 5;
para23NGC_AC2E.va = -0.1;
para23NGC_AC2E.delay = [190, 0.5];

para23NGC_AC2NGC_AC.I0 = [0.58,0.4];      
para23NGC_AC2NGC_AC.tau1 = [100,19];
para23NGC_AC2NGC_AC.tau2 = [53,6];
para23NGC_AC2NGC_AC.Plas = [0.51,0.13];
para23NGC_AC2NGC_AC.STD.U = [0.25,0.1];
para23NGC_AC2NGC_AC.STD.F = [21,9];
para23NGC_AC2NGC_AC.STD.D = [760,200];
para23NGC_AC2NGC_AC.Nb = 5;
para23NGC_AC2NGC_AC.va = -0.1;
para23NGC_AC2NGC_AC.delay = [190, 0.5];

para23NGC_AC2I.I0 = [0.58,0.4];      
para23NGC_AC2I.tau1 = [100,19];
para23NGC_AC2I.tau2 = [53,6];
para23NGC_AC2I.Plas = [0.51,0.13];
para23NGC_AC2I.STD.U = [0.25,0.1];
para23NGC_AC2I.STD.F = [21,9];
para23NGC_AC2I.STD.D = [760,200];
para23NGC_AC2I.Nb = 5;
para23NGC_AC2I.va = -0.1;
para23NGC_AC2I.delay = [190, 0.5];

%% Bit SOM
para23E2Bit_SOM.I0 = [0.25,0.2];      %Kapfer et al 2007
para23E2Bit_SOM.tau1 = [19.25,2.3];  %cv and failure rate from Bartley et al J Neurophysiology 2008
para23E2Bit_SOM.tau2 = [2.8,1.0];
para23E2Bit_SOM.Plas = [1.9,0.8];
para23E2Bit_SOM.STD.U = [0.20,0.05];
para23E2Bit_SOM.STD.F = [376,110];
para23E2Bit_SOM.STD.D = [45,11];
para23E2Bit_SOM.Nb = 10;
para23E2Bit_SOM.va = 1;
para23E2Bit_SOM.cv = 0.45;
para23E2Bit_SOM.delay = [190, 0.6];

para23Bit_SOM2E.I0 = [0.78,0.53];      
para23Bit_SOM2E.tau1 = [16.1,3.9];
para23Bit_SOM2E.tau2 = [0.7,0.3];
para23Bit_SOM2E.Plas = [0.8,0.13];
para23Bit_SOM2E.STD.U = [0.25,0.1];
para23Bit_SOM2E.STD.F = [21,9];
para23Bit_SOM2E.STD.D = [760,200];
para23Bit_SOM2E.Nb = 3;
para23Bit_SOM2E.va = -1.5;
para23Bit_SOM2E.delay = [190, 0.5];

para23Bit_SOM2Bit_SOM.I0 = [0.78,0.53];      
para23Bit_SOM2Bit_SOM.tau1 = [16.1,3.9];
para23Bit_SOM2Bit_SOM.tau2 = [0.7,0.3];
para23Bit_SOM2Bit_SOM.Plas = [0.8,0.13];
para23Bit_SOM2Bit_SOM.STD.U = [0.25,0.1];
para23Bit_SOM2Bit_SOM.STD.F = [21,9];
para23Bit_SOM2Bit_SOM.STD.D = [760,200];
para23Bit_SOM2Bit_SOM.Nb = 3;
para23Bit_SOM2Bit_SOM.va = -1.5;
para23Bit_SOM2Bit_SOM.delay = [190, 0.5];

para23Bit_SOM2I.I0 = [0.78,0.53];      
para23Bit_SOM2I.tau1 = [16.1,3.9];
para23Bit_SOM2I.tau2 = [0.7,0.3];
para23Bit_SOM2I.Plas = [0.8,0.3];
para23Bit_SOM2I.STD.U = [0.25,0.1];
para23Bit_SOM2I.STD.F = [21,9];
para23Bit_SOM2I.STD.D = [760,200];
para23Bit_SOM2I.Nb = 3;
para23Bit_SOM2I.va = -1.5;
para23Bit_SOM2I.delay = [190, 0.5];
% ------------------------------------------------------------------------
% cross talk between inhibitory neurons (most of the data is missing)
% here use the data from Holmgren et al 2003

% % ------------------------------------------------------------------------
% % here is what we need:
% 
% pare23DBC_VIP2Bip_VIP
% pare23DBC_VIP2SBC_CR
% pare23DBC_VIP2BS_PV
% pare23DBC_VIP2CH_PV
% pare23DBC_VIP2Res_PV
% pare23DBC_VIP2Mar_SOM
% pare23DBC_VIP2NGC_AC
% pare23DBC_VIP2Bit_SOM
% 
% pare23Bip_VIP2DBC_VIP
% pare23Bip_VIP2SBC_CR
% pare23Bip_VIP2BS_PV
% pare23Bip_VIP2CH_PV
% pare23Bip_VIP2Res_PV
% pare23Bip_VIP2Mar_SOM
% pare23Bip_VIP2NGC_AC
% pare23Bip_VIP2Bit_SOM
% pare23Bip_VIP2DBC_VIP
% ... ...
% %
% ------------------------------------------------------------------------

P2E2E={para23E2E};
P2E2I={para23E2BS_PV,para23E2CH_PV,para23E2Res_PV,para23E2Mar_SOM,para23E2Bit_SOM,para23E2DBC_VIP,para23E2Bip_VIP,para23E2Bip_CR,para23E2SBC_CR,para23E2NGC_AC};
P2I2E={para23BS_PV2E,para23CH_PV2E,para23Res_PV2E,para23Mar_SOM2E,para23Bit_SOM2E,para23DBC_VIP2E,para23Bip_VIP2E,para23Bip_CR2E,para23SBC_CR2E,para23NGC_AC2E};
P2I2DBC_VIP={para23BS_PV2I,para23CH_PV2I,para23Res_PV2I,para23Mar_SOM2I,para23Bit_SOM2I,para23DBC_VIP2DBC_VIP,para23Bip_VIP2I,para23Bip_CR2I,para23SBC_CR2I,para23NGC_AC2I};
P2I2Bip_VIP={para23BS_PV2I,para23CH_PV2I,para23Res_PV2I,para23Mar_SOM2I,para23Bit_SOM2I,para23DBC_VIP2I,para23Bip_VIP2Bip_VIP,para23Bip_CR2I,para23SBC_CR2I,para23NGC_AC2I};
P2I2SBC_CR={para23BS_PV2I,para23CH_PV2I,para23Res_PV2I,para23Mar_SOM2I,para23Bit_SOM2I,para23DBC_VIP2I,para23Bip_VIP2I,para23Bip_CR2I,para23SBC_CR2SBC_CR,para23NGC_AC2I};
P2I2BS_PV={para23BS_PV2BS_PV,para23CH_PV2I,para23Res_PV2BS_PV,para23Mar_SOM2I,para23Bit_SOM2I,para23DBC_VIP2BS_PV,para23Bip_VIP2BS_PV,para23Bip_CR2BS_PV,para23SBC_CR2BS_PV,para23NGC_AC2I};
P2I2CH_PV={para23BS_PV2I,para23CH_PV2CH_PV,para23Res_PV2I,para23Mar_SOM2I,para23Bit_SOM2I,para23DBC_VIP2BS_PV,para23Bip_VIP2BS_PV,para23Bip_CR2BS_PV,para23SBC_CR2BS_PV,para23NGC_AC2I};
P2I2Res_PV={para23BS_PV2I,para23CH_PV2I,para23Res_PV2I,para23Mar_SOM2I,para23Bit_SOM2I,para23DBC_VIP2I,para23Bip_VIP2I,para23Bip_CR2I,para23SBC_CR2I,para23NGC_AC2I};
P2I2Mar_SOM={para23BS_PV2I,para23CH_PV2I,para23Res_PV2I,para23Mar_SOM2Mar_SOM,para23Bit_SOM2I,para23DBC_VIP2I,para23Bip_VIP2I,para23Bip_CR2I,para23SBC_CR2I,para23NGC_AC2I};
P2I2NGC_AC={para23BS_PV2I,para23CH_PV2I,para23Res_PV2I,para23Mar_SOM2I,para23Bit_SOM2I,para23DBC_VIP2I,para23Bip_VIP2I,para23Bip_CR2I,para23SBC_CR2I,para23NGC_AC2NGC_AC};
P2I2Bit_SOM={para23BS_PV2I,para23CH_PV2I,para23Res_PV2I,para23Mar_SOM2I,para23Bit_SOM2Bit_SOM,para23DBC_VIP2I,para23Bip_VIP2I,para23Bip_CR2I,para23SBC_CR2I,para23NGC_AC2I};
P2I2Bip_CR={para23BS_PV2I,para23CH_PV2I,para23Res_PV2I,para23Mar_SOM2I,para23Bit_SOM2I,para23DBC_VIP2I,para23Bip_VIP2I,para23Bip_CR2Bip_CR,para23SBC_CR2I,para23NGC_AC2I};

P2I2I={P2I2BS_PV,P2I2CH_PV,P2I2Res_PV,P2I2Mar_SOM,P2I2Bit_SOM,P2I2DBC_VIP,P2I2Bip_VIP,P2I2Bip_CR,P2I2SBC_CR,P2I2NGC_AC};

P2to2 = {};

P2to2{1, 1} = para23E2E;       %layer 2/3 E to layer 2/3 E
P2to2{1, 2} = para23BS_PV2E;   %layer 2/3 FsPV to layer 2/3 E
P2to2{1, 3} = para23CH_PV2E;   %layer 2/3 ChPV to layer 2/3 E
P2to2{1, 4} = para23Res_PV2E;  %layer 2/3 ResPV to layer 2/3 E
P2to2{1, 5} = para23Mar_SOM2E; %layer 2/3 MarSOM to layer 2/3 E
P2to2{1, 6} = para23Mar_SOM2E; %layer 2/3 BitSOM to layer 2/3 E
P2to2{1, 7} = para23Bip_VIP2E; %layer 2/3 BipVIP to layer 2/3 E
P2to2{1, 8} = para23Bip_VIP2E; %layer 2/3 BipVIP to layer 2/3 E
P2to2{1, 9} = para23Bip_CR2E;  %layer 2/3 BipCR to layer 2/3 E
P2to2{1, 10} = para23SBC_CR2E; %layer 2/3 SbcCR to layer 2/3 E
P2to2{1, 11} = para23NGC_AC2E; %layer 2/3 NGC to layer 2/3 E

P2to2{2, 1} = para23E2BS_PV;       %layer 2/3 E to layer 2/3 FsPV
P2to2{2, 2} = para23BS_PV2BS_PV;   %layer 2/3 FsPV to layer 2/3 FsPV
P2to2{2, 3} = para23BS_PV2I;   %layer 2/3 ChPV to layer 2/3 FsPV
P2to2{2, 4} = para23Res_PV2BS_PV;  %layer 2/3 ResPV to layer 2/3 FsPV
P2to2{2, 5} = para23Mar_SOM2I; %layer 2/3 MarSOM to layer 2/3 FsPV
P2to2{2, 6} = para23Mar_SOM2I; %layer 2/3 BitSOM to layer 2/3 FsPV
P2to2{2, 7} = para23Bip_VIP2BS_PV; %layer 2/3 BipVIP to layer 2/3 FsPV
P2to2{2, 8} = para23Bip_VIP2BS_PV; %layer 2/3 BipVIP to layer 2/3 FsPV
P2to2{2, 9} = para23Bip_CR2BS_PV;  %layer 2/3 BipCR to layer 2/3 FsPV
P2to2{2, 10} = para23SBC_CR2BS_PV; %layer 2/3 SbcCR to layer 2/3 FsPV
P2to2{2, 11} = para23NGC_AC2I; %layer 2/3 NGC to layer 2/3 FsPV

P2to2{3, 1} = para23E2CH_PV;       %layer 2/3 E to layer 2/3 ChPV
P2to2{3, 2} = para23BS_PV2I;   %layer 2/3 FsPV to layer 2/3 ChPV
P2to2{3, 3} = para23CH_PV2CH_PV;   %layer 2/3 ChPV to layer 2/3 ChPV
P2to2{3, 4} = para23Res_PV2I;  %layer 2/3 ResPV to layer 2/3 ChPV
P2to2{3, 5} = para23Mar_SOM2I; %layer 2/3 MarSOM to layer 2/3 ChPV
P2to2{3, 6} = para23Mar_SOM2I; %layer 2/3 BitSOM to layer 2/3 ChPV
P2to2{3, 7} = para23Bip_VIP2BS_PV; %layer 2/3 BipVIP to layer 2/3 ChPV
P2to2{3, 8} = para23Bip_VIP2BS_PV; %layer 2/3 BipVIP to layer 2/3 ChPV
P2to2{3, 9} = para23Bip_CR2BS_PV;  %layer 2/3 BipCR to layer 2/3 ChPV
P2to2{3, 10} = para23SBC_CR2BS_PV; %layer 2/3 SbcCR to layer 2/3 ChPV
P2to2{3, 11} = para23NGC_AC2I; %layer 2/3 NGC to layer 2/3 ChPV

P2to2{4, 1} = para23E2Res_PV;       %layer 2/3 E to layer 2/3 ChPV
P2to2{4, 2} = para23BS_PV2I;   %layer 2/3 FsPV to layer 2/3 ChPV
P2to2{4, 3} = para23CH_PV2I;   %layer 2/3 ChPV to layer 2/3 ChPV
P2to2{4, 4} = para23Res_PV2I;  %layer 2/3 ResPV to layer 2/3 ChPV
P2to2{4, 5} = para23Mar_SOM2I; %layer 2/3 MarSOM to layer 2/3 ChPV
P2to2{4, 6} = para23Mar_SOM2I; %layer 2/3 BitSOM to layer 2/3 ChPV
P2to2{4, 7} = para23Bip_VIP2I; %layer 2/3 BipVIP to layer 2/3 ChPV
P2to2{4, 8} = para23Bip_VIP2I; %layer 2/3 BipVIP to layer 2/3 ChPV
P2to2{4, 9} = para23Bip_CR2I;  %layer 2/3 BipCR to layer 2/3 ChPV
P2to2{4, 10} = para23SBC_CR2I; %layer 2/3 SbcCR to layer 2/3 ChPV
P2to2{4, 11} = para23NGC_AC2I; %layer 2/3 NGC to layer 2/3 ChPV

P2to2{5, 1} = para23E2Mar_SOM;       %layer 2/3 E to layer 2/3 MarSOM
P2to2{5, 2} = para23BS_PV2I;   %layer 2/3 FsPV to layer 2/3 MarSOM
P2to2{5, 3} = para23CH_PV2I;   %layer 2/3 ChPV to layer 2/3 MarSOM
P2to2{5, 4} = para23Res_PV2I;  %layer 2/3 ResPV to layer 2/3 MarSOM
P2to2{5, 5} = para23Mar_SOM2Mar_SOM; %layer 2/3 MarSOM to layer 2/3 MarSOM
P2to2{5, 6} = para23Mar_SOM2Mar_SOM; %layer 2/3 BitSOM to layer 2/3 MarSOM
P2to2{5, 7} = para23Bip_VIP2I; %layer 2/3 BipVIP to layer 2/3 MarSOM
P2to2{5, 8} = para23Bip_VIP2I; %layer 2/3 BipVIP to layer 2/3 MarSOM
P2to2{5, 9} = para23Bip_CR2I;  %layer 2/3 BipCR to layer 2/3 MarSOM
P2to2{5, 10} = para23SBC_CR2I; %layer 2/3 SbcCR to layer 2/3 MarSOM
P2to2{5, 11} = para23NGC_AC2I; %layer 2/3 NGC to layer 2/3 MarSOM

P2to2{6, 1} = para23E2Mar_SOM;       %layer 2/3 E to layer 2/3 MarSOM
P2to2{6, 2} = para23BS_PV2I;   %layer 2/3 FsPV to layer 2/3 MarSOM
P2to2{6, 3} = para23CH_PV2I;   %layer 2/3 ChPV to layer 2/3 MarSOM
P2to2{6, 4} = para23Res_PV2I;  %layer 2/3 ResPV to layer 2/3 MarSOM
P2to2{6, 5} = para23Mar_SOM2Mar_SOM; %layer 2/3 MarSOM to layer 2/3 MarSOM
P2to2{6, 6} = para23Mar_SOM2Mar_SOM; %layer 2/3 BitSOM to layer 2/3 MarSOM
P2to2{6, 7} = para23Bip_VIP2I; %layer 2/3 BipVIP to layer 2/3 MarSOM
P2to2{6, 8} = para23Bip_VIP2I; %layer 2/3 BipVIP to layer 2/3 MarSOM
P2to2{6, 9} = para23Bip_CR2I;  %layer 2/3 BipCR to layer 2/3 MarSOM
P2to2{6, 10} = para23SBC_CR2I; %layer 2/3 SbcCR to layer 2/3 MarSOM
P2to2{6, 11} = para23NGC_AC2I; %layer 2/3 NGC to layer 2/3 MarSOM

P2to2{7, 1} = para23E2Bip_VIP;       %layer 2/3 E to layer 2/3 BipVIP
P2to2{7, 2} = para23BS_PV2I;   %layer 2/3 FsPV to layer 2/3 BipVIP
P2to2{7, 3} = para23CH_PV2I;   %layer 2/3 ChPV to layer 2/3 BipVIP
P2to2{7, 4} = para23Res_PV2I;  %layer 2/3 ResPV to layer 2/3 BipVIP
P2to2{7, 5} = para23Mar_SOM2I; %layer 2/3 MarSOM to layer 2/3 BipVIP
P2to2{7, 6} = para23Mar_SOM2I; %layer 2/3 BitSOM to layer 2/3 BipVIP
P2to2{7, 7} = para23Bip_VIP2Bip_VIP; %layer 2/3 BipVIP to layer 2/3 BipVIP
P2to2{7, 8} = para23Bip_VIP2Bip_VIP; %layer 2/3 BipVIP to layer 2/3 BipVIP
P2to2{7, 9} = para23Bip_CR2I;  %layer 2/3 BipCR to layer 2/3 BipVIP
P2to2{7, 10} = para23SBC_CR2I; %layer 2/3 SbcCR to layer 2/3 BipVIP
P2to2{7, 11} = para23NGC_AC2I; %layer 2/3 NGC to layer 2/3 BipVIP

P2to2{8, 1} = para23E2Bip_VIP;       %layer 2/3 E to layer 2/3 BipVIP
P2to2{8, 2} = para23BS_PV2I;   %layer 2/3 FsPV to layer 2/3 BipVIP
P2to2{8, 3} = para23CH_PV2I;   %layer 2/3 ChPV to layer 2/3 BipVIP
P2to2{8, 4} = para23Res_PV2I;  %layer 2/3 ResPV to layer 2/3 BipVIP
P2to2{8, 5} = para23Mar_SOM2I; %layer 2/3 MarSOM to layer 2/3 BipVIP
P2to2{8, 6} = para23Mar_SOM2I; %layer 2/3 BitSOM to layer 2/3 BipVIP
P2to2{8, 7} = para23Bip_VIP2Bip_VIP; %layer 2/3 BipVIP to layer 2/3 BipVIP
P2to2{8, 8} = para23Bip_VIP2Bip_VIP; %layer 2/3 BipVIP to layer 2/3 BipVIP
P2to2{8, 9} = para23Bip_CR2I;  %layer 2/3 BipCR to layer 2/3 BipVIP
P2to2{8, 10} = para23SBC_CR2I; %layer 2/3 SbcCR to layer 2/3 BipVIP
P2to2{8, 11} = para23NGC_AC2I; %layer 2/3 NGC to layer 2/3 BipVIP

P2to2{9, 1} = para23E2Bip_CR;       %layer 2/3 E to layer 2/3 BipCR
P2to2{9, 2} = para23BS_PV2I;   %layer 2/3 FsPV to layer 2/3 BipCR
P2to2{9, 3} = para23CH_PV2I;   %layer 2/3 ChPV to layer 2/3 BipCR
P2to2{9, 4} = para23Res_PV2I;  %layer 2/3 ResPV to layer 2/3 BipCR
P2to2{9, 5} = para23Mar_SOM2I; %layer 2/3 MarSOM to layer 2/3 BipCR
P2to2{9, 6} = para23Mar_SOM2I; %layer 2/3 BitSOM to layer 2/3 BipCR
P2to2{9, 7} = para23Bip_VIP2I; %layer 2/3 BipVIP to layer 2/3 BipCR
P2to2{9, 8} = para23Bip_VIP2I; %layer 2/3 BipVIP to layer 2/3 BipCR
P2to2{9, 9} = para23Bip_CR2Bip_CR;  %layer 2/3 BipCR to layer 2/3 BipCR
P2to2{9, 10} = para23SBC_CR2I; %layer 2/3 SbcCR to layer 2/3 BipCR
P2to2{9, 11} = para23NGC_AC2I; %layer 2/3 NGC to layer 2/3 BipCR

P2to2{10, 1} = para23E2SBC_CR;       %layer 2/3 E to layer 2/3 SbcCR
P2to2{10, 2} = para23BS_PV2I;   %layer 2/3 FsPV to layer 2/3 SbcCR
P2to2{10, 3} = para23CH_PV2I;   %layer 2/3 ChPV to layer 2/3 SbcCR
P2to2{10, 4} = para23Res_PV2I;  %layer 2/3 ResPV to layer 2/3 SbcCR
P2to2{10, 5} = para23Mar_SOM2I; %layer 2/3 MarSOM to layer 2/3 SbcCR
P2to2{10, 6} = para23Mar_SOM2I; %layer 2/3 BitSOM to layer 2/3 SbcCR
P2to2{10, 7} = para23Bip_VIP2I; %layer 2/3 BipVIP to layer 2/3 SbcCR
P2to2{10, 8} = para23Bip_VIP2I; %layer 2/3 BipVIP to layer 2/3 SbcCR
P2to2{10, 9} = para23Bip_CR2I;  %layer 2/3 BipCR to layer 2/3 SbcCR
P2to2{10, 10} = para23SBC_CR2SBC_CR; %layer 2/3 SbcCR to layer 2/3 SbcCR
P2to2{10, 11} = para23NGC_AC2I; %layer 2/3 NGC to layer 2/3 SbcCR

P2to2{11, 1} = para23E2NGC_AC;       %layer 2/3 E to layer 2/3 SbcCR
P2to2{11, 2} = para23BS_PV2I;   %layer 2/3 FsPV to layer 2/3 SbcCR
P2to2{11, 3} = para23CH_PV2I;   %layer 2/3 ChPV to layer 2/3 SbcCR
P2to2{11, 4} = para23Res_PV2I;  %layer 2/3 ResPV to layer 2/3 SbcCR
P2to2{11, 5} = para23Mar_SOM2I; %layer 2/3 MarSOM to layer 2/3 SbcCR
P2to2{11, 6} = para23Mar_SOM2I; %layer 2/3 BitSOM to layer 2/3 SbcCR
P2to2{11, 7} = para23Bip_VIP2I; %layer 2/3 BipVIP to layer 2/3 SbcCR
P2to2{11, 8} = para23Bip_VIP2I; %layer 2/3 BipVIP to layer 2/3 SbcCR
P2to2{11, 9} = para23Bip_CR2I;  %layer 2/3 BipCR to layer 2/3 SbcCR
P2to2{11, 10} = para23SBC_CR2I; %layer 2/3 SbcCR to layer 2/3 SbcCR
P2to2{11, 11} = para23NGC_AC2NGC_AC; %layer 2/3 NGC to layer 2/3 SbcCR