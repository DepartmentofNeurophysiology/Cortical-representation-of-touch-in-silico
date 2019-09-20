function ParaMat = simcolumn_ParaMat_ver2(SynPar, preinfo, postinfo, CM, DM)
% generate synapses parameter matrix
% Am is log-normal distributed
% CV and Fail are fit with bernoli distribution, with similar emperical
% means compared with recorded data


% get number of different cell types in L4 as well as in L2/3
pretype = unique(preinfo(:,4));
Nt_pre = [];
for i = 1:length(pretype)
    Nt_pre(pretype(i)) = length(find(preinfo(:,4) == pretype(i)));
end
posttype = unique(postinfo(:,4));
Nt_post = [];
for i = 1:length(posttype)
    Nt_post(posttype(i)) = length(find(postinfo(:,4) == posttype(i)));
end
% keyboard
%generate synaptic efficacy matrix
% amplitude, rise and decay time constant, short-term plasticity
[Am, Trise, Tfall, Plas, CV, Fail, Delay] = ...
    simcolumn_generateSynparas_ver2(SynPar, Nt_pre, Nt_post, CM, DM);

ParaMat.Am = Am;
ParaMat.Trise = Trise;
ParaMat.Tfall = Tfall;
ParaMat.Plas = Plas;
ParaMat.CV = CV;
ParaMat.Fail = Fail;
ParaMat.Delay = Delay;

% keyboard
% % synaptic delay
% Delay = simcolumn_generateDelayMat(CM, DM, Lpre, Lpost, ...
%     Nt_pre, Nt_post);
