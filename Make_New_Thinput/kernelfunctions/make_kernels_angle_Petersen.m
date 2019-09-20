function [Kernels, Kernels_SD, Kernels_grid] = make_kernels_angle_Petersen(Nkernel, binsize, params)
% Example use
% params_akernels.distance_grid = 3;              % (mm) from Petersen 2008
% params_akernels.STDstim = 70/1000;              % (mm) from Petersen 2008
% params_akernels.nkeep = round(100/binsize);     % how many data points to keep from Petersen kernels (sets kerneltime)
% params_akernels.noiseamp = 0;                   % activation functions
% [Kernels_ba, ~, ~] = make_kernels_angle_Petersen(Nkernel, binsize, params_akernels);

Kernels_SD.ActivationFunction.function = @activation_Petersen_2008; % sigmoid
Kernels.ActivationFunction.function = @activation_Petersen_2008; % sigmoid

%% Base angle kernels from Chao (from grid displacement kernels, Petersen 2008)
ks.Para = Petersen2008_FG(Nkernel, binsize);

% reorganize the kernels
kernels_temp.F = [ks.Para.MonoUnits.F; ks.Para.Bi_VolUnits.F; ...
    ks.Para.Bi_VolPosUnits.F; ks.Para.PolyUnits.F];
kernels_temp.G.v = [ks.Para.MonoUnits.G.v; ks.Para.Bi_VolUnits.G.v; ...
    ks.Para.Bi_VolPosUnits.G.v; ks.Para.PolyUnits.G.v];
kernels_temp.G.k = [ks.Para.MonoUnits.G.k; ks.Para.Bi_VolUnits.G.k; ...
    ks.Para.Bi_VolPosUnits.G.k; ks.Para.PolyUnits.G.k];
kernels_temp.G.b = [ks.Para.MonoUnits.G.b; ks.Para.Bi_VolUnits.G.b; ...
    ks.Para.Bi_VolPosUnits.G.b; ks.Para.PolyUnits.G.b];
kernels_temp.G.q = [ks.Para.MonoUnits.G.q; ks.Para.Bi_VolUnits.G.q; ...
    ks.Para.Bi_VolPosUnits.G.q; ks.Para.PolyUnits.G.q];

Kernels_SD.Kernels = cell(Nkernel,1);
Kernels_SD.ActivationFunction.Params = cell(Nkernel,1);
% y(t) = k./(1+q.*exp(-b.*t)).^v
% NB v, k and b are scaled so the right firing rates are found in a
% 'Petersen-like' experiment (grid displacement, white noise stim, STD = 1,
% Gauss filter sigma = 1.6 ms)
if isfield(params, 'vfac')
    vfac = params.vfac;
else
    vfac = 1;
end
if isfield(params, 'kfac')
    kfac = params.kfac;
else
    kfac = 1;
end
if isfield(params, 'bfac')
    bfac = params.bfac;
else
    bfac = 1;
end
if isfield(params, 'qfac')
    qfac = params.qfac;
else
    qfac = 1;
end
    
for nk = 1:Nkernel
    Kernels_SD.Kernels{nk} = fliplr(kernels_temp.F(nk,end-params.nkeep:end))';  % Chao used xcorr instead of conv
    Kernels_SD.ActivationFunction.Params{nk}.v = kernels_temp.G.v(nk,:)*vfac;   % independent of STDstim
    % NB Somewhere along the lines Chao increased k, so undo this (/4)
    Kernels_SD.ActivationFunction.Params{nk}.k = kernels_temp.G.k(nk,:)*kfac/4;      % independent of STDstim, independent of binsize, max firing rate
    Kernels_SD.ActivationFunction.Params{nk}.b = kernels_temp.G.b(nk,:)*bfac;        % scales with STDstim
    Kernels_SD.ActivationFunction.Params{nk}.q = kernels_temp.G.q(nk,:)*qfac;      % independent of STDstim
    Kernels_SD.ActivationFunction.Params{nk}.noiseamp = params.noiseamp;        % Chao: 0.1;
end

Kernels_SD.kerneltime = (binsize*(0:length(Kernels_SD.Kernels{1})-1))'; 

%% Transform to mm
Kernels_grid = Kernels_SD;
for nk = 1:Nkernel
    Kernels_grid.Kernels{nk} = Kernels_grid.Kernels{nk}*params.STDstim;
    Kernels_grid.ActivationFunction.Params{nk}.b = Kernels_grid.ActivationFunction.Params{nk}.b/(params.STDstim^2);
end


%% Transform to base angle
WK = Whisker_Recording;
WK.grid_displacement_recordings = Kernels_grid.Kernels;
WK.grid_displacement_absrel = 'rel';
WK.grid_displacement_unit = 'mm';
WK.distance_grid = params.distance_grid; % (mm)
WK = WK.transform('base_angle');

% Put in Kernelstruct
Kernels.kerneltime = Kernels_SD.kerneltime; 
for nk = 1:Nkernel
    Kernels.Kernels{nk,1} = WK.base_angle_recordings{nk}; % NB Chao used xcorr, not conv

    Kernels.ActivationFunction.Params{nk} = Kernels_SD.ActivationFunction.Params{nk};
    Kernels.ActivationFunction.Params{nk}.b = Kernels.ActivationFunction.Params{nk}.b/(params.STDstim/WK.distance_grid)^2;
    Kernels.ActivationFunction.Params{nk}.noiseamp = params.noiseamp; % Chao: 0.1;
    
end

%% Reorganize
new_position_vec = randperm(Nkernel);
Kernels.ActivationFunction.Params = Kernels.ActivationFunction.Params(new_position_vec);
Kernels.Kernels = Kernels.Kernels(new_position_vec);
Kernels_SD.ActivationFunction.Params = Kernels_SD.ActivationFunction.Params(new_position_vec);
Kernels_SD.Kernels = Kernels_SD.Kernels(new_position_vec);
Kernels_grid.ActivationFunction.Params = Kernels_grid.ActivationFunction.Params(new_position_vec);
Kernels_grid.Kernels = Kernels_grid.Kernels(new_position_vec);

end
%% Helper functions 
% Original from Chao Huang
function Para = Petersen2008_FG(Ntot, binsize, seed)

% fs and gs for the LNP model, from Petersen et al 2008 (Neuron)
% values are estimated from the paper, without accurate statistics;
% especially for the multi-feature group
% sigmas should be log-normally distributed; asymmatry index alpha should
% also be
% Ntot = 100; % total number of thalamic neurons in one barrel loids
if nargin == 2
    seed = [];
end
if isempty(seed)
    rng('shuffle');
else
    rng(seed);
end
seeds = randi(10000, [8,1]);

t = -100:binsize:0; % time vector for the kernels, with 1ms resolution

Para = [];
Para.Ntot = Ntot;
%% monophsic units
% sigma 2.0+-1.0, highly directional selective with
% alpha=0.93+-0.09; low firing rates; about 11% of the population
% first generate sigma
P = [];
Nc = round(Ntot*0.14);
P.NG = 1;
P.NC = Nc;
P.Sig = simcolumn_generatePara_lognorm([2.0, 1.0], [Nc, 1], seeds(1));
% peak value: seems around 0.4SD, which translate to 0.01 radius; assuming
% a CV of 0.15
P.Pk = 0.4 + 0.4*0.15*randn(size(P.Sig));
% prefered angle, assuming half protraction, half retraction
RI = randperm(Nc);
P.Pk(RI < (Nc/2)) = P.Pk(RI < Nc/2)*-1;
% peak location; assuming -2*sig + a random negtive number
P.Pl = -2*P.Sig - 3*1.8*rand(size(P.Sig));
% generate kernel f
P.F = GaussianKernelF_func(P, t);
% kernel g is modeled with sigmoid function
% consider logistic function form y(t) = k./(1+q.*exp(-b.*t)).^v
% in which k is the maximum, q related to y(0), b determines slope, and v
% affects slope around asymtote; 
% v should be around 2; b is good in range of 1.2-2; k is 10-40; q should
% in range of 2.2-3
P.G.v = 2+0.2*randn(size(P.Sig));
P.G.k = 5*(20+5*rand(size(P.Sig)));
P.G.b = 1.2+0.8*rand(size(P.Sig));
P.G.q = 2.4+0.8*rand(size(P.Sig));
% generate tuning function g (not needed, but see if anything goes wrong)
G = TuningFuncG(P, -4:0.05:4);

Para.MonoUnits = P;
%% biphasic units
%% 1st group: volecity units, 0.25 of total population
% with 2 roughly equal but opposite gaussian; relatively fast sigma
% ~1.6ms; high firing rates; highly directional
P = [];
Nc = round(Ntot*0.31);
P.NG = 2;
P.NC = Nc;
P.Sig(:,1) = simcolumn_generatePara_lognorm([1.8, 0.8], [Nc, 1], seeds(2));
P.Pk(:,1) = 0.32 + 0.32*0.15*randn(size(P.Sig));
% prefered angle, assuming half protraction, half retraction
RI = randperm(Nc);
P.Pk(RI < (Nc/2)) = P.Pk(RI < Nc/2)*-1;
% peak location; assuming -2*sig + a random negtive number
P.Pl(:,1) = -2*P.Sig - 3*1.8*rand(size(P.Sig));
% the second gaussian; relatively the same parameters as the first one
P.Sig(:,2) = P.Sig(:,1)+0.05*P.Sig(:,1).*rand(size(P.Sig(:,1)));
P.Pk(:,2) = -1*P.Pk(:,1)+0.05*P.Pk(:,1).*rand(size(P.Sig(:,1)));
% location need to be checked
P.Pl(:,2) = P.Pl(:,1)-2*P.Sig(:,2) - 3*1.2*rand(size(P.Sig(:,2)));
% generate kernel f
P.F = GaussianKernelF_func(P, t);
% tuning function can be roughly the same as in monophasic units; however
% the firing rate is at least twice larger
P.G.v = 2+0.2*randn(size(P.Sig(:,1)));
P.G.k = 4.5*(65+15*randn(size(P.Sig(:,1))));
P.G.b = 1.2+0.8*rand(size(P.Sig(:,1)));
P.G.q = 2.4+0.8*rand(size(P.Sig(:,1)));
% generate tuning function g (not needed, but see if anything goes wrong)
G = TuningFuncG(P, -4:0.05:4);
Para.Bi_VolUnits = P;

%% 2nd group: Pos/Vol units, 0.33 of population; two gaussian with different
% size, smaller peaks
P = [];
Nc = round(Ntot*0.41);
P.NG = 2;
P.NC = Nc;
P.Sig(:,1) = simcolumn_generatePara_lognorm([1.8, 0.8], [Nc, 1], seeds(3));
P.Pk(:,1) = 0.3 + 0.3*0.3*randn(size(P.Sig));
% prefered angle, assuming half protraction, half retraction
RI = randperm(Nc);
P.Pk(RI < (Nc/2)) = P.Pk(RI < Nc/2)*-1;
% peak location; assuming -2*sig + a random negtive number
P.Pl(:,1) = -2*P.Sig(:,1) - 3*1.7*rand(size(P.Sig(:,1)));
% second peak
P.Sig(:,2) = simcolumn_generatePara_lognorm([2, 0.9], [Nc, 1], seeds(4));
% Pk value is still opposite to the 1st peak; however the value should
% differ 
Pk = P.Pk;
idx = find(abs(P.Pk(:,1)) < 0.3);
Pk(idx) = -1*Pk(idx).*(1.3+0.7*rand(size(idx))); 
idx = find(abs(P.Pk(:,1)) >= 0.3);
Pk(idx) = -1*Pk(idx).*(0.3+0.5*rand(size(idx))); 
P.Pk(:,2) = Pk;
% peak location; assuming -2*sig + a random negtive number
P.Pl(:,2) = P.Pl(:,1) - 2*P.Sig(:,2) - 3*1.3*rand(size(P.Sig(:,2)));
% generate kernel f
P.F = GaussianKernelF_func(P, t);
% tuning function has slightly lower slope, with possible bi-directional
% tuning
P.G.v = 2+0.2*randn(size(P.Sig(:,1)));
P.G.k = 4.5*(75+20*randn(size(P.Sig(:,1))));
P.G.b = 1+0.8*rand(size(P.Sig(:,1)));
P.G.q = 2.4+0.8*rand(size(P.Sig(:,1)));
% assuming 20% of the units have significant bidirectional tuning
RI = randperm(Nc);
RI = RI(1:round(0.2*Nc));
P.G.v(:, 2) = 3+0.2*randn(size(P.Sig(:,1)));
P.G.k(RI, 2) = 4.5*(40+20*rand(size(P.Sig(RI,1))));
P.G.b(:, 2) = -0.4+0.2*rand(size(P.Sig(:,1)));
P.G.q(:, 2) = 4+0.8*rand(size(P.Sig(:,1)));
G = TuningFuncG(P, -4:0.05:4);
Para.Bi_VolPosUnits = P;

%% polyphasic units
% about 11% of the population, with 3-5 gaussian kernels
P = [];
% Nc = round(Ntot*0.11);
Nc = Para.Ntot - Para.MonoUnits.NC - Para.Bi_VolUnits.NC - Para.Bi_VolPosUnits.NC;
P.NG = 2;
P.NC = Nc;
P.Sig(:,1) = simcolumn_generatePara_lognorm([1.8, 0.8], [Nc, 1], seeds(5));
P.Pk(:,1) = 0.25 + 0.25*0.3*randn(size(P.Sig));
% prefered angle, assuming half protraction, half retraction
RI = randperm(Nc);
P.Pk(RI < (Nc/2)) = P.Pk(RI < Nc/2)*-1;
% peak location; assuming -2*sig + a random negtive number
P.Pl(:,1) = -2*P.Sig(:,1) - 3*1.5*rand(size(P.Sig(:,1)));
% second peak
P.Sig(:,2) = simcolumn_generatePara_lognorm([2, 0.9], [Nc, 1], seeds(6));
% Pk value is still opposite to the 1st peak; however the value should
% differ 
Pk = P.Pk;
idx = find(abs(P.Pk(:,1)) < 0.25);
Pk(idx) = -1*Pk(idx).*(1.3+0.7*rand(size(idx))); 
idx = find(abs(P.Pk(:,1)) >= 0.25);
Pk(idx) = -1*Pk(idx).*(0.3+0.5*rand(size(idx))); 
P.Pk(:,2) = Pk;
% peak location; assuming -2*sig + a random negtive number
P.Pl(:,2) = P.Pl(:,1) - 2*P.Sig(:,2) - 3*1.3*rand(size(P.Sig(:,2)));
% third peak
P.Sig(:,3) = simcolumn_generatePara_lognorm([2, 1], [Nc, 1], seeds(7));
% Pk value is still opposite to the 1st peak; however the value should
% differ 
P.Pk(:,3) = 0.07 + 0.07*0.5*randn(size(P.Sig(:,1)));
% peak location; assuming -2*sig + a random negtive number
P.Pl(:,3) = P.Pl(:,2) - 2*P.Sig(:,3) - 3*1.1*rand(size(P.Sig(:,3)));
% fourth peak
P.Sig(:,4) = simcolumn_generatePara_lognorm([2, 1.1], [Nc, 1], seeds(8));
% Pk value is still opposite to the 1st peak; however the value should
% differ 
P.Pk(:,4) = 0.04 + 0.04*0.5*randn(size(P.Sig(:,1)));
RI = randperm(Nc);
P.Pk(RI(1:round(0.5*RI)), 4) = 0;
% peak location; assuming -2*sig + a random negtive number
P.Pl(:,4) = P.Pl(:,3) - 2*P.Sig(:,4) - 3*1.3*rand(size(P.Sig(:,4)));
% generate kernel f
P.F = GaussianKernelF_func(P, t);
% tuning function has slightly lower slope, with possible bi-directional
% tuning
P.G.v = 2+0.2*randn(size(P.Sig(:,1)));
P.G.k = 4.5*(70+15*randn(size(P.Sig(:,1))));
P.G.b = 1.0+0.6*rand(size(P.Sig(:,1)));
P.G.q = 2.4+0.8*rand(size(P.Sig(:,1)));
% assuming 30% of the units have significant bidirectional tuning
RI = randperm(Nc);
RI = RI(1:round(0.3*Nc));
P.G.v(:, 2) = 3+0.2*randn(size(P.Sig(:,1)));
P.G.k(RI, 2) = 4.5*(40+10*rand(size(P.Sig(RI,1))));
P.G.b(:, 2) = -0.4+0.2*rand(size(P.Sig(:,1)));
P.G.q(:, 2) = 4+0.8*rand(size(P.Sig(:,1)));
G = TuningFuncG(P, -4:0.05:4);
Para.PolyUnits = P;

%% Re-evaluate Ntot
Para.Ntot = Para.MonoUnits.NC+Para.Bi_VolUnits.NC+Para.Bi_VolPosUnits.NC+Para.PolyUnits.NC;
%% check the size of kernel G
Para.MonoUnits.G.v(:,2) = Para.MonoUnits.G.v(:,1); 
Para.MonoUnits.G.k(:,2) = 0; 
Para.MonoUnits.G.b(:,2) = Para.MonoUnits.G.b(:,1); 
Para.MonoUnits.G.q(:,2) = Para.MonoUnits.G.q(:,1); 

Para.Bi_VolUnits.G.v(:,2) = Para.Bi_VolUnits.G.v(:,1); 
Para.Bi_VolUnits.G.k(:,2) = 0; 
Para.Bi_VolUnits.G.b(:,2) = Para.Bi_VolUnits.G.b(:,1); 
Para.Bi_VolUnits.G.q(:,2) = Para.Bi_VolUnits.G.q(:,1); 


end


function ParaMat = simcolumn_generatePara_lognorm(p, Msize, seed)
    if nargin == 2
        seed = [];
    end
    [mu, sigma] = convert_normtolognorm(p);
    
    ParaMat = para_lognorm(mu, sigma, Msize(1), Msize(2), p, seed);
end

function [mu, sigma] = convert_normtolognorm(p)
    %convert parameter for normal distribution to parameter for log-normal
    %distribution. the parameters generated can be used in lognrnd function to
    %generate log-normal distributions with mean and std specified in input p

    m = p(1);
    v = p(2)^2;
    mu = log((m^2)/sqrt(v+m^2));
    sigma = sqrt(log(v/(m^2)+1));
end

function M = para_lognorm(mu, sigma, Npost, Npre, p, seed)
    %generate log-normal parameter based on mu, sigma information 
    if nargin == 5
        seed = [];
    end
    if isempty(seed)
        rng('shuffle');
    else
        rng(seed);
    end
    
    M = lognrnd(mu,sigma,Npost,Npre);
    cutoff = max(p(1)+3*p(2), 5);
    % cutoff(cutoff>13) = 13;
    M(M>cutoff) = cutoff;
end

function k = GaussianKernelF_func(P, t)
    % get numeric kernel from kernel parameters P
    % 

    k = [];

    for i = 1:length(t)
        k(:,i) = sum(P.Pk.*exp(-(t(i)-P.Pl).^2./P.Sig.^2), 2);
    end
    
end

function k = TuningFuncG(P, t)
    % get numeric kernel from kernel parameters P
    % 

    k = [];

    for i = 1:length(t)
        k(:,i) = sum(P.G.k./(1+P.G.q.*exp(-P.G.b.*t(i))).^P.G.v, 2);
    end
end

