function Kernels = make_kernels_angle_curvature_Svobodaexp(Nkernel_ba, Nkernel_c, Nkernel_m, binsize, savename)


if nargin == 4
    savename = [];
end


Nkernel = Nkernel_ba+Nkernel_c+Nkernel_m;
Ndim = 2;                                                           % dim 1 = base angle, dim 2 = curvature
Kernels.ActivationFunction.function = @activation_Petersen_2008_ND; % sigmoid

% NB % r(t) = k./(1+q.*exp(-b.*t)).^v
%% Base angle kernels
params_akernels.distance_grid = 3;              % (mm) from Petersen 2008
params_akernels.STDstim = 70/1000;              % (mm) from Petersen 2008
params_akernels.nkeep = round(100/binsize);     % how many data points to keep from Petersen kernels (sets kerneltime)
params_akernels.noiseamp = 0;                   % activation functions
params_akernels.vfac = 1;
params_akernels.kfac = 1;
params_akernels.bfac = 1/15;
params_akernels.qfac = 2;

%% Curvature kernels and activation functions
params_ckernels.mf = 10;
params_ckernels.spreadf = 40;

params_ckernels.mv = 3; % normal
params_ckernels.spreadv = 0.2;
params_ckernels.mk = 120; % normal
params_ckernels.spreadk = params_ckernels.mk/5;
params_ckernels.mb = 700; % uniform
params_ckernels.spreadb = params_ckernels.mb/5;
params_ckernels.mq = 4; % uniform
params_ckernels.spreadq = .5;

params_ckernels.noiseamp = 0;                   % activation functions
%% Base angle kernels from Chao (from grid displacement kernels, Petersen 2008)
% 'Petersen' kernels: assume they work on base angles
if Nkernel_ba>1
    [Kernels_ba, ~, ~] = make_kernels_angle_Petersen(Nkernel_ba, binsize, params_akernels);

    % Put in Kernelstruct
    Kernels.kerneltime = Kernels_ba.kerneltime; 
    for nk = 1:Nkernel_ba
        Kernels.ActivationFunction.Params{nk}=cell(1,Ndim);
        for ndim = 1:Ndim
            if ndim == 1
                Kernels.Kernels{nk, ndim} = Kernels_ba.Kernels{nk}; 
                Kernels.ActivationFunction.Params{nk}{ndim} = Kernels_ba.ActivationFunction.Params{nk};
                Kernels.ActivationFunction.Params{nk}{ndim}.dimfactor = 1;
                

            else
                Kernels.Kernels{nk, ndim} = zeros(size(Kernels.kerneltime));
                Kernels.ActivationFunction.Params{nk}{ndim}.v = [0,0];
                Kernels.ActivationFunction.Params{nk}{ndim}.k = [0,0];
                Kernels.ActivationFunction.Params{nk}{ndim}.b = [0,0];
                Kernels.ActivationFunction.Params{nk}{ndim}.q = [0,0];
                Kernels.ActivationFunction.Params{nk}{ndim}.noiseamp = 0;
                Kernels.ActivationFunction.Params{nk}{ndim}.dimfactor = 0;
            end
        end
    end
else
    Kernels.kerneltime = (0:100)*binsize;
end
%% Curvature Kernels
params_ckernels.kerneltime = Kernels.kerneltime;
if Nkernel_c>1
    Kernels_c = make_kernels_curvature(Nkernel_c, binsize, params_ckernels);
    for nk = Nkernel_ba+1:Nkernel_ba+Nkernel_c;
        Kernels.ActivationFunction.Params{nk}=cell(1,Ndim);
        for ndim = 1:Ndim
            if ndim == 1
                Kernels.Kernels{nk, ndim} = zeros(size(Kernels.kerneltime));
                Kernels.ActivationFunction.Params{nk}{ndim}.v = [0,0];
                Kernels.ActivationFunction.Params{nk}{ndim}.k = [0,0];
                Kernels.ActivationFunction.Params{nk}{ndim}.b = [0,0];
                Kernels.ActivationFunction.Params{nk}{ndim}.q = [0,0];
                Kernels.ActivationFunction.Params{nk}{ndim}.noiseamp = 0;
                Kernels.ActivationFunction.Params{nk}{ndim}.dimfactor = 0;
            else
                Kernels.Kernels{nk, ndim} = Kernels_c.Kernels{nk-Nkernel_ba};
                Kernels.ActivationFunction.Params{nk}{ndim} = Kernels_c.ActivationFunction.Params{nk-Nkernel_ba};
                Kernels.ActivationFunction.Params{nk}{ndim}.dimfactor = 1;                
            end
        end
    end
end
%% Mixed Kernels
if Nkernel_m>1
    Kernels_c2 = make_kernels_curvature(Nkernel_m, binsize, params_ckernels);
    [Kernels_ba2, ~, ~] = make_kernels_angle_Petersen(Nkernel_m, binsize, params_akernels);
    for nk = Nkernel_ba+Nkernel_c+1:Nkernel;
        Kernels.ActivationFunction.Params{nk}=cell(1,Ndim);
        for ndim = 1:Ndim
            if ndim == 1
                % base angle
                Kernels.Kernels{nk, ndim} = Kernels_ba2.Kernels{nk-Nkernel_ba-Nkernel_c}; 
                Kernels.ActivationFunction.Params{nk}{ndim} = Kernels_ba.ActivationFunction.Params{nk-Nkernel_ba-Nkernel_c};
                Kernels.ActivationFunction.Params{nk}{ndim}.dimfactor = rand;
            else
                % curvature
                Kernels.Kernels{nk, ndim} = Kernels_c2.Kernels{nk-Nkernel_ba-Nkernel_c};
                Kernels.ActivationFunction.Params{nk}{ndim} = Kernels_c.ActivationFunction.Params{nk-Nkernel_ba-Nkernel_c};
                Kernels.ActivationFunction.Params{nk}{ndim}.dimfactor = 1-Kernels.ActivationFunction.Params{nk}{1}.dimfactor;
            end 
            

        end
        
    end
end

%% Reorganize
new_position_vec = randperm(Nkernel);
Kernels.ActivationFunction.Params = Kernels.ActivationFunction.Params(new_position_vec);
Kernels.Kernels(:,:) = Kernels.Kernels(new_position_vec, :);
%% Save
if ~isempty(savename)
    save(savename,'Kernels')
end
