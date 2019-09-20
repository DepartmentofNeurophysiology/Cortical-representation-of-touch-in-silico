function Kernels = make_kernels_curvature(Nkernel, binsize, params)
% Example use
% params_ckernels.mf = 10;
% params_ckernels.spreadf = 40;
% params_ckernels.mv = 2; % normal
% params_ckernels.spreadv = 0.2;
% params_ckernels.mk = 80; % normal
% params_ckernels.spreadk = 40;
% params_ckernels.mb = 40; % uniform
% params_ckernels.spreadb = 1;
% params_ckernels.mq = 0; % uniform
% params_ckernels.spreadq = 2;
% params_ckernels.noiseamp = 0;                   % activation functions    
% params_ckernels.kerneltime = (0:100)*binsize;
% Kernels = make_kernels_curvature(Nkernel, binsize, params_ckernels);


Kernels.ActivationFunction.function = @activation_Petersen_2008;

for nk = 1:Nkernel;
                
    % Make sin kernel
    fkernel = (params.mf+params.spreadf*rand)/1000; 
    taukernel = (2*rand)*(1/(1*fkernel));
    if rand<0.5
        signkernel = -1;
    else
        signkernel = 1;
    end
    kernel = signkernel*(heaviside(params.kerneltime).*(1-heaviside(params.kerneltime-1/fkernel))).*exp(-params.kerneltime/taukernel).*sin(2*pi*fkernel*params.kerneltime);
    kernel = kernel./sqrt(sum(kernel.^2)*binsize);
    Kernels.Kernels{nk,1} = kernel;

    Kernels.ActivationFunction.Params{nk}.v = [params.mv+params.spreadv*randn,0];
    Kernels.ActivationFunction.Params{nk}.k = [params.mk+params.spreadk*randn,0];
    Kernels.ActivationFunction.Params{nk}.b = [params.mb+params.spreadb*rand,0];
    Kernels.ActivationFunction.Params{nk}.q = [params.mq+params.spreadq*rand,0];
    Kernels.ActivationFunction.Params{nk}.noiseamp = params.noiseamp;

                
end

Kernels.kerneltime = params.kerneltime;

%% Reorganize
new_position_vec = randperm(Nkernel);
Kernels.ActivationFunction.Params = Kernels.ActivationFunction.Params(new_position_vec);
Kernels.Kernels = Kernels.Kernels(new_position_vec);

end