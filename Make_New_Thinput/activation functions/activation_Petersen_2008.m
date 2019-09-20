function PSTH = activation_Petersen_2008(ConvTrace, Params)
% Input: 
% * ConvTrace: ConvTrace = 1xNtime array with a grid displacement recording convolved with a grid displacement kernel
% * Params: Parameters of the activation function for the
% kernels for this dimension
% Output:
% PSTH = 1xNtime array with # spikes / sec (so independent of binsize!)


[Ndimc, Ntime] = size(ConvTrace);
[~, Ndimp ] = size(Params);

if Ndimc>1 
    error('Number of dimensions of the convolved recordings should be 1')
end

if Ndimp>1
    error('Number of dimensions of the parameters should be 1')
end

PSTH = zeros(1,Ntime);
for ndp = 1:2
    PSTH = PSTH+Params.k(ndp)./((1+Params.q(ndp).*exp(-Params.b(ndp).*ConvTrace)).^Params.v(ndp));
end
PSTH = (PSTH + Params.noiseamp.*PSTH.*randn(size(PSTH))); % unit: # spikes / sec

end