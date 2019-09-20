function PSTH = activation_sigmoid_linear_ND(ConvTrace, Params)
% Input: 
% * ConvTrace: ConvTrace = NdimxNtime array with a grid displacement recording convolved with a grid displacement kernel
% * Params{dimension}: Parameters of the activation function for the
% kernels for this dimension
% * binsize (ms)
% * plotyn (0 or 1): (optional) whether to plot or not


[Ndimc, Ntime] = size(ConvTrace);
[~, Ndimp ] = size(Params);

if ~(Ndimc == Ndimp)
    error('Number of dimensions of the convolved recordings and relevant parameters for the activation function should be the same')
end


ConvTrace_combined = zeros(size(ConvTrace(1,:)));
for nd = 1:Ndimp
    % Combine the dimensions into a single PSTH
    ConvTrace_combined = ConvTrace_combined + Params{nd}.dimfactor*ConvTrace(nd,:);
end
PSTH = Params{1}.k./((1+Params{1}.q.*exp(-Params{1}.b.*ConvTrace_combined)).^Params{nd}.v);
end