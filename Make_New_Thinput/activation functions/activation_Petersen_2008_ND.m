function PSTH = activation_Petersen_2008_ND(ConvTrace, Params)
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

PSTH = zeros(1,Ntime);

PSTHperdim = cell(1,Ndimp);
for nd = 1:Ndimp
    % Calculate a PSTH per dimension
    PSTHperdim{nd} = zeros(1,Ntime);
    for ndp = 1:2
        PSTHperdim{nd} = PSTHperdim{nd}+Params{nd}.k(ndp)./((1+Params{nd}.q(ndp).*exp(-Params{nd}.b(ndp).*ConvTrace(nd,:))).^Params{nd}.v(ndp));
    end
    PSTHperdim{nd} = (PSTHperdim{nd} + Params{nd}.noiseamp.*(PSTHperdim{nd}).*randn(size(PSTHperdim{nd})));
    
    % Combine the dimensions into a single PSTH
    PSTH = PSTH + Params{nd}.dimfactor*PSTHperdim{nd};
end

end