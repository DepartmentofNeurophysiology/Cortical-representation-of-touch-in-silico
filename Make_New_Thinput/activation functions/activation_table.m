function PSTH = activation_table(ConvTrace, Table)
% Input: 
% * ConvTrace: ConvTrace = 1xNtime array with a grid displacement recording convolved with a grid displacement kernel
% * Table: Lookup-table with for each input value an output (Table.
% Output:
% PSTH = 1xNtime array with # spikes / sec (so independent of binsize!)


[Ndimc, Ntime] = size(ConvTrace);
[Nvaluei, Ndimpi ] = size(Table.input);
[Nvalueo, Ndimpo ] = size(Table.output);

if Ndimc>1
    error('Number of dimensions of the convolved recordings should be 1')
end

if Ndimpi>1 || Ndimpo>1
    error('Number of dimensions of the lookup table should be 1')
end

if ~Nvaluei == Nvalueo
    error('Number of inputvalues should equal number of output values')
end

PSTH = zeros(1,Ntime);
for nt = 1:Ntime
    [~, mindi] = min(abs(ConvTrace(nt)-Table.input));
    PSTH(nt) = Table.output(mindi);
end

end