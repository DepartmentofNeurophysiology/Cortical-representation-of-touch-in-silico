function PSTH = activation_table_ND(ConvTrace, Table)
% Input: 
% * ConvTrace: ConvTrace = NdimxNtime array with a grid displacement recording convolved with a grid displacement kernel
% * Table: Lookup-table with for each input value an output 
%       * Table.input =(1xNdim) cell array, with each a (1xNvaluethisdim) array with input values for this dimension
%       * Table.output =(Nvaluex x Nvaluey x ...) array, with output values for each combination of input values  
% Output:
% PSTH = 1xNtime array with # spikes / sec (so independent of binsize!)

% Example (3D sigmoid):
% Table.input{1} = -10:0.1:10;
% Table.input{2} = -5:0.1:5;
% Table.input{3} = -12:0.1:12;
% Table.output = nan*ones(length(Table.input{1}),length(Table.input{2}), length(Table.input{3}));
% for nx = 1:length(Table.input{1})
%     for ny = 1:length(Table.input{2})
%         for nz = 1:length(Table.input{3})
%             Table.output(nx,ny,nz) = 1./(1+exp(-(Table.input{1}(nx)+2.*Table.input{2}(ny) - 1.5*Table.input{3}(nz))));
%         end
%     end
% end
% CT = randn(3,1000);
% PSTH = activation_table_ND(CT, Table);



[Ndimc, Ntime] = size(ConvTrace);
[X, Ndimi] = size(Table.input);
Ndimo = ndims(Table.output);
if (Ndimi == 1) && (X >1)
    Table.input = Table.input';
    [~, Ndimi] = size(Table.input);
end
if ~(Ndimc == Ndimi)
    error('Number of input dimensions should match the number of dimensions of the recordings')
end
if ~(Ndimi == Ndimo)
    error('Number of input dimensions should match the number of output dimensions')
end

evalstr = 'Table.output(';
for nd = 1:Ndimo
    if nd == Ndimo
        evalstr = [evalstr, 'mindi(', num2str(nd), '));'];
    else
        evalstr = [evalstr, 'mindi(', num2str(nd), '),'];
    end
end

PSTH = zeros(1,Ntime);
for nt = 1:Ntime
    % find input values
    mindi = nan*ones(1,Ndimo);
    for nd = 1:Ndimo
        [~, mindi(nd)] = min(abs(ConvTrace(nd,nt)-Table.input{nd}));
    end
    % calculate output
    try
        PSTH(nt) = eval(evalstr);
    catch
        keyboard
    end
end

end