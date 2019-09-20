function ParaMat = simcolumn_assemble_ParaMat(cellinfo, Mats)
% assemble parameter matrix to be used to run simulation
% this function is specifically used to assemble All-to-All matrix (matrix
% not include input layer cells)
% cellinfo provide a reference to number of cells in each compartment of
% the ParaMat; Mats is compartment part of the overall ParaMat, which is
% expected to be a cell with length equals (length(cellinfo))^2
% the order of the Mats is expected to be the same as matlab linear
% indexing, i.e. row-to-column, with size of (length(cellinfo), length(cellinfo))

Mat_size = length(cellinfo);
% check if Mats is compatible with cellinfo
if length(Mats) ~= Mat_size^2
    error('the number of input Matrix is not compatible with input cell information')
end

Ncell = [];
for i = 1:length(cellinfo)
    Ncell(i) = size(cellinfo{i}, 1);
end
Nt = [0, cumsum(Ncell)];
% ParaMat is a sqaure matrix
ParaMat = zeros(sum(Ncell));

% assemble ParaMat
for i = 1:length(Mats)
    % matrix indexing of current Mats part
    [r, c] = ind2sub(Mat_size, i);
    % indexing inside the ParaMat
    if ~isempty(Mats{i})
        ParaMat(Nt(r)+1:Nt(r+1), Nt(c)+1:Nt(c+1)) = Mats{i};
    end
end

ParaMat = sparse(ParaMat);