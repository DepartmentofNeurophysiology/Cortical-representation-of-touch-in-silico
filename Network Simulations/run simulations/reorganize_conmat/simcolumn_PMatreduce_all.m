function PMmat_new = simcolumn_PMatreduce_all(PMat_old, CM_new)
% further reduce the demionsion of parameter matrix need to run the simulation
% the previous version have full matrix; now only use sparse matrix, in
% which non-connected pairs are dropped; the new matrix will be
% Nconn-by-1 size, with each non-zero entry as values of synapses
% keyboard

idx = abs(CM_new);
PMmat_new = PMat_old(idx);

if length(PMmat_new) ~= length(CM_new)
    error('dimension mismatch between parameter matrix and connectivity matrix')
end