function [CM_new, preCell, preidx, postidx] = simcolumn_CMreduce_all(CM_old, cellinfo_pre)
% further reduce the demionsion of matrix need to run the simulation
% the previous version have full matrix; now only use sparse matrix, in
% which non-connected pairs are dropped; the new matrix will be
% Nconn-by-1 size, with each non-zero entry as linear index of connected
% pairs neurons in the full connectivity matrix
% also a new, preCell mat is generated, which makes assigning spikes easier

% create two connecvitivy lists, in which store idx of all post (pre)
% synaptic parteners for each pre (post) synaptic neurons
% the indexs are in reference to reduced N-by-1 matrix, in which each entry
% is a synapse

idx = find(CM_old);

CM_new = idx.*CM_old(idx);

[idxpost, idxpre] = find(CM_old);

preCell = [idxpost, idxpre, cellinfo_pre(idxpre, 6)];

postidx = {};
for i = 1:size(CM_old, 2)
    postidx{i} = find(idxpre == i);
end

preidx = {};
for i = 1:size(CM_old, 1)
    preidx{i} = find(idxpost == i);
end

% keyboard