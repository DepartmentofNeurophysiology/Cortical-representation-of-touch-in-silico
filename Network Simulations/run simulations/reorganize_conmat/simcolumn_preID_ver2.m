function preID = simcolumn_preID_ver2(CM, cellinfo_pre, cellinfo_post)
% generate a array with the same size as CM, in which each entry is the
% type of presynaptic neurons for each connection
% note: this function only works with N-by-1 connectivity matrix

Npre = size(cellinfo_pre, 1);
% Ntpre = [];
% Nt = unique(cellinfo_pre(:,4));
% for i = 1:length(Nt)
%     Ntpre(Nt(i)) = length(find(cellinfo_pre(:,4) == Nt(i)));
% end

Npost = size(cellinfo_post, 1);

% CM contains linear index of each connection in full matrix (size of
% Npost-by-Npre
[idxpost, idxpre] = ind2sub([Npost, Npre], abs(CM));
preID = [idxpost, cellinfo_pre(idxpre, 4), cellinfo_pre(idxpre, 6)];