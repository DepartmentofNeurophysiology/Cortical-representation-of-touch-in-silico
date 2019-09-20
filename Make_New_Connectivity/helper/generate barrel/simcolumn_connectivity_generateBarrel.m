function [l23info_1, l23para_1, l4info_1, l4para_1, l23type, l4type] = ...
    simcolumn_connectivity_generateBarrel(l23info, l4info, xy_shift, barrelID)

% get a copy of original cellinfo
temp = l23info(:, 1:3);
% shafle l23info
idx = randperm(size(temp, 1));
temp = temp(idx, :);
% modify x, y location
temp(:,1) = temp(:,1) - xy_shift(1);
temp(:,2) = temp(:,2) - xy_shift(2);
% assign cell type 
[l23para_1,l23info_1,l23type] = celltype_assign(temp,0.84,'l23');
% assign barrel identifier
l23info_1(:, 5) = barrelID;
% assign E/I identifier
l23info_1(l23info_1(:,4) == 1,6) = 1;
l23info_1(l23info_1(:,4) > 1, 6) = -1;
% l4 info
temp = l4info(:, 1:3);
% shafle l23info
idx = randperm(size(temp, 1));
temp = temp(idx, :);
% modify x, y location
temp(:,1) = temp(:,1) - xy_shift(1);
temp(:,2) = temp(:,2) - xy_shift(2);
% assign cell type 
[l4para_1,l4info_1,l4type] = celltype_assign(temp,0.88,'l4');
% barrel identifier
l4info_1(:,5) = barrelID;
% assign E/I identifier
l4info_1(l4info_1(:,4) < 3, 6) = 1;
l4info_1(l4info_1(:,4) > 2, 6) = -1;