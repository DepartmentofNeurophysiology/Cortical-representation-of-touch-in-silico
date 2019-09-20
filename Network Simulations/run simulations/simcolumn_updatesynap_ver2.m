function [s_mat, plas_mat, F_list] = simcolumn_updatesynap_ver2(s_mat, plas_mat, plas_ValueMat, F_list, tm)
% update synaptic current and STD terms at each time step using synaptic 
% input list generated from function simcolumn_Matform_synlist

% find the synapse inputs that should reach post-synaptic side at this 
% time point 
idx = find(F_list(:,3) <= tm);

if ~isempty(idx)
%     count = length(idx) - length(unique(F_list(idx, 4)));
    
    % add those event into the s matrix
%     s_mat(F_list(idx, 4)) = s_mat(F_list(idx, 4)) + F_list(idx, 2);
    for i = 1:length(idx)
        s_mat(F_list(idx(i), 4)) = s_mat(F_list(idx(i), 4)) + F_list(idx(i), 2);
    end
    % update STD running variable
    temp = gather(plas_mat(F_list(idx, 1))).*plas_ValueMat(F_list(idx, 1));
    temp(temp < 0.4) = 0.4;
    temp(temp > 5) = 5;
    plas_mat(F_list(idx, 1)) = temp;
%     try
%         plas_mat(F_list(idx, 1)) = temp;
%     catch
%         [~, I, ~] = unique(F_list(idx, 1));
%         plas_mat((F_list(idx(I), 1))) = temp(I);
%         disp(length(idx) - length(I))
%     end
    
%     plas_mat(F_list(idx, 1)) = plas_mat(F_list(idx, 1)).*plas_ValueMat(F_list(idx, 1));
    % need to remove these events from F_list
    F_list(idx, :) = [];
    
%     % set limits on plas mat
%     plas_mat(plas_mat < 0.4) = 0.4;
%     plas_mat(plas_mat > 5) = 5;
end