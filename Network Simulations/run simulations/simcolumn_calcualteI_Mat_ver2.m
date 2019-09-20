function i = simcolumn_calcualteI_Mat_ver2(g_mat, EItype, vm)
% calculate synaptic current from running variable g_mat
% keyboard
% first remove NaN values
% g_mat(isnan(g_mat)) = 0;

% check  = g_mat(:, EItype == -1);

% seperate excitatory and inhibitory inputs
iE = sum(g_mat(:, EItype == 1), 2);
iI = sum(g_mat(:, EItype == -1), 2);

% modulation of current by post-synaptic membrane potential
iE = iE.*(vm/-70);
% iI = iI.*(vm+80)./25;

% total current
i = iE + iI;