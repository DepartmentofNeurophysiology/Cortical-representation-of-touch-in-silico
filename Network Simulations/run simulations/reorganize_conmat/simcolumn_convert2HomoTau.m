function Tau_mat_homo = simcolumn_convert2HomoTau(Tau_mat, Nt_pre, Nt_post)
% generate homogeneous time constrant metrix from connection-specifiy time
% constant matrix
% Tau_mat: connection-specifiy time constant matrix
% Nt_pre : an array contains number of neurons in each pre-synaptic neuron
% group

N_pre= [0, cumsum(Nt_pre)];
N_post = [0, cumsum(Nt_post)];

Tau_mat_homo = zeros(size(Tau_mat, 1), length(Nt_pre));

for j = 1:length(N_post) - 1
    for i = 1:length(N_pre) - 1
        %     keyboard
        inter = Tau_mat(N_post(j) + 1 : N_post(j+1), N_pre(i) + 1 : N_pre(i+1));
        Tau_mat_homo(N_post(j) + 1 : N_post(j+1), i) = ...
            mean(inter(find(inter)));
    end
end

Tau_mat_homo(isnan(Tau_mat_homo)) = 0;
