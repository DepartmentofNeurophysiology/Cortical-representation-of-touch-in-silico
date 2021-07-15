function Neuron = simcolumn_connectivity_maxConn(Neuron)
% set maximum connection number of each neuron, calculated as percentage of
% axon (defines output number) and dendrite (defines input number) in the
% model region
% calculate percentage of axon and dendrite in the model space

% global Model_Space;
Model_Space = make_Model_Space();

step_size = Model_Space.stepsize;
row = Model_Space.row;
arc = Model_Space.arc;

% X, Y, Z is used to calculate total distribution
[X, Y, Z] = meshgrid((row(1) - 300):step_size:(row(2) + 300), ...
    (arc(1) - 300):step_size:(arc(2) + 300), -700:step_size:700);

% reference point is set to the center of the model space
Ref_point = [mean(row), mean(arc), mean(Model_Space.L2)];

% Xm, Ym, Zm is used to calculate distribution in the model region
[Xm, Ym, Zm] = meshgrid(row(1):step_size:row(2), arc(1):step_size:arc(2), -700:step_size:700);
% keyboard
% different axon distribution
for i = 1:length(Neuron.Axon)
    if ~strcmpi(Neuron.AxonRestriction{i}, 'home') 
        dist_all = simcolumn_connectivity_AxonFunc(Ref_point, Neuron.Axon{i}, X, Y, Z);
        dist_mod = simcolumn_connectivity_AxonFunc(Ref_point, Neuron.Axon{i}, Xm, Ym, Zm);
        Neuron.AxonConnRatio(i) = sum(dist_mod(:))/sum(dist_all(:));
    else
        Neuron.AxonConnRatio(i) = 1;
    end
%     keyboard
end
Neuron.AxonConnRatio(Neuron.AxonConnRatio > 1) = 1;

% different dendrite distribution
for i = 1:length(Neuron.Dendrite)
    if ~strcmpi(Neuron.DendriteRestriction{i}, 'home') 
        dist_all = simcolumn_connectivity_DendFunc(Ref_point, Neuron.Dendrite{i}, X, Y, Z);
        dist_mod = simcolumn_connectivity_DendFunc(Ref_point, Neuron.Dendrite{i}, Xm, Ym, Zm);
        Neuron.DendConnRatio(i) = sum(dist_mod(:))/sum(dist_all(:));
    else
        Neuron.DendConnRatio(i) = 1;
    end
end
Neuron.DendConnRatio(Neuron.DendConnRatio > 1) = 1;