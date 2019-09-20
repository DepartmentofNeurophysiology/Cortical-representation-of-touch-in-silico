function [F_list, preSpikes, timestepseed] = simcolumn_synapse_onestep(ALL_In, Am_Mat, Plas_Mat, CV_Mat, Pf_Mat,...
    Delay_Mat, preID_mat, Npost, F_list, tm, timestepseed)
% generate a list of synaptic input parameters, which contains:
% 1st column: linear index of the synapse
% 2nd column: synaptic weight, used to update s in the synaptic model
% 3rd column: synaptic delay
% 4th column: linear index of the synapse event in the homogeneous time
% constant matrix
% tm is the current time on the simulation
% generate one instance of synaptic weight based on Am and CV

% Plas_Mat is running in GPU
if exist('timestepseed','var') 
    rng(timestepseed)
else
    rng('shuffle')
    scurr = rng;
    timestepseed = scurr.Seed;
end

idx = (ALL_In);
Am = abs(Am_Mat(idx).*(1+CV_Mat(idx).*randn(size(idx))));

%clipping the response
RI = find(abs(Am)>3*abs(Am_Mat(idx)));
Am(RI) = abs(3*Am_Mat(idx(RI)));
% now add sign to the Am 
Am = Am.*sign(Am_Mat(idx));
% generate failure
% first apply STD on failure rate
STD_value = gather(Plas_Mat(idx));
Pf = Pf_Mat(idx)./STD_value;
test = rand(size(idx));
fail_index = find((Pf - test) > 0);
% apply STD
Am = Am.*STD_value;
% remove failed inputs
idx(fail_index) = [];
Am(fail_index) = [];
% get synaptic delay
% Delay = Delay_Mat(idx) + 0.15*randn(size(idx));
Delay = Delay_Mat(idx);
% Delay(Delay < 0.3) = 0.3;
Delay = Delay + tm;

% remove NaN value from Am
Am(isnan(Am)) = 0;

% calcualte the linear index for the new homogeneous time constant matrix
r = preID_mat(idx, 1);
c = preID_mat(idx, 2);
%update F list
inter = [idx, Am, Delay, (c-1)*Npost + r];
F_list = [F_list; inter];

%% preSpikes is used to run STDP; contains spike timing for each spikes,
% including synaptic delay
preSpikes = [idx, Delay];





