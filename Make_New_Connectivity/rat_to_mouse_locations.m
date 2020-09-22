function cellinfo_mouse = rat_to_mouse_locations(cellinfo_rat)
%% Funcion to get back to relevant mouse cell locations. 
% The measured cell densities are in mice, but most axon/dendritic distribution patterns, connectivity and synaptic efficacy in the literature are from rats. 
% To make these fit, the measured cell densities were scaled to fit rat data. So the resulting cell locations are appropriate for rat. 
% Once the connectivity data are generated, these cell locations are not used in the network simulations (as neurons are simulated as point neurons). 

disp('convert cell locations to mouse locations')
xyz_rat = cellinfo_rat(:,1:3);
xyz_mouse(:,1) = 2*xyz_rat(:,1)/3;       % x location (~rostral-caudal)
xyz_mouse(:,2) = 2*xyz_rat(:,2)/3;       % y location (~dorsal-ventral)
xyz_mouse(:,3) = xyz_rat(:,3)/2;         % z location (depth, distance from pia)

cellinfo_mouse = cellinfo_rat;
cellinfo_mouse(:,1:3) = xyz_mouse;
