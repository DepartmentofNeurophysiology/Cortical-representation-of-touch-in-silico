function Model_Space = make_Model_Space()

% the global parameter model_space
% defines the following values:
% 1: starting and ending points in x-y domain
% 2: border of each layer

% global Model_Space;

Model_Space.row = [-450, 450]; % this is a x coordinate of the cellinfo matrix
Model_Space.arc = [-450, 450];
Model_Space.L2 = [0, 430];
Model_Space.L4 = [430, 630];
Model_Space.stepsize = 20; % step size (resolution) used to calculate axon or dendrite density; in micron

% information for each barrel
Model_Space.barrels.row = {[-450, -150], [-150, 150], [150, 450], [150, 450], [-150, 150], [-450, -150], [-450, -150], [-150, 150], [150, 450]};
Model_Space.barrels.arc = {[-150, 150], [-150, 150], [-150, 150], [150, 450], [150, 450], [150, 450], [-450, -150], [-450, -150], [-450, -150]};
Model_Space.barrels.lable = [1, 2, 3, 4, 5, 6, 7, 8, 9];