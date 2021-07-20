function Denddist = simcolumn_connectivity_DendFunc(Ref_point, DendPara, X, Y, Z)
% axon distribution function, a modified gaussian function
% AxonPara is determined by function simcolumn_connectivity_assignDenAxon
% Ref_point is the zero point on which the distribution is calculated

% Model_Space = make_Model_Space();

% stepsize = Model_Space.stepsize;

% [X, Y, Z] = meshgrid(-800:stepsize:800, -800:stepsize:800, -1000:stepsize:1000);
DendCenter = DendPara(1:3);
DendSpan = DendPara(4:6);

% transform the coordinate so the reference point is at (0, 0, 0)
DendCenter = DendCenter - Ref_point;

% % distribution in X, Y, Z
% Denddist = (exp(-(X - DendCenter(1)).^2/DendSpan(1)^2) - 0.2*exp(-(X - DendCenter(1)).^2/(DendSpan(1)+20)^2))...
%     .*(exp(-(Y - DendCenter(2)).^2/DendSpan(2)^2) - 0.2*exp(-(Y - DendCenter(2)).^2/(DendSpan(2)+20)^2))...
%     .*(exp(-(Z - DendCenter(3)).^2/DendSpan(3)^2) - 0.2*exp(-(Z - DendCenter(3)).^2/(DendSpan(3)+20)^2));

% use matlab build-in mvnpdf function
mu = DendCenter;
sigma = [DendSpan(1)^2 0 0; 0 DendSpan(2)^2 0; 0 0 DendSpan(3)^2];
sigma1 = [(DendSpan(1)-5)^2 0 0; 0 (DendSpan(2)-5)^2 0; 0 0 (DendSpan(3)-5)^2];

Denddist = mvnpdf([X(:) Y(:) Z(:)], mu, sigma) - 0.1*mvnpdf([X(:) Y(:) Z(:)], mu, sigma1);

Denddist = reshape(Denddist, size(X));