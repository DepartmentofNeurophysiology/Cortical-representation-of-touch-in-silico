function Axondist = simcolumn_connectivity_AxonFunc(Ref_point, AxonPara, X, Y, Z)
% axon distribution function, a modified gaussian function
% AxonPara is determined by function simcolumn_connectivity_assignDenAxon
% Ref_point is the zero point on which the distribution is calculated

% global Model_Space;

% stepsize = Model_Space.stepsize;
% keyboard
% [X, Y, Z] = meshgrid(-800:stepsize:800, -800:stepsize:800, -1000:stepsize:1000);
AxonCenter = AxonPara(1:3);
AxonSpan = AxonPara(4:6);

% transform the coordinate so the reference point is at (0, 0, 0)
AxonCenter = AxonCenter - Ref_point;
% keyboard
% % distribution in X, Y, Z
% Axondist = (exp(-(X - AxonCenter(1)).^2/AxonSpan(1)^2) - 0.5*exp(-(X - AxonCenter(1)).^2/(AxonSpan(1)+20)^2))...
%     .*(exp(-(Y - AxonCenter(2)).^2/AxonSpan(2)^2) - 0.5*exp(-(Y - AxonCenter(2)).^2/(AxonSpan(2)+20)^2))...
%     .*(exp(-(Z - AxonCenter(3)).^2/AxonSpan(3)^2) - 0.5*exp(-(Z - AxonCenter(3)).^2/(AxonSpan(3)+20)^2));
% 
% % need to write the function with a normalization factor so each
% % distribution is constant 
% Axondist = exp(-(X - AxonCenter(1)).^2/AxonSpan(1)^2);

% use matlab build-in mvnpdf function
mu = AxonCenter;
sigma = [AxonSpan(1)^2 0 0; 0 AxonSpan(2)^2 0; 0 0 AxonSpan(3)^2];
sigma1 = [(AxonSpan(1)-0.2*AxonSpan(1))^2 0 0; 0 (AxonSpan(2)-0.2*AxonSpan(2))^2 0; 0 0 (AxonSpan(3)-0.2*AxonSpan(3))^2];

Axondist = mvnpdf([X(:) Y(:) Z(:)], mu, sigma) - 0.1*mvnpdf([X(:) Y(:) Z(:)], mu, sigma1);
% Axondist(Axondist > 0.9 * max(Axondist)) = 0.9*max(Axondist);
% keyboard
Axondist = reshape(Axondist, size(X));
