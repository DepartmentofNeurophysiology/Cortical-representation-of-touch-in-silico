function [mu, sigma] = convert_normtolognorm(p)
%convert parameter for normal distribution to parameter for log-normal
%distribution. the parameters generated can be used in lognrnd function to
%generate log-normal distributions with mean and std specified in input p

m = p(1);
v = p(2)^2;
mu = log((m^2)/sqrt(v+m^2));
sigma = sqrt(log(v/(m^2)+1));