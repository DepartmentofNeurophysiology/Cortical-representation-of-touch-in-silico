function M = para_lognorm(mu, sigma, Npost, Npre, p)
%generate log-normal parameter based on mu, sigma information 

M = lognrnd(mu,sigma,Npost,Npre);
cutoff = max(p(1)+3*p(2), 5);
% cutoff(cutoff>13) = 13;
M(M>cutoff) = cutoff;