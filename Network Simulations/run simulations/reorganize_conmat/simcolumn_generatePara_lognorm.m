function ParaMat = simcolumn_generatePara_lognorm(p, Msize)
% keyboard
[mu, sigma] = convert_normtolognorm(p);
ParaMat = para_lognorm(mu, sigma, Msize(1), Msize(2), p);