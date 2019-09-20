function Plas_mat = simcolumn_setPlaspara(Paras, Npostpre)
% set short-term dynamics matrix for the simulation program
% need to clip the matrix based on the mean value of parameter

Plas_mat = simcolumn_generatePara_lognorm(Paras, Npostpre);
if Paras(1)< 0.9
    Plas_mat(Plas_mat>0.95) = 0.95;
    Plas_mat(Plas_mat<0.35) = 0.35;
elseif Paras(1) > 1.25
    Plas_mat(Plas_mat<1.05) = 1.05;
else
    Plas_mat(Plas_mat<0.8) = 0.8;
    Plas_mat(Plas_mat>1.2) = 1.2;
end