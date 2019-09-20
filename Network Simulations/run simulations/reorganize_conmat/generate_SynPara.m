function [Qs, Pr, Pf, cv] = generate_SynPara(Am, Para)
%generate Qs, Pr, Pf and cv from Am, Nb, variation
% keyboard
if isfield(Para,'Nb') == 0
    Nb = 3;
else
    Nb = Para.Nb;
end
if isfield(Para,'va') == 0
    va = 0.2*mean(Am(find(Am)));
else
    va = Para.va;
end

ind = find(Am);
if isempty(ind) == 0
    Qs_temp = Am(ind)/Nb; %quantum content assuming Pr = 1
    Qs_temp = Qs_temp + va*rand(length(ind),1);
    Pr_temp = Am(ind)./Qs_temp/Nb;
    Pr_temp(Pr_temp<0.05) = 0.05;
    Pf_temp = (1-Pr_temp).^Nb;
    cv_temp = ((1-Pr_temp)./Pr_temp./Nb).^0.5;
    Pf_temp(Pf_temp>0.8) = 0.8;
    
    
    Qs = zeros(size(Am));
    Qs(ind) = Qs_temp;
    Pr = zeros(size(Am));
    Pr(ind) = Pr_temp;
    Pf = zeros(size(Am));
    Pf(ind) = Pf_temp;
    cv_temp(cv_temp>1.5) = 1.5;
    if isfield(Para, 'cv') == 1 %modulate cv
        cv_temp = cv_temp./mean(cv_temp).*Para.cv;
    end
    cv_temp(cv_temp>1.5) = 1.5;
    cv = zeros(size(Am));
    cv(ind) = cv_temp;
    
    if isreal(cv_temp) == 0
        keyboard
    end
else
    Qs = zeros(size(Am));
    Pr = zeros(size(Am));
    Pf = zeros(size(Am));
    cv = zeros(size(Am));
end

% keyboard

return

mean(Pf(find(Pf)))
std(Pf(find(Pf)))

mean(cv(find(cv)))
std(cv(find(cv)))