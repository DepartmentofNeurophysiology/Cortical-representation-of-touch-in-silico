function [Am, Trise, Tfall, Plas, CV, Fail, Delay] = ...
    simcolumn_generateSynparas_ver2(SynPar, Ntypes_pre, Ntypes_post, CMl23tol23, DM)
% keyboard

% add zero in to Ntypes
Nt_pre = [0, cumsum(Ntypes_pre)];
Nt_post = [0, cumsum(Ntypes_post)];

Am=[];
Trise = [];
Tfall = [];
Plas = [];
CV = [];
Fail = [];
Delay = [];

for pre = 1:length(Nt_pre) - 1
    for post = 1:length(Nt_post) - 1
        
        Am_temp = simcolumn_generatePara_lognorm(SynPar{post, pre}.I0, [Ntypes_post(post), Ntypes_pre(pre)]); % amplitude is log-normal distributed
        % set maximum of Am
        Am_temp(Am_temp > min([5*SynPar{post, pre}.I0(1), 8])) = min([5*SynPar{post, pre}.I0(1), 8]);
        Am(Nt_post(post) + 1: Nt_post(post + 1), Nt_pre(pre) + 1: Nt_pre(pre + 1)) = Am_temp;
        
        Trise_temp = abs(SynPar{post, pre}.tau2(1) + SynPar{post, pre}.tau2(2)*randn(Ntypes_post(post), Ntypes_pre(pre))); % normally distributed
        Tfall_temp = abs(SynPar{post, pre}.tau1(1) + SynPar{post, pre}.tau1(2)*randn(Ntypes_post(post), Ntypes_pre(pre))); % normally distributed
        %make sure tau_fall is larger than tau_rise
        Tfall_temp(Tfall_temp<Trise_temp) = Trise_temp(Tfall_temp<Trise_temp) + (SynPar{post, pre}.tau1(1)-SynPar{post, pre}.tau2(1));
        Plas(Nt_post(post) + 1: Nt_post(post + 1), Nt_pre(pre) + 1: Nt_pre(pre + 1)) = ...
            simcolumn_setPlaspara(SynPar{post, pre}.Plas, [Ntypes_post(post), Ntypes_pre(pre)]); % normally distributed
        
        Trise(Nt_post(post) + 1: Nt_post(post + 1), Nt_pre(pre) + 1: Nt_pre(pre + 1)) = Trise_temp;
        Tfall(Nt_post(post) + 1: Nt_post(post + 1), Nt_pre(pre) + 1: Nt_pre(pre + 1)) = Tfall_temp;
% keyboard        

        %% CV and failure rate
        % to keep it consistent with earlier version, add the sign of
        % synaptic weight matrix, i.e. excitatory with +1 and inhibitory
        % with -1
        Am(Nt_post(post) + 1: Nt_post(post + 1), Nt_pre(pre) + 1: Nt_pre(pre + 1)) = ...
            Am(Nt_post(post) + 1: Nt_post(post + 1), Nt_pre(pre) + 1: Nt_pre(pre + 1)) .* ...
            CMl23tol23(Nt_post(post) + 1: Nt_post(post + 1), Nt_pre(pre) + 1: Nt_pre(pre + 1));
        
        [~, ~, Fail(Nt_post(post) + 1: Nt_post(post + 1), Nt_pre(pre) + 1: Nt_pre(pre + 1)), ...
            CV(Nt_post(post) + 1: Nt_post(post + 1), Nt_pre(pre) + 1: Nt_pre(pre + 1))] = ...
            generate_SynPara(Am(Nt_post(post) + 1: Nt_post(post + 1), Nt_pre(pre) + 1: Nt_pre(pre + 1)), SynPar{post, pre});
        
        %% delay matrix
        Delay(Nt_post(post) + 1: Nt_post(post + 1), Nt_pre(pre) + 1: Nt_pre(pre + 1)) = ...
            DM(Nt_post(post) + 1: Nt_post(post + 1), Nt_pre(pre) + 1: Nt_pre(pre + 1))/SynPar{post, pre}.delay(1) ... 
            + SynPar{post, pre}.delay(2) + 0.15*SynPar{post, pre}.delay(2)*randn(Ntypes_post(post), Ntypes_pre(pre));
    end
end
% keyboard
%no parameters should be lower than 0
Am(Am<0)=-1*Am(Am<0);
Am(isnan(Am)) = 0;

%convert the maxtrix to sparse form
Am=sparse(Am.*CMl23tol23);
Trise=sparse(abs(Trise.*CMl23tol23));
Tfall=sparse(abs(Tfall.*CMl23tol23));
Plas(Plas>2) = 2;
Plas=sparse(abs(Plas.*CMl23tol23));
Delay(Delay < 0.4) = 0.401;
Delay=sparse(abs(Delay.*CMl23tol23));
CV = sparse(CV.*abs(CMl23tol23));
Fail = sparse(Fail.*abs(CMl23tol23));
