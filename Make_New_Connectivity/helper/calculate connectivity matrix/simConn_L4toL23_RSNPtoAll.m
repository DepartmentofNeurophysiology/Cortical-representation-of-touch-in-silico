function CM = simConn_L4toL23_RSNPtoAll(ADO_Mat, postIND, postType, AxonRatio, DendRatio, DMs)
% ADO_Mat: axon-dendritic overlap matrix
% postIND: indexing to seperate different type of post synaptic cells
% postType: numbers specify the type of post-synaptic cell
% no connection is known at this stage, set all to zero

CM = zeros(size(ADO_Mat));
return

% keyboard
% connectivity matrix need to be calculated for each type of post synaptic
% cell
for i = 1:length(postIND) - 1
    temp = ADO_Mat(postIND(i)+1:postIND(i+1), :);

    ADO_temp = temp;
    DM = DMs(postIND(i)+1:postIND(i+1), :);
    DR = DendRatio(postIND(i)+1:postIND(i+1));

    switch postType(i)
        case(1) % post-synaptic cell is pyramidal
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_RSNPtoPyr(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{2, 3} % post-synaptic cell is RSNP or chandiler
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_RSNPtoFsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(4) % post-synaptic cell is RsPV
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_RSNPtoRsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{5, 6} % post-synaptic cell is MarSOM
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_RSNPtoMarSOM(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{7, 8} % post-synaptic cell is BipVIP
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_RSNPtoBipVIP(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(9) % post-synaptic cell is BipCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_RSNPtoBipCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(10) % post-synaptic cell in SbcCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_RSNPtoSbcCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(11) % post-synaptic cell in NG cells
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_RSNPtoNG(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
    end
end

%% nested function
    function CMs = simConn_L4toL23_RSNPtoPyr(ADO_temp, AR, DR)
        % map axon-dendrite overlapping index into connectivity matrix
        % the value is normalized to meausre average connection probability
        % from experiment; for Pyr-Pyr connection, within 50 micron p =
        % 0.19
        % also total number of connection is controlled by the axon and
        % dendrite ratio in the model region
        % assuming a total convengenc/divergence rate of 140 for total
        % connection
        
        pConn = ADO_temp/1.6454e-12*0.2;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 300*AR', 180*DR');
    end


    function CMs = simConn_L4toL23_RSNPtoFsPV(ADO_temp, AR, DR)
        % for Pyr-FsPV connection, within 50 micron p =
        % 0.85
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/1.2171e-12*0.45;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 50*AR', 210*DR');
    end


    function CMs = simConn_L4toL23_RSNPtoRsPV(ADO_temp, AR, DR)
        % for Pyr-FsPV connection, within 150 micron p =
        % 0.20
        % assuming a total convengenc rate of 240 and divergence rate of 40
        
        pConn = ADO_temp/1.2930e-12*0.4;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 45*AR', 210*DR');
    end


    function CMs = simConn_L4toL23_RSNPtoMarSOM(ADO_temp, AR, DR)
        % for Pyr-MarSOM connection, no connection
        % assuming a total convengenc rate of 700 and divergence rate of 90
%         
%         pConn = ADO_temp/3.1268e-12*0.7;
%         % convert pConn into binary connectivity matrix
%         RI = rand(size(pConn));
        CMs = zeros(size(ADO_temp));
%         CMs(RI < pConn) = 1;
%         
%         CMs = conn_reduction(CMs, 70*AR', 700*DR');
    end

    function CMs = simConn_L4toL23_RSNPtoBipVIP(ADO_temp, AR, DR)
        % for Pyr-RSNP connection, within 150 micron p =
        % 0.21
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/1.2968e-12*0.4;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 40*AR', 210*DR');
    end

    function CMs = simConn_L4toL23_RSNPtoBipCR(ADO_temp, AR, DR)
        % for Pyr-RSNP connection, within 150 micron p =
        % 0.183
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/1.2486e-12*0.4;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 35*AR', 180*DR');
    end

    function CMs = simConn_L4toL23_RSNPtoSbcCR(ADO_temp, AR, DR)
        % for Pyr-RSNP connection, within 150 micron p =
        % 0.174
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/1.2328e-12*0.4;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 30*AR', 180*DR');
    end

    function CMs = simConn_L4toL23_RSNPtoNG(ADO_temp, AR, DR)
        % for Pyr-RSNP connection, within 150 micron p =
        % 0.6
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/1.2429e-12*0.3;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 30*AR', 180*DR');
    end
end

% %%
% %check connection probability as a function of distance
% bin = 0:25:500;
% N_all = [0, cumsum(NiN2)];
% Pconn = {};
% for post = 1:length(N_all) - 1
%     Pconn{post,1} = PConn_pairs(CM_modified(N_all(post)+1:N_all(post+1),:), ...
%         bin, DM(N_all(post)+1:N_all(post+1), :));
% end