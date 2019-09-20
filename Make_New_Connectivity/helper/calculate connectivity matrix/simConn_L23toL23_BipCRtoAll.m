function CM = simConn_L23toL23_BipCRtoAll(ADO_Mat, postIND, postType, AxonRatio, DendRatio, DMs)
% ADO_Mat: axon-dendritic overlap matrix
% postIND: indexing to seperate different type of post synaptic cells
% postType: numbers specify the type of post-synaptic cell
% keyboard
% connectivity matrix need to be calculated for each type of post synaptic
% cell
for i = 1:length(postIND) - 1
    temp = ADO_Mat(postIND(i)+1:postIND(i+1), :);
    
%     ADO_temp = temp;
%     DM = DMs(postIND(i)+1:postIND(i+1), :);
%     DR = DendRatio(postIND(i)+1:postIND(i+1));
    
    switch postType(i)
        case(1) % post-synaptic cell is pyramidal
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_BipCRtoPyr(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{2, 3} % post-synaptic cell is FsPV or chandiler
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_BipCRtoFsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(4) % post-synaptic cell is RsPV
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_BipCRtoRsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{5, 6} % post-synaptic cell is MarSOM
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_BipCRtoMarSOM(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{7, 8} % post-synaptic cell is BipVIP
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_BipCRtoBipVIP(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(9) % post-synaptic cell is BipCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_BipCRtoBipCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(10) % post-synaptic cell in SbcCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_BipCRtoSbcCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(11) % post-synaptic cell in NG cells
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_BipCRtoNG(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
    end
end

%% nested function
    function CMs = simConn_L23toL23_BipCRtoPyr(ADO_temp, AR, DR)
        % map axon-dendrite overlapping index into connectivity matrix
        % the value is normalized to meausre average connection probability
        % from experiment; for BipCR-Pyr connection, within 150 micron p =
        % 0.112
        % also total number of connection is controlled by the axon and
        % dendrite ratio in the model region
        % assuming a total convengenc/divergence rate of 140 for total
        % connection
        
        pConn = ADO_temp/2.8954e-12*0.112;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 160*AR', 15*DR');
    end


    function CMs = simConn_L23toL23_BipCRtoFsPV(ADO_temp, AR, DR)
        % for BipCR-FsPV connection, within 150 micron p =
        % 0.297
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/3.1171e-12*0.297;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 18*AR', 15*DR');
    end


    function CMs = simConn_L23toL23_BipCRtoRsPV(ADO_temp, AR, DR)
        % for BipCR-RsPV connection, 0.297
        % assuming a total convengenc rate of 240 and divergence rate of 40
        
        pConn = ADO_temp/2.6-12*0.297;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 18*AR', 15*DR');
    end


    function CMs = simConn_L23toL23_BipCRtoMarSOM(ADO_temp, AR, DR)
        % for BipCR-MarSOM connection, 0.8 with 150
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/3.2468e-12*0.8;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 34*AR', 21*DR');
    end

    function CMs = simConn_L23toL23_BipCRtoBipVIP(ADO_temp, AR, DR)
        % for BipCR-BipVIP connection, within 150 micron p =
        % 0.38
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/3.2968e-12*0.38;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 21*AR', 15*DR');
    end

    function CMs = simConn_L23toL23_BipCRtoBipCR(ADO_temp, AR, DR)
        % for BipCR-BipCR connection, 0.38
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/3.0668e-12*0.38;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 18*AR', 15*DR');
    end

    function CMs = simConn_L23toL23_BipCRtoSbcCR(ADO_temp, AR, DR)
        % for BipCR-SbcCR connection, within 150 micron p =
        % 0.38
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/3.0968e-12*0.38;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 18*AR', 15*DR');
    end

    function CMs = simConn_L23toL23_BipCRtoNG(ADO_temp, AR, DR)
        % for BipCR-NG connection, within 150 micron p =
        % 0.38
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.6968e-12*0.38;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 21*AR', 18*DR');
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