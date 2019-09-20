function CM = simConn_L23toL23_SbcCRtoAll(ADO_Mat, postIND, postType, AxonRatio, DendRatio, DMs)
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
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_SbcCRtoPyr(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{2, 3} % post-synaptic cell is FsPV or chandiler
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_SbcCRtoFsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(4) % post-synaptic cell is RsPV
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_SbcCRtoRsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{5, 6} % post-synaptic cell is MarSOM
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_SbcCRtoMarSOM(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{7, 8} % post-synaptic cell is BipVIP
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_SbcCRtoBipVIP(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(9) % post-synaptic cell is BipCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_SbcCRtoBipCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(10) % post-synaptic cell in SbcCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_SbcCRtoSbcCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(11) % post-synaptic cell in NG cells
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_SbcCRtoNG(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
    end
end

%% nested function
    function CMs = simConn_L23toL23_SbcCRtoPyr(ADO_temp, AR, DR)
        % map axon-dendrite overlapping index into connectivity matrix
        % the value is normalized to meausre average connection probability
        % from experiment; for SbcCR-Pyr connection, within 150 micron p =
        % 0.142
        % also total number of connection is controlled by the axon and
        % dendrite ratio in the model region
        % assuming a total convengenc/divergence rate of 140 for total
        % connection
        
        pConn = ADO_temp/2.8154e-12*0.142;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 210*AR', 15*DR');
    end


    function CMs = simConn_L23toL23_SbcCRtoFsPV(ADO_temp, AR, DR)
        % for SbcCR-FsPV connection, within 150 micron p =
        % 0.2
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.6071e-12*0.2;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 18*AR', 15*DR');
    end


    function CMs = simConn_L23toL23_SbcCRtoRsPV(ADO_temp, AR, DR)
        % for SbcCR-RsPV connection, 0.5
        % assuming a total convengenc rate of 240 and divergence rate of 40
        
        pConn = ADO_temp/2.6930e-12*0.5;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 21*AR', 18*DR');
    end


    function CMs = simConn_L23toL23_SbcCRtoMarSOM(ADO_temp, AR, DR)
        % for SbcCR-MarSOM connection, 0.5 with 150
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.7468e-12*0.5;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 24*AR', 24*DR');
    end

    function CMs = simConn_L23toL23_SbcCRtoBipVIP(ADO_temp, AR, DR)
        % for SbcCR-BipVIP connection, within 150 micron p =
        % 0.5
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.8968e-12*0.5;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 21*AR', 21*DR');
    end

    function CMs = simConn_L23toL23_SbcCRtoBipCR(ADO_temp, AR, DR)
        % for SbcCR-BipCR connection, 0.5
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.8668e-12*0.5;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 21*AR', 15*DR');
    end

    function CMs = simConn_L23toL23_SbcCRtoSbcCR(ADO_temp, AR, DR)
        % for SbcCR-SbcCR connection, within 150 micron p =
        % 0.1
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.5968e-12*0.15;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 13*AR', 12*DR');
    end

    function CMs = simConn_L23toL23_SbcCRtoNG(ADO_temp, AR, DR)
        % for SbcCR-NG connection, within 150 micron p =
        % 0.45
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.4968e-12*0.45;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 21*AR', 15*DR');
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