function CM = simConn_L23toL23_RsPVtoAll(ADO_Mat, postIND, postType, AxonRatio, DendRatio)
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
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_RsPVtoPyr(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{2, 3} % post-synaptic cell is FsPV or chandiler
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_RsPVtoFsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(4) % post-synaptic cell is RsPV
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_RsPVtoRsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{5, 6} % post-synaptic cell is MarSOM
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_RsPVtoMarSOM(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{7, 8} % post-synaptic cell is BipVIP
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_RsPVtoBipVIP(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(9) % post-synaptic cell is BipCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_RsPVtoBipCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(10) % post-synaptic cell in SbcCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_RsPVtoSbcCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(11) % post-synaptic cell in NG cells
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_RsPVtoNG(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
    end
end

%% nested function
    function CMs = simConn_L23toL23_RsPVtoPyr(ADO_temp, AR, DR)
        % map axon-dendrite overlapping index into connectivity matrix
        % the value is normalized to meausre average connection probability
        % from experiment; for RsPV-Pyr connection, within 100 micron p =
        % 0.41
        % also total number of connection is controlled by the axon and
        % dendrite ratio in the model region
        % assuming a total convengenc/divergence rate of 140 for total
        % connection
        
        pConn = ADO_temp/3.2596e-12*0.41;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 480*AR', 24*DR');
    end


    function CMs = simConn_L23toL23_RsPVtoFsPV(ADO_temp, AR, DR)
        % for RsPV-FsPV connection, within 50 micron p =
        % 0.26
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.9e-12*0.26;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 18*AR', 15*DR');
    end


    function CMs = simConn_L23toL23_RsPVtoRsPV(ADO_temp, AR, DR)
        % for FsPV-FsPV connection, within 100 micron 0.41
        % assuming a total convengenc rate of 240 and divergence rate of 40
        
        pConn = ADO_temp/2.7258e-12*0.41;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 18*AR', 18*DR');
    end


    function CMs = simConn_L23toL23_RsPVtoMarSOM(ADO_temp, AR, DR)
        % for RsPV-MarSOM connection, 0.41 with 150
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/3.0604e-12*0.41;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 21*AR', 18*DR');
    end

    function CMs = simConn_L23toL23_RsPVtoBipVIP(ADO_temp, AR, DR)
        % for RsPV-BipVIP connection, within 150 micron p =
        % 0.41
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/3.1129e-12*0.41;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 17*AR', 17*DR');
    end

    function CMs = simConn_L23toL23_RsPVtoBipCR(ADO_temp, AR, DR)
        % for RsPV-BipCR connection, no connection
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
%         pConn = ADO_temp/2.9968e-12*0.7;
%         % convert pConn into binary connectivity matrix
%         RI = rand(size(pConn));
        CMs = zeros(size(ADO_temp));
%         CMs(RI < pConn) = 1;
%         
%         CMs = conn_reduction(CMs, 700*AR', 90*DR');
    end

    function CMs = simConn_L23toL23_RsPVtoSbcCR(ADO_temp, AR, DR)
        % for RsPV-SbcCR connection, within 150 micron p =
        % 0.33
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.9968e-12*0.33;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 17*AR', 15*DR');
    end

    function CMs = simConn_L23toL23_RsPVtoNG(ADO_temp, AR, DR)
        % for RsPV-NG connection, within 150 micron p =
        % 0.33
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.6646e-12*0.33;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 13*AR', 13*DR');
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