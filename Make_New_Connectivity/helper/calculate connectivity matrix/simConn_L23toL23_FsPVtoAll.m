function CM = simConn_L23toL23_FsPVtoAll(ADO_Mat, postIND, postType, AxonRatio, DendRatio, DMs)
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
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_FsPVtoPyr(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{2, 3} % post-synaptic cell is FsPV or chandiler
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_FsPVtoFsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(4) % post-synaptic cell is RsPV
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_FsPVtoRsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{5, 6} % post-synaptic cell is MarSOM
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_FsPVtoMarSOM(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{7, 8} % post-synaptic cell is BipVIP
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_FsPVtoBipVIP(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(9) % post-synaptic cell is BipCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_FsPVtoBipCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(10) % post-synaptic cell in SbcCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_FsPVtoSbcCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(11) % post-synaptic cell in NG cells
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_FsPVtoNG(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
    end
end

%% nested function
    function CMs = simConn_L23toL23_FsPVtoPyr(ADO_temp, AR, DR)
        % map axon-dendrite overlapping index into connectivity matrix
        % the value is normalized to meausre average connection probability
        % from experiment; for FsPV-Pyr connection, within 50 micron p =
        % 0.8
        % also total number of connection is controlled by the axon and
        % dendrite ratio in the model region
        % assuming a total convengenc/divergence rate of 140 for total
        % connection
        
        ADO_temp(ADO_temp < 4e-10) = 0.01* ADO_temp(ADO_temp < 4e-10);
        pConn = ADO_temp/2.4488e-8*0.9;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 920*AR', 70*DR');
    end


    function CMs = simConn_L23toL23_FsPVtoFsPV(ADO_temp, AR, DR)
        % for FsPV-FsPV connection, within 50 micron p =
        % 0.8
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.4151e-8*0.8;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 40*AR', 40*DR');
    end


    function CMs = simConn_L23toL23_FsPVtoRsPV(ADO_temp, AR, DR)
        % for FsPV-FsPV connection, no connection reported
        % assuming a total convengenc rate of 240 and divergence rate of 40
        
%         pConn = ADO_temp/2.6930e-12*0.20;
%         % convert pConn into binary connectivity matrix
%         RI = rand(size(pConn));
        CMs = zeros(size(ADO_temp));
%         CMs(RI < pConn) = 1;
%         
%         CMs = conn_reduction(CMs, 240*AR', 40*DR');
    end


    function CMs = simConn_L23toL23_FsPVtoMarSOM(ADO_temp, AR, DR)
        % for FsPV-FsPV connection, no connection reported
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
%         pConn = ADO_temp/2.9968e-12*0.7;
%         % convert pConn into binary connectivity matrix
%         RI = rand(size(pConn));
        CMs = zeros(size(ADO_temp));
%         CMs(RI < pConn) = 1;
%         
%         CMs = conn_reduction(CMs, 700*AR', 90*DR');
    end

    function CMs = simConn_L23toL23_FsPVtoBipVIP(ADO_temp, AR, DR)
        % for FsPV-BipVIP connection, within 150 micron p =
        % 0.297
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.5364e-8*0.297;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 15*AR', 25*DR');
    end

    function CMs = simConn_L23toL23_FsPVtoBipCR(ADO_temp, AR, DR)
        % for FsPV-BipCR connection, within 150 micron p =
        % 0.297
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.6177e-8*0.297;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 15*AR', 25*DR');
    end

    function CMs = simConn_L23toL23_FsPVtoSbcCR(ADO_temp, AR, DR)
        % for FsPV-SbcCR connection, within 150 micron p =
        % 0.2
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.5851e-8*0.2;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 15*AR', 20*DR');
    end

    function CMs = simConn_L23toL23_FsPVtoNG(ADO_temp, AR, DR)
        % for FsPV-NG connection, within 150 micron p =
        % 0.6
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/2.7929e-8*0.6;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 20*AR', 30*DR');
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