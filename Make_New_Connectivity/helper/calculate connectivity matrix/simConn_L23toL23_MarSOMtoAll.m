function CM = simConn_L23toL23_MarSOMtoAll(ADO_Mat, postIND, postType, AxonRatio, DendRatio)
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
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_MarSOMtoPyr(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{2, 3} % post-synaptic cell is FsPV or chandiler
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_MarSOMtoFsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(4) % post-synaptic cell is RsPV
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_MarSOMtoRsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{5, 6} % post-synaptic cell is MarSOM
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_MarSOMtoMarSOM(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{7, 8} % post-synaptic cell is BipVIP
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_MarSOMtoBipVIP(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(9) % post-synaptic cell is BipCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_MarSOMtoBipCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(10) % post-synaptic cell in SbcCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_MarSOMtoSbcCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(11) % post-synaptic cell in NG cells
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_MarSOMtoNG(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
    end
end

%% nested function
    function CMs = simConn_L23toL23_MarSOMtoPyr(ADO_temp, AR, DR)
        % map axon-dendrite overlapping index into connectivity matrix
        % the value is normalized to meausre average connection probability
        % from experiment; for MarSOM-Pyr connection, within 150 micron p =
        % 0.8
        % also total number of connection is controlled by the axon and
        % dendrite ratio in the model region
        % assuming a total convengenc/divergence rate of 140 for total
        % connection
        
        pConn = ADO_temp/3.9466e-12*0.8;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 760*AR', 30*DR');
    end


    function CMs = simConn_L23toL23_MarSOMtoFsPV(ADO_temp, AR, DR)
        % for MarSOM-FsPV connection, within 150 micron p =
        % 0.8
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/3.3171e-12*0.8;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 40*AR', 27*DR');
    end


    function CMs = simConn_L23toL23_MarSOMtoRsPV(ADO_temp, AR, DR)
        % for MarSOM-MarSOM connection, dense connection
        % assuming a total convengenc rate of 240 and divergence rate of 40
        
        pConn = ADO_temp/3.1930e-12*0.8;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(ADO_temp));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 40*AR', 27*DR');
    end


    function CMs = simConn_L23toL23_MarSOMtoMarSOM(ADO_temp, AR, DR)
        % for MarSOM-MarSOM connection, no connection
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
%         pConn = ADO_temp/2.9968e-12*0.7;
%         % convert pConn into binary connectivity matrix
%         RI = rand(size(pConn));
        CMs = zeros(size(ADO_temp));
%         CMs(RI < pConn) = 1;
%         
%         CMs = conn_reduction(CMs, 700*AR', 90*DR');
    end

    function CMs = simConn_L23toL23_MarSOMtoBipVIP(ADO_temp, AR, DR)
        % for MarSOM-BipVIP connection, within 150 micron p =
        % 0.7
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/3.6168e-12*0.8;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 30*AR', 27*DR');
    end

    function CMs = simConn_L23toL23_MarSOMtoBipCR(ADO_temp, AR, DR)
        % for MarSOM-BipCR connection, 0.7
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/3.6568e-12*0.8;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 24*AR', 30*DR');
    end

    function CMs = simConn_L23toL23_MarSOMtoSbcCR(ADO_temp, AR, DR)
        % for MarSOM-SbcCR connection, within 150 micron p =
        % 0.7
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/3.2968e-12*0.8;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 27*AR', 25*DR');
    end

    function CMs = simConn_L23toL23_MarSOMtoNG(ADO_temp, AR, DR)
        % for MarSOM-NG connection, within 150 micron p =
        % 0.7
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/3.3968e-12*0.7;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 25*AR', 24*DR');
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