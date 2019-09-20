function CM = simConn_L4toL23_SPyrtoAll(ADO_Mat, postIND, postType, AxonRatio, DendRatio)
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
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_SPyrtoPyr(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{2, 3} % post-synaptic cell is FsPV or chandiler
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_SPyrtoFsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(4) % post-synaptic cell is RsPV
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_SPyrtoRsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{5, 6} % post-synaptic cell is MarSOM
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_SPyrtoMarSOM(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case{7, 8} % post-synaptic cell is BipVIP
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_SPyrtoBipVIP(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(9) % post-synaptic cell is BipCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_SPyrtoBipCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(10) % post-synaptic cell in SbcCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_SPyrtoSbcCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
        case(11) % post-synaptic cell in NG cells
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L4toL23_SPyrtoNG(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)));
    end
end

%% nested function
    function CMs = simConn_L4toL23_SPyrtoPyr(ADO_temp, AR, DR)
        % map axon-dendrite overlapping index into connectivity matrix
        % the value is normalized to meausre average connection probability
        % from experiment; for Pyr-Pyr connection, within 50 micron p =
        % 0.19
        % also total number of connection is controlled by the axon and
        % dendrite ratio in the model region
        % assuming a total convengenc/divergence rate of 140 for total
        % connection
        
        pConn = ADO_temp/1.2454e-12*0.22;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
%         CMs = conn_reduction(CMs, 360*AR', 130*DR');
        CMs = conn_reduction(CMs, 380*AR', 160*DR');
    end


    function CMs = simConn_L4toL23_SPyrtoFsPV(ADO_temp, AR, DR)
        % for Pyr-FsPV connection, within 50 micron p =
        % 0.85
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/1.3771e-12*0.4;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 70*AR', 170*DR');
    end


    function CMs = simConn_L4toL23_SPyrtoRsPV(ADO_temp, AR, DR)
        % for Pyr-FsPV connection, within 150 micron p =
        % 0.20
        % assuming a total convengenc rate of 240 and divergence rate of 40
        
        pConn = ADO_temp/1.5930e-12*0.4;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 40*AR', 140*DR');
    end


    function CMs = simConn_L4toL23_SPyrtoMarSOM(ADO_temp, AR, DR)
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

    function CMs = simConn_L4toL23_SPyrtoBipVIP(ADO_temp, AR, DR)
        % for Pyr-FsPV connection, within 150 micron p =
        % 0.21
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/1.5968e-12*0.4;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 40*AR', 140*DR');
    end

    function CMs = simConn_L4toL23_SPyrtoBipCR(ADO_temp, AR, DR)
        % for Pyr-FsPV connection, within 150 micron p =
        % 0.183
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/1.5486e-12*0.4;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 40*AR', 140*DR');
    end

    function CMs = simConn_L4toL23_SPyrtoSbcCR(ADO_temp, AR, DR)
        % for Pyr-FsPV connection, within 150 micron p =
        % 0.174
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/1.3328e-12*0.4;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 30*AR', 150*DR');
    end

    function CMs = simConn_L4toL23_SPyrtoNG(ADO_temp, AR, DR)
        % for Pyr-FsPV connection, within 150 micron p =
        % 0.6
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        pConn = ADO_temp/1.5429e-12*0.3;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMs = zeros(size(pConn));
        CMs(RI < pConn) = 1;
        
        CMs = conn_reduction(CMs, 40*AR', 140*DR');
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