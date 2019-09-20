function CM = simConn_L23toL23_PyrtoAll(ADO_Mat, postIND, postType, AxonRatio, DendRatio, DMs, subtype)
% ADO_Mat: axon-dendritic overlap matrix
% postIND: indexing to seperate different type of post synaptic cells
% postType: numbers specify the type of post-synaptic cell
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
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_PyrtoPyr(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)), subtype(1:size(ADO_Mat,2)));
        case{2, 3} % post-synaptic cell is FsPV or chandiler
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_PyrtoFsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)), subtype(1:size(ADO_Mat,2)));
        case(4) % post-synaptic cell is RsPV
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_PyrtoRsPV(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)), subtype(1:size(ADO_Mat,2)));
        case{5, 6} % post-synaptic cell is MarSOM
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_PyrtoMarSOM(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)), subtype(1:size(ADO_Mat,2)));
        case{7, 8} % post-synaptic cell is BipVIP
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_PyrtoBipVIP(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)), subtype(1:size(ADO_Mat,2)));
        case(9) % post-synaptic cell is BipCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_PyrtoBipCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)), subtype(1:size(ADO_Mat,2)));
        case(10) % post-synaptic cell in SbcCR
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_PyrtoSbcCR(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)), subtype(1:size(ADO_Mat,2)));
        case(11) % post-synaptic cell in NG cells
            CM(postIND(i)+1:postIND(i+1), :) = simConn_L23toL23_PyrtoNG(temp, AxonRatio, DendRatio(postIND(i)+1:postIND(i+1)), subtype(1:size(ADO_Mat,2)));
    end
end

%% nested function
    function CMs = simConn_L23toL23_PyrtoPyr(ADO_temp, AR, DR, st)
        % map axon-dendrite overlapping index into connectivity matrix
        % the value is normalized to meausre average connection probability
        % from experiment; for Pyr-Pyr connection, within 50 micron p =
        % 0.19
        % also total number of connection is controlled by the axon and
        % dendrite ratio in the model region
        % assuming a total convengenc/divergence rate of 140 for total
        % connection
        CMs = zeros(size(ADO_temp));
        
        % need to calculate L2 and L3 pyramidal neurons separatly
        % L2 pyramidal
        idx = find(st == 11);
        pConn = ADO_temp(:, idx)/1.5454e-12*0.24;
        pConn(pConn > 0.18) = 0.18;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 440*AR(idx)', 440*DR');
        CMs(:, idx) = CMtemp;
        
        % L3 pyramidal
        idx = find(st == 12);
        pConn = ADO_temp(:, idx)/4.0454e-12*0.23;
        pConn(pConn > 0.24) = 0.24;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 300*AR(idx)', 300*DR');
        CMs(:, idx) = CMtemp;
        
%         CMs = conn_reduction(CMs, 240*AR', 240*DR');
    end


    function CMs = simConn_L23toL23_PyrtoFsPV(ADO_temp, AR, DR, st)
        % for Pyr-FsPV connection, within 50 micron p =
        % 0.80
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        CMs = zeros(size(ADO_temp));
        
        % need to calculate L2 and L3 pyramidal neurons separatly
        % L2 pyramidal
        idx = find(st == 11);
        pConn = ADO_temp(:, idx)/1.8454e-12*0.80;
        pConn(pConn > 0.96) = 0.96;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 100*AR(idx)', 780*DR');
        CMs(:, idx) = CMtemp;
        
        % L3 pyramidal
        idx = find(st == 12);
        pConn = ADO_temp(:, idx)/3.5454e-12*0.80;
        pConn(pConn > 0.96) = 0.96;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 72*AR(idx)', 720*DR');
        CMs(:, idx) = CMtemp;
    end


    function CMs = simConn_L23toL23_PyrtoRsPV(ADO_temp, AR, DR, st)
        % for Pyr-FsPV connection, within 150 micron p =
        % 0.20
        % assuming a total convengenc rate of 240 and divergence rate of 40
        
        CMs = zeros(size(ADO_temp));
        
        % need to calculate L2 and L3 pyramidal neurons separatly
        % L2 pyramidal
        idx = find(st == 11);
        pConn = ADO_temp(:, idx)/1.2454e-12*0.22;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 64*AR(idx)', 300*DR');
        CMs(:, idx) = CMtemp;
        
        % L3 pyramidal
        idx = find(st == 12);
        pConn = ADO_temp(:, idx)/2.8454e-12*0.20;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 42*AR(idx)', 270*DR');
        CMs(:, idx) = CMtemp;
        
        CMs = conn_reduction(CMs, 40*AR', 240*DR');
    end


    function CMs = simConn_L23toL23_PyrtoMarSOM(ADO_temp, AR, DR, st)
        % for Pyr-MarSOM connection, within 150 micron p =
        % 0.70
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        CMs = zeros(size(ADO_temp));
        
        % need to calculate L2 and L3 pyramidal neurons separatly
        % L2 pyramidal
        idx = find(st == 11);
        pConn = ADO_temp(:, idx)/1.8454e-12*0.7;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 810*AR(idx)', 740*DR');
        CMs(:, idx) = CMtemp;
        
        % L3 pyramidal
        idx = find(st == 12);
        pConn = ADO_temp(:, idx)/3.2454e-12*0.70;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 720*AR(idx)', 700*DR');
        CMs(:, idx) = CMtemp;
    end

    function CMs = simConn_L23toL23_PyrtoBipVIP(ADO_temp, AR, DR, st)
        % for Pyr-FsPV connection, within 150 micron p =
        % 0.21
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        CMs = zeros(size(ADO_temp));
        
        % need to calculate L2 and L3 pyramidal neurons separatly
        % L2 pyramidal
        idx = find(st == 11);
        pConn = ADO_temp(:, idx)/1.6454e-12*0.21;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 36*AR(idx)', 280*DR');
        CMs(:, idx) = CMtemp;
        
        % L3 pyramidal
        idx = find(st == 12);
        pConn = ADO_temp(:, idx)/3.0454e-12*0.21;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 30*AR(idx)', 260*DR');
        CMs(:, idx) = CMtemp;
    end

    function CMs = simConn_L23toL23_PyrtoBipCR(ADO_temp, AR, DR, st)
        % for Pyr-FsPV connection, within 150 micron p =
        % 0.183
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        CMs = zeros(size(ADO_temp));
        
        % need to calculate L2 and L3 pyramidal neurons separatly
        % L2 pyramidal
        idx = find(st == 11);
        pConn = ADO_temp(:, idx)/1.6454e-12*0.19;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 36*AR(idx)', 270*DR');
        CMs(:, idx) = CMtemp;
        
        % L3 pyramidal
        idx = find(st == 12);
        pConn = ADO_temp(:, idx)/3.1454e-12*0.19;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 30*AR(idx)', 250*DR');
        CMs(:, idx) = CMtemp;
    end

    function CMs = simConn_L23toL23_PyrtoSbcCR(ADO_temp, AR, DR, st)
        % for Pyr-FsPV connection, within 150 micron p =
        % 0.174
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        CMs = zeros(size(ADO_temp));
        
        % need to calculate L2 and L3 pyramidal neurons separatly
        % L2 pyramidal
        idx = find(st == 11);
        pConn = ADO_temp(:, idx)/1.5454e-12*0.18;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 36*AR(idx)', 210*DR');
        CMs(:, idx) = CMtemp;
        
        % L3 pyramidal
        idx = find(st == 12);
        pConn = ADO_temp(:, idx)/3.0454e-12*0.18;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 270*AR(idx)', 240*DR');
        CMs(:, idx) = CMtemp;
    end

    function CMs = simConn_L23toL23_PyrtoNG(ADO_temp, AR, DR, st)
        % for Pyr-FsPV connection, within 150 micron p =
        % 0.6
        % assuming a total convengenc rate of 700 and divergence rate of 90
        
        CMs = zeros(size(ADO_temp));
        
        % need to calculate L2 and L3 pyramidal neurons separatly
        % L2 pyramidal
        idx = find(st == 11);
        pConn = ADO_temp(:, idx)/1.2454e-12*0.6;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 41*AR(idx)', 490*DR');
        CMs(:, idx) = CMtemp;
        
        % L3 pyramidal
        idx = find(st == 12);
        pConn = ADO_temp(:, idx)/2.7454e-12*0.6;
        % convert pConn into binary connectivity matrix
        RI = rand(size(pConn));
        CMtemp = zeros(size(pConn));
        CMtemp(RI < pConn) = 1;
        CMtemp = conn_reduction(CMtemp, 30*AR(idx)', 410*DR');
        CMs(:, idx) = CMtemp;
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